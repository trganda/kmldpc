#include "ldpclinearsystem.h"

LDPCLinearSystem::LDPCLinearSystem(toml::value arguments)
    : arguments_(std::move(arguments)),
      codec_(lab::XORSegCodec(arguments_)),
      codec_data_(codec_.uu_len(), codec_.cc_len()),
      modem_linear_system_(lab::ModemLinearSystem(arguments_, codec_.cc_len())) {
    const auto range = toml::find(arguments_, "range");
    min_snr_ = toml::find<double>(range, "minimum_snr");
    max_snr_ = toml::find<double>(range, "maximum_snr");
    step_snr_ = toml::find<double>(range, "step_snr");
    max_err_blk_ = toml::find<int>(range, "maximum_error_number");
    max_num_blk_ = toml::find<int>(range, "maximum_block_number");
    LOG(lab::logger::Info, true) << '[' << std::fixed << std::setprecision(3)
                                 << min_snr_ << ','
                                 << step_snr_ << ','
                                 << max_snr_
                                 << ']' << std::endl;
    LOG(lab::logger::Info, true) << '['
                                 << "MAX_ERROR_BLK = " << max_err_blk_ << ','
                                 << "MAX_BLK = " << max_num_blk_ << ']'
                                 << std::endl;
}

void LDPCLinearSystem::Simulator() {
    // Threads number
    const auto max_threads = (unsigned long) ((max_snr_ - min_snr_) / step_snr_ + 1);
    // Save simulation results
    std::vector<std::pair<double, double>> ber_result(max_threads);
    std::vector<std::pair<double, double>> fer_result(max_threads);
    std::vector<std::future<std::pair<double, double>>> ber_and_fer(max_threads);
    // Record metric for histogram
    const auto histogram = toml::find(arguments_, "histogram");
    const bool histogram_enable = toml::find<bool>(histogram, "enable");
    lab::ThreadsPool threads_pool;
    for (unsigned long i = 0; i < max_threads; i++) {
        ber_and_fer[i] = threads_pool.submit(
            std::bind(
                &LDPCLinearSystem::Run, this, this->codec_, this->modem_linear_system_,
                this->codec_data_, (min_snr_ + step_snr_ * i), histogram_enable
        ));
    }

    for (size_t i=0; i<max_threads; i++) {
        auto result = ber_and_fer[i].get();
        ber_result[i] = std::pair<double, double>((min_snr_ + step_snr_ * i), result.first);
        fer_result[i] = std::pair<double, double>((min_snr_ + step_snr_ * i), result.second);
    }

    LOG(lab::logger::Info, true) << "BER Result" << std::endl;
    for (auto item : ber_result) {
        LOG(lab::logger::Info, true) << std::fixed << std::setprecision(3) << std::setfill('0')
                                     << std::setw(3) << std::right
                                     << item.first << ' '
                                     << std::setprecision(14)
                                     << item.second << std::endl;
    }
    LOG(lab::logger::Info, true) << "FER Result" << std::endl;
    for (auto item : fer_result) {
        LOG(lab::logger::Info, true) << std::fixed << std::setprecision(3) << std::setfill('0')
                                     << std::setw(3) << std::right
                                     << item.first << ' '
                                     << std::setprecision(14)
                                     << item.second << std::endl;
    }
}

std::pair<double, double> LDPCLinearSystem::Run(
    lab::XORSegCodec codec,
    lab::ModemLinearSystem mls,
    const CodecData cdata,
    double snr, bool histogram_enable
) const {
    //var_ = pow(10.0, -0.1 * (snr)) / (codec_.m_coderate * modem_linear_system_.modem_.input_len_);
    double var = pow(10.0, -0.1 * (snr));
    double sigma = sqrt(var);
    mls.set_sigma(sigma);
    mls.set_var(var);
    lab::SourceSink ssink = lab::SourceSink();
    ssink.ClrCnt();
    std::fstream out;
    if (histogram_enable) {
        std::string histfilename = "histogram_" + std::to_string(snr) + ".txt";
        out = std::fstream(histfilename, std::ios::out);
    }
    while ((ssink.num_tot_blk() < max_num_blk_
        && ssink.num_err_blk() < max_err_blk_)) {
        ssink.GetBitStr(cdata.uu_, cdata.uu_len_);
        codec.Encoder(cdata.uu_, cdata.cc_);

        // Generate H
        double real;
        double imag;
        lab::CLCRandNum::Get().Normal(&real, 1);
        lab::CLCRandNum::Get().Normal(&imag, 1);
        std::complex<double> true_h(real, imag);
        true_h *= sqrt(0.5);
        LOG(lab::logger::Info, false) << "Generated H = " << true_h << std::endl;
        std::vector<std::complex<double>> generated_h(1);
        for (auto &i : generated_h) {
            i = true_h;
        }

        // Modulation and pass through the channel
        mls.PartitionModemLSystem(cdata.cc_, generated_h);
        // Get constellation
        auto constellations = mls.constellations();
        // Get received symbols
        auto received_symbols = mls.GetRecvSymbol();
        // KMeans
        kmldpc::KMeans kmeans = kmldpc::KMeans(received_symbols, constellations, 20);
        kmeans.Run();
        auto clusters = kmeans.clusters();
        auto idx = kmeans.idx();

        // Get H hat
        std::complex<double> h_hat = clusters[0] / constellations[0];
        std::vector<std::complex<double>> h_hats(4);
        for (size_t i = 0; i < h_hats.size(); i++) {
            h_hats[i] = h_hat * exp(std::complex<double>(0, (lab::kPi / 2) * i));
        }
        LOG(lab::logger::Info, false) << std::fixed << std::setprecision(0) << std::setfill('0')
                                      << "Current Block Number = "
                                      << std::setw(7) << std::right << (ssink.num_tot_blk() + 1)
                                      << std::endl;
        if (histogram_enable) {
            auto metrics = codec.GetHistogramData(mls, h_hats, cdata.uu_hat_);
            auto idx_of_min = std::distance(
                metrics.begin(),
                min_element(metrics.begin(), metrics.end()));
            for (size_t i = idx_of_min; i < idx_of_min + metrics.size(); i++) {
                out << metrics[i % metrics.size()] << ' ';
            }
            out << std::endl;
        } else {
            codec.Decoder(mls, h_hats, cdata.uu_hat_);
        }
        ssink.CntErr(cdata.uu_, cdata.uu_hat_, codec.uu_len(), 1);
        if (int(ssink.num_tot_blk()) > 0 && int(ssink.num_tot_blk()) % 100 == 0) {
            ssink.PrintResult(snr);
        }
    }
    if (histogram_enable) {
        out.close();
    }
    ssink.PrintResult(snr);
    // BER and FER
    auto ret = std::pair<double, double>(ssink.ber(), ssink.fer());
    return ret;
}
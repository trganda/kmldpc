#include "ldpclinearsystem.h"

LDPCLinearSystem::LDPCLinearSystem(toml::value arguments)
    : arguments_(std::move(arguments)),
      source_sink_(lab::CSourceSink()),
      codec_(lab::XORSegCodec(arguments_)),
      modem_linear_system_(lab::ModemLinearSystem(arguments_, codec_.GetCcLen())) {
  const auto range = toml::find(arguments_, "range");

  min_snr_ = toml::find<double>(range, "minimum_snr");
  max_snr_ = toml::find<double>(range, "maximum_snr");
  step_snr_ = toml::find<double>(range, "step_snr");
  max_err_blk_ = toml::find<int>(range, "maximum_error_number");
  max_num_blk_ = toml::find<int>(range, "maximum_block_number");

  uu_len_ = codec_.GetUuLen();
  cc_len_ = codec_.GetCcLen();

  uu_ = new int[uu_len_];
  uu_hat_ = new int[uu_len_];
  cc_ = new int[cc_len_];
  cc_hat_ = new int[cc_len_];

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

LDPCLinearSystem::~LDPCLinearSystem() {
  delete[] uu_;
  delete[] uu_hat_;
  delete[] cc_;
  delete[] cc_hat_;
}

void LDPCLinearSystem::Simulator() {
  // Save simulation results
  std::vector<std::pair<double, double>> ber_result;
  std::vector<std::pair<double, double>> fer_result;
  // Record metric for histogram
  const auto histogram = toml::find(arguments_, "histogram");
  const bool histogram_enable = toml::find<bool>(histogram, "enable");

  auto hardware_threads = std::thread::hardware_concurrency();

  double snr = min_snr_;
  while (snr <= max_snr_) {
    Run(snr, histogram_enable);
    // BER
    ber_result.emplace_back(snr, source_sink_.GetBer());
    // FER
    fer_result.emplace_back(snr, source_sink_.GetFer());
    snr += step_snr_;
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

void LDPCLinearSystem::Run(double snr, bool histogram_enable) {
  //var_ = pow(10.0, -0.1 * (snr)) / (codec_.m_coderate * modem_linear_system_.modem_.input_len_);
  double var = pow(10.0, -0.1 * (snr));
  double sigma = sqrt(var);

  modem_linear_system_.SetSigma(sigma);
  modem_linear_system_.SetVar(var);
  source_sink_.ClrCnt();

  std::fstream out;
  if (histogram_enable) {
    std::string histfilename = "histogram_" + std::to_string(snr) + ".txt";
    out = std::fstream(histfilename, std::ios::out);
  }
  while ((source_sink_.GetNumTotBlk() < max_num_blk_
      && source_sink_.GetNumErrBlk() < max_err_blk_)) {

    source_sink_.GetBitStr(uu_, codec_.GetUuLen());
    codec_.Encoder(uu_, cc_);

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
    modem_linear_system_.PartitionModemLSystem(cc_, generated_h);
    // Get constellation
    auto constellations = modem_linear_system_.GetConstellations();
    // Get received symbols
    auto received_symbols = modem_linear_system_.GetRecvSymbol();
    // KMeans
    kmldpc::KMeans kmeans = kmldpc::KMeans(received_symbols, constellations, 20);
    kmeans.Run();
    auto clusters = kmeans.GetClusters();
    auto idx = kmeans.GetIdx();

    // Get H hat
    std::complex<double> h_hat = clusters[0] / constellations[0];
    std::vector<std::complex<double>> h_hats(4);
    for (size_t i = 0; i < h_hats.size(); i++) {
      h_hats[i] = h_hat * exp(std::complex<double>(0, (lab::kPi / 2) * i));
    }

    LOG(lab::logger::Info, false) << std::fixed << std::setprecision(0) << std::setfill('0')
                                  << "Current Block Number = "
                                  << std::setw(7) << std::right << (source_sink_.GetNumTotBlk() + 1)
                                  << std::endl;

    if (histogram_enable) {
      auto metrics = codec_.GetHistogramData(modem_linear_system_, h_hats, uu_hat_);
      auto idx_of_min = std::distance(metrics.begin(),
                                      min_element(metrics.begin(), metrics.end()));
      for (size_t i = idx_of_min; i < idx_of_min + metrics.size(); i++) {
        out << metrics[i % metrics.size()] << ' ';
      }
      out << std::endl;
    } else {
      codec_.Decoder(modem_linear_system_, h_hats, uu_hat_);
    }

    source_sink_.CntErr(uu_, uu_hat_, codec_.GetUuLen(), 1);

    if (int(source_sink_.GetNumTotBlk()) > 0 && int(source_sink_.GetNumTotBlk()) % 100 == 0) {
      source_sink_.PrintResult(snr);
    }
  }
  if (histogram_enable) {
    out.close();
  }
  source_sink_.PrintResult(snr);
}
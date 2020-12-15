#include "ldpclinearsystem.h"

LDPCLinearSystem::LDPCLinearSystem()
    : source_sink_(lab::CSourceSink()),
      codec_(lab::XORSegCodec()), modem_linear_system_(lab::ModemLinearSystem()),
      min_snr_(0.0), max_snr_(0.0), step_snr_(0.0),
      max_err_blk_(0), max_num_blk_(0),
      uu_(nullptr), uu_hat_(nullptr), uu_len_(0),
      cc_(nullptr), cc_hat_(nullptr), cc_len_(0), sym_prob_(nullptr) {}

LDPCLinearSystem::~LDPCLinearSystem() {
    delete[] uu_;
    delete[] uu_hat_;
    delete[] cc_;
    delete[] cc_hat_;
    delete[] sym_prob_;
}

void LDPCLinearSystem::InitSimulator()
{
	char codec_file[80];
	char modem_file[80];
	char temp_str[80];
	FILE* fp;

	std::string config_file = "Setup_of_LDPC_Linear_System0.txt";
	if ((fp = fopen(config_file.c_str(), "r")) == nullptr)
	{
		fprintf(stderr, "\nCan't Open the %s file!\n", config_file.c_str());
		exit(-1);
	}

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%lf", &min_snr_);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%lf", &max_snr_);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%lf", &step_snr_);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &max_err_blk_);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &max_num_blk_);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", codec_file);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", modem_file);
	fclose(fp);

	codec_.Malloc(0, codec_file);

    uu_len_ = codec_.GetUuLen();
    cc_len_ = codec_.GetCcLen();

	modem_linear_system_.Malloc(cc_len_, 0, modem_file);

    uu_ = new int[uu_len_];
    uu_hat_ = new int[uu_len_];

    cc_ = new int[cc_len_];
    cc_hat_ = new int[cc_len_];

    sym_prob_ = new double[cc_len_ / modem_linear_system_.GetModem().GetInputLen()
                           * modem_linear_system_.GetModem().GetNumSymbol()];

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

void LDPCLinearSystem::Simulator()
{
	double sigma, var;

    InitSimulator();

    // Save simulation results
	std::vector<std::pair<double, double>> ber_result;
	std::vector<std::pair<double, double>> fer_result;

	double snr = min_snr_;
	while (snr <= max_snr_) {
        //var_ = pow(10.0, -0.1 * (snr)) / (codec_.m_coderate * modem_linear_system_.modem_.input_len_);
        var = pow(10.0, -0.1 * (snr));

        sigma = sqrt(var);

        modem_linear_system_.GetLinearSystem().SetSigma(sigma);
        modem_linear_system_.GetLinearSystem().SetVar(var);

        source_sink_.ClrCnt();

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
            for (auto & i : generated_h) {
                i = true_h;
            }

            // Modulation and pass through the channel
            modem_linear_system_.MLSystemPartition(cc_, generated_h);
            // Get constellation
            auto constellations = modem_linear_system_.GetLinearSystem().GetMModem()->GetConstellations();
            // Get received symbols
            auto received_symbols = modem_linear_system_.GetRSymbol();
            // KMeans
            kmldpc::KMeans kmeans = kmldpc::KMeans(received_symbols, constellations, 20);
            kmeans.Run();
            auto clusters = kmeans.GetClusters();
            auto idx = kmeans.GetIdx();

            // Get H hat
            // std::complex<double> h_hat = true_h;
            std::complex<double> h_hat = clusters[0] / constellations[0];
            std::vector<std::complex<double>> h_hats(4);
            for (size_t i = 0; i < h_hats.size(); i++) {
                h_hats[i] = h_hat * exp(std::complex<double>(0, (lab::kPi / 2) * i));
            }

            LOG(lab::logger::Info, false) << std::fixed << std::setprecision(0) << std::setfill('0')
                                          << "Current Block Number = "
                                          << std::setw(7) << std::right << (source_sink_.GetNumTotBlk() + 1)
                                          << std::endl;

            codec_.Decoder(modem_linear_system_, h_hats, uu_hat_);

            source_sink_.CntErr(uu_, uu_hat_, codec_.GetUuLen(), 1);

            if (int(source_sink_.GetNumTotBlk()) > 0 && int(source_sink_.GetNumTotBlk()) % 100 == 0) {
                source_sink_.PrintResult(snr);
            }
        }
        source_sink_.PrintResult(snr);

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

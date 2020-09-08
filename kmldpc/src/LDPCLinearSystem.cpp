#include "LDPCLinearSystem.h"
#include "Log.h"

extern CLCRandNum rndGen0;
extern CWHRandNum rndGen1;


void LDPC_Linear_System::StartSimulator()
{
	int setup_no;
	char file_name[80];
	char codec_file[80];
	char modem_file[80];
	char temp_str[80];
	FILE* fp;

	setup_no = 0;
	sprintf(file_name, "Setup_of_LDPC_Linear_System%d.txt", setup_no);

	if ((fp = fopen(file_name, "r")) == nullptr)
	{
		fprintf(stderr, "\nCan't open the %s file!\n", file_name);
		exit(1);
	}

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%lf", &m_min_snr);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%lf", &m_max_snr);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%lf", &m_inc_snr);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_max_blk_err);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_max_blk_num);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", codec_file);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", modem_file);
	fclose(fp);

	if ((fp = fopen("snrresult.txt", "a+")) == nullptr)
	{
		fprintf(stderr, "\n Cannot open the file!!!\n");
		exit(1);
	}

	fprintf(fp, "\n\n*** *** *** *** *** *** ***");
	fprintf(fp, "\nStart Time: ");
	loctime = time(nullptr);
	ptr = localtime(&loctime);
	fprintf(fp, asctime(ptr));
	fprintf(fp, "\n The results correspond to %s", file_name);
	fclose(fp);

	m_codec.Malloc(0, codec_file);

	m_total_angle = 2;
	m_len_uu = m_codec.m_extrabits_len + m_codec.m_len_uu;
	m_len_cc = m_codec.m_len_cc;

	Modem_Lin_Sym.Malloc(m_len_cc, 0, modem_file);

	m_uu = new int[m_len_uu];
	m_uu_hat = new int[m_len_uu];

	m_cc = new int[m_len_cc];
	m_cc_hat = new int[m_len_cc];

	m_sym_prob = new double[m_len_cc / Modem_Lin_Sym.m_modem.len_input * Modem_Lin_Sym.m_modem.num_symbol];

    LOG(kmldpc::Info, true) << '[' << std::fixed << std::setprecision(3)
                           << m_min_snr << ','
                           << m_inc_snr << ','
                           << m_max_snr
                           << ']' << std::endl;
    LOG(kmldpc::Info, true) << '['
                           <<"MAX_ERROR_BLK = " << m_max_blk_err << ','
                           << "MAX_BLK = " << m_max_blk_num << ']'
                           << std::endl;

	//输出仿真参数到文件
	if ((fp = fopen("snrresult.txt", "a+")) == nullptr) {
		fprintf(stderr, "\n Cannot open the file!!!\n");
		exit(1);
	}
	fprintf(fp, "\n SNR: [%2.2lf : %2.2lf : %2.2lf]", m_min_snr, m_inc_snr, m_max_snr);
	fprintf(fp, "\n max_blk_err : %d", m_max_blk_err);
	fprintf(fp, "\n max_blk_num : %d", m_max_blk_num);
	Modem_Lin_Sym.m_modem.PrintCodeParameter(fp);
	fclose(fp);

	if ((fp = fopen("snrber.txt", "a+")) == nullptr) {
		fprintf(stderr, "\n Cannot open the file!!!\n");
		exit(1);
	}
	fprintf(fp, "\n SNR: [%2.2lf : %2.2lf : %2.2lf]", m_min_snr, m_inc_snr, m_max_snr);
	fprintf(fp, "\n max_blk_err : %d", m_max_blk_err);
	fprintf(fp, "\n max_blk_num : %d", m_max_blk_num);
	Modem_Lin_Sym.m_modem.PrintCodeParameter(fp);
	fclose(fp);

	if ((fp = fopen("snrfer.txt", "a+")) == nullptr) {
		fprintf(stderr, "\n Cannot open the file!!!\n");
		exit(1);
	}
	fprintf(fp, "\n SNR: [%2.2lf : %2.2lf : %2.2lf]", m_min_snr, m_inc_snr, m_max_snr);
	fprintf(fp, "\n max_blk_err : %d", m_max_blk_err);
	fprintf(fp, "\n max_blk_num : %d", m_max_blk_num);
	Modem_Lin_Sym.m_modem.PrintCodeParameter(fp);
	fclose(fp);
}

void LDPC_Linear_System::EndSimulator()
{
	delete[]m_uu;
	delete[]m_uu_hat;

	delete[]m_cc;
	delete[]m_cc_hat;

	delete[]m_sym_prob;

	m_codec.Free();
	Modem_Lin_Sym.Free();
}

void LDPC_Linear_System::Simulator()
{
	double snr, sigma, var;
	FILE* fp = nullptr;

	StartSimulator();

	std::vector<std::complex<double>> totalH;
	std::vector<std::complex<double>> um;

	for (int i = 0; i < m_total_angle; i++) {
		std::complex<double> H(0.0, ((m_PI / 2) / (m_total_angle)) * (i));
		H = exp(H);
		totalH.push_back(H);
	}

	for (snr = m_min_snr; snr < m_max_snr + 0.5 * m_inc_snr; snr += m_inc_snr)
	{
		//var = pow(10.0, -0.1 * (snr)) / (m_codec.m_coderate * Modem_Lin_Sym.m_modem.len_input);
		var = pow(10.0, -0.1 * (snr));

		sigma = sqrt(var);

		Modem_Lin_Sym.Lin_Sym.sigma = sigma;
		Modem_Lin_Sym.Lin_Sym.var = var; 

		LOG(kmldpc::Info, true) << std::fixed << std::setprecision(2) << "SNR = " << snr << std::endl;

		m_source_sink.ClrCnt();
		m_source_extrabits.ClrCnt();

		while ((m_source_sink.m_num_tot_blk < m_max_blk_num 
			&& m_source_sink.m_num_err_blk < m_max_blk_err)) {

			int current_error_blk = m_source_sink.m_num_err_blk;

			m_source_sink.GetBitStr(m_uu, m_codec.m_len_uu);
			m_source_extrabits.GetBitStr(m_uu + m_codec.m_len_uu, m_codec.m_extrabits_len);
			m_codec.Encoder(m_uu, m_cc);

			// Generate H
			double real;
			double imag;
			rndGen0.Normal(&real, 1);
			rndGen0.Normal(&imag, 1);

			std::complex<double> trueH(real, imag);
			trueH *= sqrt(0.5);
			LOG(kmldpc::Info, false) << "Generated H = " << trueH << std::endl;
			std::vector<std::complex<double>> selecth(1);
			for (auto & i : selecth) {
				i = trueH;
			}

			// Modulation and pass through the channel
			Modem_Lin_Sym.modem_linear_system_parition(m_cc, selecth);
			// Get constellation
			auto constellations = Modem_Lin_Sym.Lin_Sym.m_modem->getConstellation();
			// Get received symbols
			auto receivedSymbols = Modem_Lin_Sym.getRSymbol();
			// KMeans
			kmldpc::KMeans kmeans = kmldpc::KMeans(receivedSymbols, constellations, 20);
			kmeans.run();
			auto clusters = kmeans.getClusters();
			auto idx = kmeans.getIdx();

			// Get H hat
			// std::complex<double> hHat = trueH;
			std::complex<double> hHat = clusters[0] / constellations[0];
			std::vector<std::complex<double>> hHats(4);
			for (int i = 0; i < hHats.size(); i++) {
				hHats[i] = hHat * exp(std::complex<double>(0, (m_PI / 2) * i));
			}

			m_codec.Decoder(Modem_Lin_Sym, hHats, m_uu_hat);

			m_source_sink.CntErr(m_uu, m_uu_hat, m_codec.m_len_uu, 1);

			if (m_source_sink.m_num_err_blk > current_error_blk) {
			    std::string dir = "records";
				std::string filename = dir + "/" + "RECEIVED_SYMBOL_SNR_" +
					std::to_string(snr) + "_ID_" +
					std::to_string(int(m_source_sink.m_num_tot_blk)) + ".mat";
				hHats.push_back(trueH);

                kmeans.dumpToMat(filename, hHats);
                LOG(kmldpc::Info, false) << "Wrote error case to " << filename << std::endl;
			}

			if (int(m_source_sink.m_num_tot_blk) > 0 && int(m_source_sink.m_num_tot_blk) % 100 == 0) {
                m_source_sink.PrintResult();
			}
		}

		m_source_sink.PrintResult();

		if ((fp = fopen("snrresult.txt", "a+")) == nullptr)
		{
			fprintf(stderr, "\n Cannot open the file!!!\n");
			exit(1);
		}
		fprintf(fp, "\n\nsnr =  %lf: var = %lf: ", snr, var);
		m_source_sink.PrintResult(fp);
		fclose(fp);

		//BER
		if ((fp = fopen("snrber.txt", "a+")) == nullptr)
		{
			fprintf(stderr, "\n Cannot open the file!!!\n");
			exit(1);
		}
		fprintf(fp, "\n%lf %12.10lf", snr, m_source_sink.m_ber);
		fclose(fp);

		//FER
		if ((fp = fopen("snrfer.txt", "a+")) == nullptr)
		{
			fprintf(stderr, "\n Cannot open the file!!!\n");
			exit(1);
		}
		fprintf(fp, "\n%lf %12.10lf", snr, m_source_sink.m_fer);
		fclose(fp);

		if ((fp = fopen("snrresult.txt", "a+")) == nullptr)
		{
			fprintf(stderr, "\n Cannot open the file!!!\n");
			exit(1);
		}
		fclose(fp);

	}

	EndSimulator();
}
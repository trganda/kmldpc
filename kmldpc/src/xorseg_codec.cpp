#include "xorseg_codec.h"
#include "randnum.h"
#include <vector>
#include <complex>
#include <algorithm>

extern CLCRandNum rndGen0;

void XORSegCodec::Malloc(int code_no, const char *file_name)
{
	// setup from "Setup_of_Codec.txt"
	char temp_str[80] = {' '};
	char mark[80];

	FILE *fp = fopen(file_name, "r");
	if (nullptr == fp)
	{
		fprintf(stderr, "\nCannot open %s", file_name);
		exit(-1);
	}

	sprintf(mark, "RanXORLDPC***%d***PARAMETERS", code_no);
	while (strcmp(temp_str, mark))
	{
		fscanf(fp, "%s", temp_str);
	}

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_extrabits_len);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_list_count);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", temp_str);

	fclose(fp);

	// setup from "LDPC.txt"
	m_LDPC_codec.Malloc(code_no, temp_str);
	m_len_uu = m_LDPC_codec.m_codedim;
	m_len_cc = m_LDPC_codec.m_codelen;
	m_coderate = double(m_len_uu) / m_len_cc;

	// random generator matrix
	m_generator = new int *[m_len_cc];
	for (int i = 0; i < m_extrabits_len; i++)
	{
		m_generator[i] = new int[m_len_cc];
		for (int j = 0; j < m_len_cc; j++)
		{
			m_generator[i][j] = (rndGen0.Uniform() < 0.5 ? 0 : 1);
		}
	}

	m_uu_Ran = new int[m_extrabits_len];
	m_cc_Ran = new int[m_len_cc];
	m_rr = new int[m_len_cc];
	bitLin = new double[m_len_cc];
	bitLout = new double[m_len_cc];
	//m_p0_cc = new double[m_len_cc];
}

void XORSegCodec::Free()
{
	delete[] m_uu_Ran;
	delete[] m_cc_Ran;
	delete[] m_rr;
	delete[] bitLin;
	delete[] bitLout;
	for (int i = 0; i < m_extrabits_len; i++)
	{
		delete[] m_generator[i];
	}
	delete[] m_generator;
	//delete[] m_p0_cc;
}

void XORSegCodec::Encoder(int *uu, int *cc)
{
	m_LDPC_codec.Encoder(uu, cc);
	//MatrixProd(uu + m_len_uu, m_cc_Ran, m_generator, m_extrabits_len, m_len_cc);
	//for (int i = 0; i < m_len_cc; i++) {
	//	cc[i] ^= m_cc_Ran[i];
	//}
}

void XORSegCodec::Encoder_Ran(int *uu, int *cc)
{
	MatrixProd(uu, cc, m_generator, m_extrabits_len, m_len_cc);
}

void XORSegCodec::Decoder(Modem_Linear_System &modem_linear_system, std::vector<std::complex<double>> &totalH, int *uu_hat, std::vector<int> &difflisterr)
{
	Decoder_Ran(modem_linear_system, totalH, uu_hat, difflisterr);

	//std::vector<std::pair<int, std::complex<double>>> thetaList(m_extrabits_len);
	//for (int i = 0; i < thetaList.size(); i++) {
	//	thetaList[i] = std::pair<int, std::complex<double>>(uu_hat[m_len_uu + i], totalH[uu_hat[m_len_uu + i]]);
	//}
	//Demmaping(modem_linear_system, thetaList);
}

void XORSegCodec::Decoder(Modem_Linear_System &modem_linear_system, std::vector<std::complex<double>> &hHats, int *uu_hat)
{
	std::vector<int> parityResults(hHats.size());
    std::vector<std::pair<int, std::complex<double>>> temp;
	for (int i = 0; i < parityResults.size(); i++)
	{
		temp = {std::pair<int, std::complex<double>>(0, hHats[i])};
		Demmaping(modem_linear_system, temp);
		parityResults[i] = getParityCheck();
		LOG(kmldpc::Info) << "Hhat = " << hHats[i] << "Parity Count = " << parityResults[i] << std::endl;
	}

	auto minIndex = std::distance(parityResults.begin(), min_element(parityResults.begin(), parityResults.end()));
    LOG(kmldpc::Info) << "minIndex = " << minIndex << std::endl;
	temp = {std::pair<int, std::complex<double>>(0, hHats[minIndex])};
	Demmaping(modem_linear_system, temp);

	m_LDPC_codec.Decoder(bitLout, uu_hat);
}

void XORSegCodec::Decoder_Ran(Modem_Linear_System &modem_linear_system, std::vector<std::complex<double>> &totalH, int *uu_hat, std::vector<int> &difflisterr)
{
	std::vector<std::pair<int, std::complex<double>>> thetaList(m_extrabits_len);
	// get gray modulation constellation
	std::vector<std::complex<double>> graySymbol(modem_linear_system.m_modem.num_symbol);
	for (int i = 0; i < graySymbol.size(); i++)
	{
		graySymbol[i] = std::complex<double>(modem_linear_system.m_modem.output_symbol[i][0],
											 modem_linear_system.m_modem.output_symbol[i][1]);
	}
	// get received symbol
	std::vector<std::complex<double>> um;
	um = modem_linear_system.getRSymbol();
	// calculate probability
	int totalangle = totalH.size();
	std::vector<std::vector<double>> optsets(thetaList.size(), std::vector<double>(totalH.size()));
	int num_of_symbol = modem_linear_system.m_len_xx / 2;
	int symbolPerPart = num_of_symbol / m_extrabits_len;
	for (int i = 0; i < thetaList.size(); i++)
	{
		for (int j = 0; j < optsets[i].size(); j++)
		{
			std::vector<std::complex<double>> curdata(um.begin() + i * symbolPerPart, um.begin() + (i + 1) * symbolPerPart);
			optsets[i][j] = adaptFunc(curdata,
									  graySymbol, totalH[j], modem_linear_system.Lin_Sym.var);
		}
		//convertToPb(optsets[i]);
		auto maxidx = std::distance(optsets[i].begin(), std::max_element(optsets[i].begin(), optsets[i].end()));
		thetaList[i] = std::pair<int, std::complex<double>>(maxidx, totalH[maxidx]);
	}

	std::vector<std::pair<int, double>> differ(optsets.size());
	for (int i = 0; i < differ.size(); i++)
	{
		differ[i] = std::pair<int, double>(i, std::abs(optsets[i][0] - optsets[i][1]));
	}
	auto cmp = [](std::pair<int, double> &a, std::pair<int, double> &b) { return a.second < b.second; };
	std::sort(differ.begin(), differ.end(), cmp);

	int maxlistsize = 1 << m_list_count;

	for (int i = 0; i < difflisterr.size(); i++)
	{
		auto tempthetalist = thetaList;
		if (!isCorrectInList(tempthetalist, differ, totalH, uu_hat, i))
		{
			difflisterr[i]++;
		}
		else
		{
			break;
		}
	}
}

int XORSegCodec::Decoder_Ran_Hist(Modem_Linear_System &modem_linear_system, std::vector<std::complex<double>> &totalH, int *uu)
{
	std::vector<std::pair<int, std::complex<double>>> thetaList(m_extrabits_len);
	// get gray modulation constellation
	std::vector<std::complex<double>> graySymbol(modem_linear_system.m_modem.num_symbol);
	for (int i = 0; i < graySymbol.size(); i++)
	{
		graySymbol[i] = std::complex<double>(modem_linear_system.m_modem.output_symbol[i][0],
											 modem_linear_system.m_modem.output_symbol[i][1]);
	}
	int totalangle = totalH.size();
	// get received symbol
	std::vector<std::complex<double>> um;
	um = modem_linear_system.getRSymbol();

	for (int i = 0; i < thetaList.size(); i++)
	{
		thetaList[i] = std::pair<int, std::complex<double>>(i, totalH[i]);
	}

	Demmaping(modem_linear_system, thetaList);
	int result = getParityCheck(thetaList, totalangle);
	return result;
}

void XORSegCodec::Demmaping(Modem_Linear_System &modem_linear_system, std::vector<std::pair<int, std::complex<double>>> &thetaList)
{
	modem_linear_system.Soft_Demodulation(thetaList);

	// demapping to get soft information
	for (int i = 0; i < m_len_cc; i++)
	{
		bitLin[i] = 0.5;
	}

	modem_linear_system.m_modem.Demapping(bitLin, modem_linear_system.m_sym_prob, bitLout, m_len_cc);
}

double XORSegCodec::adaptFunc(std::vector<std::complex<double>> &data, std::vector<std::complex<double>> &graySymbol, std::complex<double> &h, double var)
{
	double result = 0.0;
	double co = 1.0 / sqrt(2 * m_PI * var);

	auto tempSymbol = graySymbol;
	for (auto &i : tempSymbol)
	{
		i *= h;
	}

	for (auto i : data)
	{
		double logsum = 0.0;
		for (auto j : tempSymbol)
		{
			double alpha = 1.0 / tempSymbol.size();
			std::complex<double> differ = i - j;
			logsum += alpha * co * exp(-((pow(differ.real(), 2) + pow(differ.imag(), 2)) / var));
		}
		result += log(logsum);
	}
	result /= double(data.size());
	return result;
}

void XORSegCodec::convertToPb(std::vector<double> &opts)
{
	double sum = 0.0;
	auto maxidx = max_element(opts.begin(), opts.end());
	for (auto &j : opts)
	{
		j = exp(j - *maxidx);
		sum += j;
	}
	for (auto &j : opts)
	{
		j /= sum;
	}
}

int XORSegCodec::getParityCheck(std::vector<std::pair<int, std::complex<double>>> &thetaList, int totalangle)
{
	// hard decision
	for (int i = 0; i < m_len_cc; i++)
	{
		if (bitLout[i] > 0.5)
		{
			m_rr[i] = 1;
		}
		else
		{
			m_rr[i] = 0;
		}
	}
	// xor
	for (int i = 0; i < m_extrabits_len; i++)
	{
		int bitanglelen = totalangle >> 1;
		for (int j = 0; j < bitanglelen; j++)
		{
			m_uu_Ran[i * bitanglelen + j] = (thetaList[i].first >> j) % 2;
		}
	}
	MatrixProd(m_uu_Ran, m_cc_Ran, m_generator, m_extrabits_len, m_len_cc);
	for (int i = 0; i < m_len_cc; i++)
	{
		m_rr[i] ^= m_cc_Ran[i];
	}

	return m_LDPC_codec.parityCheck(m_rr);
}

int XORSegCodec::getParityCheck()
{
	// hard decision
	for (int i = 0; i < m_len_cc; i++)
	{
		if (bitLout[i] > 0.5)
		{
			m_rr[i] = 1;
		}
		else
		{
			m_rr[i] = 0;
		}
	}

	return m_LDPC_codec.parityCheck(m_rr);
}

bool XORSegCodec::isCorrectInList(std::vector<std::pair<int, std::complex<double>>> &thetaList,
								  std::vector<std::pair<int, double>> &differ,
								  std::vector<std::complex<double>> &totalH,
								  int *uu, int list_len)
{
	int listsize = 1 << list_len;
	if (1 == listsize)
	{
		int temperr = 0;
		for (int i = 0; i < m_extrabits_len; i++)
		{
			if (thetaList[i].first != uu[m_len_uu + i])
			{
				temperr++;
			}
		}
		if (temperr > 0)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	std::vector<int> debug(list_len);
	for (int i = 0; i < listsize; i++)
	{
		for (int j = 0; j < list_len; j++)
		{
			debug[j] = (i >> j) % 2;
		}
		for (int j = 0; j < list_len; j++)
		{
			thetaList[differ[j].first] = std::pair<int, std::complex<double>>(debug[j], totalH[debug[j]]);
		}
		int temperr = 0;
		for (int j = 0; j < m_extrabits_len; j++)
		{
			if (thetaList[j].first != uu[m_len_uu + j])
			{
				temperr++;
			}
		}
		if (temperr > 0)
		{
			continue;
		}
		else
		{
			return true;
		}
	}

	return false;
}

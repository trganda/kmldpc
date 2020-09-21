﻿#include "XorSegCodec.h"
#include "RandNum.h"
#include <vector>
#include <complex>
#include <algorithm>
#include <iomanip>

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
    fscanf(fp, "%u", &m_using_5G_LDPC);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_iter_cnt);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%u", &m_using_Syndrom_Metric);

    fscanf(fp, "%s", temp_str);
    fscanf(fp, "%u", &m_histogram);

    fscanf(fp, "%s", temp_str);
    fscanf(fp, "%s", temp_str);

	fclose(fp);

	// setup from "LDPC.txt"
	if (m_using_5G_LDPC == 1) {
	    LOG(kmldpc::Info, true) << "Using 5G LDPC." << std::endl;
        m_5GLDPC_codec.Malloc(code_no, temp_str);
        m_len_uu = m_5GLDPC_codec.m_codedim;
        m_len_cc = m_5GLDPC_codec.m_codelen_puncture;
    } else {
        LOG(kmldpc::Info, true) << "Using traditional LDPC." << std::endl;
	    m_LDPC_codec.Malloc(code_no, temp_str);
	    m_len_uu = m_LDPC_codec.m_codedim;
	    m_len_cc = m_LDPC_codec.m_codelen;
	}

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
    if (m_using_5G_LDPC == 1) {
        m_5GLDPC_codec.Encoder_5G(uu, cc);
    } else {
        m_LDPC_codec.Encoder(uu, cc);
    }
	//MatrixProd(uu + m_len_uu, m_cc_Ran, m_generator, m_extrabits_len, m_len_cc);
	//for (int i = 0; i < m_len_cc; i++) {
	//	cc[i] ^= m_cc_Ran[i];
	//}
}

void XORSegCodec::Encoder_Ran(int *uu, int *cc)
{
	MatrixProd(uu, cc, m_generator, m_extrabits_len, m_len_cc);
}

void XORSegCodec::Decoder(Modem_Linear_System &modem_linear_system,
                          std::vector<std::complex<double>> &hHats, int *uu_hat)
{
	std::vector<double> metricResults(hHats.size(), 0);
	std::vector<std::vector<double>> softSyndromsData(hHats.size());
    std::vector<std::pair<int, std::complex<double>>> temp;
	for (auto i = 0; i < metricResults.size(); i++)
	{
		temp = {std::pair<int, std::complex<double>>(0, hHats[i])};
		Demmaping(modem_linear_system, temp);
		if (m_using_5G_LDPC == 1) {
            m_5GLDPC_codec.Decoder_5G(bitLout, uu_hat, m_iter_cnt);
		    if (m_using_Syndrom_Metric == 1) {
		        softSyndromsData[i] = std::vector<double> (m_5GLDPC_codec.m_syndromsoft,
                                                           m_5GLDPC_codec.m_syndromsoft+m_5GLDPC_codec.m_num_row);
                for (auto j = 0; j < m_5GLDPC_codec.m_num_row; j++) {
                    metricResults[i] += log(softSyndromsData[i][j]);
                }
		    } else {
                metricResults[i] = getParityCheckAfterDecoding();
		    }
		} else {
		    if (m_using_Syndrom_Metric == 1) {
                m_LDPC_codec.Decoder(bitLout, uu_hat, m_iter_cnt);
                softSyndromsData[i] = std::vector<double> (m_LDPC_codec.m_syndromsoft,
                                                           m_LDPC_codec.m_syndromsoft+m_LDPC_codec.m_codedim);
                for (auto j = 0; j < m_LDPC_codec.m_codedim; j++) {
                    metricResults[i] += log(softSyndromsData[i][j]);
                }
		    } else {
                metricResults[i] = getParityCheck();
		    }
		}

		LOG(kmldpc::Info, false) << std::fixed << std::setprecision(14)
                                 << "Hhat = " << hHats[i]
                                 << " Metric = "
                                 << std::setw(5) << std::right
                                 << metricResults[i] << std::endl;
		metricResults[i] = abs(metricResults[i]);
	}

	auto minIndex = std::distance(metricResults.begin(), min_element(metricResults.begin(), metricResults.end()));
    LOG(kmldpc::Info, false) << "hatIndex = " << minIndex << std::endl;
	temp = {std::pair<int, std::complex<double>>(0, hHats[minIndex])};
	Demmaping(modem_linear_system, temp);

	if (m_using_5G_LDPC == 1) {
	    m_5GLDPC_codec.Decoder_5G(bitLout, uu_hat, m_5GLDPC_codec.m_max_iter);
	} else {
	    m_LDPC_codec.Decoder(bitLout, uu_hat, m_LDPC_codec.m_max_iter);
	}
}

double XORSegCodec::Histogram(Modem_Linear_System &modem_linear_system, std::vector<std::complex<double>> &hHats,
                              int *uu_hat) {
    std::vector<double> metricResults(hHats.size(), 0);
    std::vector<std::vector<double>> softSyndromsData(hHats.size());
    std::vector<std::pair<int, std::complex<double>>> temp;
    for (auto i = 0; i < metricResults.size(); i++)
    {
        temp = {std::pair<int, std::complex<double>>(0, hHats[i])};
        Demmaping(modem_linear_system, temp);
        if (m_using_5G_LDPC == 1) {
            m_5GLDPC_codec.Decoder_5G(bitLout, uu_hat, m_iter_cnt);
            if (m_using_Syndrom_Metric == 1) {
                softSyndromsData[i] = std::vector<double> (m_5GLDPC_codec.m_syndromsoft,
                                                           m_5GLDPC_codec.m_syndromsoft+m_5GLDPC_codec.m_num_row);
                for (auto j = 0; j < m_5GLDPC_codec.m_num_row; j++) {
                    metricResults[i] += log(softSyndromsData[i][j]);
                }
            } else {
                metricResults[i] = getParityCheckAfterDecoding();
            }
        } else {
            if (m_using_Syndrom_Metric == 1) {
                m_LDPC_codec.Decoder(bitLout, uu_hat, m_iter_cnt);
                softSyndromsData[i] = std::vector<double> (m_LDPC_codec.m_syndromsoft,
                                                           m_LDPC_codec.m_syndromsoft+m_LDPC_codec.m_codedim);
                for (auto j = 0; j < m_LDPC_codec.m_codedim; j++) {
                    metricResults[i] += log(softSyndromsData[i][j]);
                }
            } else {
                metricResults[i] = getParityCheck();
            }
        }

        LOG(kmldpc::Info, false) << std::fixed << std::setprecision(14)
                                 << "Hhat = " << hHats[i]
                                 << " Metric = "
                                 << std::setw(5) << std::right
                                 << metricResults[i] << std::endl;
        metricResults[i] = abs(metricResults[i]);
    }

    return metricResults[0];
}

void XORSegCodec::Demmaping(Modem_Linear_System &modem_linear_system,
                            std::vector<std::pair<int, std::complex<double>>> &thetaList)
{
	modem_linear_system.Soft_Demodulation(thetaList);

	// demapping to get soft information
	for (int i = 0; i < m_len_cc; i++)
	{
		bitLin[i] = 0.5;
	}

	modem_linear_system.m_modem.Demapping(bitLin, modem_linear_system.m_sym_prob, bitLout, m_len_cc);
}

double XORSegCodec::adaptFunc(std::vector<std::complex<double>> &data,
                              std::vector<std::complex<double>> &graySymbol,
                              std::complex<double> &h, double var)
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

int XORSegCodec::getParityCheckAfterDecoding() {
    return m_5GLDPC_codec.parityCheck(m_5GLDPC_codec.m_cc_hat);
}


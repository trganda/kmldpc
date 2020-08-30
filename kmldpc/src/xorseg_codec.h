#ifndef XOR_SEG_CODEC_H
#define XOR_SEG_CODEC_H

#include "binary_ldpc_codec.h"
#include "modem_linear_system.h"

class XORSegCodec
{
public:
	XORSegCodec(){};
	~XORSegCodec(){};

public:
	void Malloc(int code_no, const char *file_name);
	void Free();

public:
	void Encoder(int *uu, int *cc);
	void Encoder_Ran(int *uu, int *cc);
	void Decoder(Modem_Linear_System &modem_linear_system, std::vector<std::complex<double>> &totalH, int *uu_hat, std::vector<int> &difflisterr);
	void Decoder(Modem_Linear_System &modem_linear_system, std::vector<std::complex<double>> &hHats, int *uu_hat);
	void Decoder_Ran(Modem_Linear_System &modem_linear_system, std::vector<std::complex<double>> &totalH, int *uu_hat, std::vector<int> &difflisterr);
	int Decoder_Ran_Hist(Modem_Linear_System &modem_linear_system, std::vector<std::complex<double>> &totalH, int *uu);

private:
	void Demmaping(Modem_Linear_System &modem_linear_system, std::vector<std::pair<int, std::complex<double>>> &thetaList);
	double adaptFunc(std::vector<std::complex<double>> &data, std::vector<std::complex<double>> &graySymbol, std::complex<double> &h, double var);
	void convertToPb(std::vector<double> &optsets);
	int getParityCheck(std::vector<std::pair<int, std::complex<double>>> &thetaList, int totalangle);
	int getParityCheck();

private:
	bool isCorrectInList(std::vector<std::pair<int, std::complex<double>>> &thetaList,
						 std::vector<std::pair<int, double>> &differ, std::vector<std::complex<double>> &totalH, int *uu, int list_len);

public:
	CBinaryLDPCCodec m_LDPC_codec;
	int **m_generator;	 // random generator matrix
	int *m_cc_Ran;		 // matrix multi results
	int *m_uu_Ran;		 // extra bits
	int *m_rr;			 // hard decision
	int m_extrabits_len; // extra bits len
	int m_list_count;
	int m_len_uu;	   // length of uu for LDPC
	int m_len_cc;	   // length of cc for LDPC
	double m_coderate; // code rate
	double *m_p0_cc;   // soft information
	double *bitLin;	   // bit input probability
	double *bitLout;   // bit out probability
};

#endif
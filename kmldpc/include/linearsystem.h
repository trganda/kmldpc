#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

#include <complex>

#include "modem.h"
#include "randnum.h"
#include "utility.h"

class Linear_System
{
public:
	double sigma;
	double var;

	double * noise;
	double * m_yy;
	double * m_Nyy;
	double *m_symbol_prob;

	int m_len_xx;
	int num_symbol;
		
	CModem * m_modem;
		
	void Malloc(CModem * modem, int len_xx, int code_no, char *file_name);
	void Free();
	
	void Soft_AWGN_Demodulation(double *yy, double *sym_prob);
	void AWGN_linear_system(double *xx, double *sym_prob);
	void Parition_HAWGN_system(double* xx, std::vector<std::complex<double>>& selectH);
};

#endif
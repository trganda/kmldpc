#ifndef MODEM_LINEAR_SYSTEM_H
#define MODEM_LINEAR_SYSTEM_H

#include "modem.h"
#include "utility.h"
#include "linear_system.h"
#include <complex>

class Modem_Linear_System
{
public:
	Modem_Linear_System(void);
	~Modem_Linear_System(void);

	int m_len_cc;
	int m_len_xx;

	double * m_xx;
	double * m_sym_prob;
	
	CModem m_modem;
	Linear_System Lin_Sym;
	char m_constellation_file[255];
	
	void Malloc(int len_cc, int code_no, char *file_name);
	void Free();	
	void modem_linear_system(int *cc, double *m_sym_prob);
	void modem_linear_system_parition(int* cc, std::vector<std::complex<double>>& selectH);
	void Soft_Demodulation(std::vector<std::pair<int, std::complex<double>>>& thetaList);
	std::vector<std::complex<double>> getRSymbol();
};

#endif
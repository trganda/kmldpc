#ifndef MODEM_H
#define MODEM_H

#include <iostream>
#include <complex>
#include <vector>

class CModem
{
public:
	int num_symbol{};
	int len_input{};
	int len_output{};

	int **input_symbol{};
	double **output_symbol{};

	double *m_symbol_prob{};

	double m_Es{};

	void Malloc(int code_no, char *file_name);
	void PrintCodeParameter(FILE *fp) const;
	void Free();

	void Mapping(int *cc, double *xx, int num_bit_in_blk);
	void Demapping(double *bitLin, double *symRin, double *bitLou, int num_bit_in_blk);

	std::vector<std::complex<double>> getConstellation();
};

#endif

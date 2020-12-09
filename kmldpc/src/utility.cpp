/*******************************************************************************
Programmed by Xiao Ma, Sun Yat-sen University. 
If you have any suggestions, please contact me at maxiao@mail.sysu.edu.cn
The program can only be employed for academic research.
******************************************************************************/

#include "utility.h"
#include "randnum.h"

extern CLCRandNum rndGen0;

void MatrixProd(int* uu, int* cc, int** G, int dim, int len) {
	for (int n = 0; n < len; ++n) {
		cc[n] = 0;
	}
	for (int k = 0; k < dim; ++k) {
		for (int n = 0; n < len; ++n) {
			cc[n] ^= (uu[k] & G[k][n]);
		}
	}
}

double Seqmax(double *x, int num)
{
	int i;
	double temp;

	temp = x[0];
	for (i = 1; i < num; i++)
	{
		if (x[i] > temp)
			temp = x[i];
	}
	return temp;
}

void ProbClip(double *xx, int len_xx)
{
	int i;

	for (i = 0; i < len_xx; i++){
		if (xx[i] < SMALLPROB)
			xx[i] = SMALLPROB;
		else if (xx[i] > 1.0 - SMALLPROB)
			xx[i] = 1.0 - SMALLPROB;
	}
}

std::ostream & operator<<( std::ostream & os,const std::complex<double> & c) {
    os << std::fixed << std::setprecision(14)
       << std::resetiosflags(std::ios::showpos)
       << '('
       << std::showpos << c.real() << std::showpos << c.imag() << 'i'
       << ')'
       << std::resetiosflags(std::ios::showpos);
    return os;
}
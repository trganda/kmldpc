#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>
#include <complex>
#include <iomanip>

#define m_PI 3.14159265358979
#define SMALLPROB 1.0e-12
#define DOUBLE_MAX 1.0e+300
#define SQRT2 1.4142135623730950488016

#define SMALL_LLR -28.0
#define LARGE_LLR 28.0

void MatrixProd(int* uu, int* cc, int** G, int dim, int len);

void ProbClip(double *xx, int len_xx);

double Seqmax(double *x, int num);

template <class ForwardIterator>
void ProbClip(ForwardIterator first, ForwardIterator last);

template<class ForwardIterator>
void ProbClip(ForwardIterator first, ForwardIterator last)
{
	while (first < last) {
		if (*first < SMALLPROB)
			*first = SMALLPROB;
		else if (*first > (1.0 - SMALLPROB))
			*first = 1.0 - SMALLPROB;
		first++;
	}
}

std::ostream & operator<<( std::ostream & os, const std::complex<double> & c);

#endif
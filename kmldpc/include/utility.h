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

#ifndef MAX
#define MAX(x, y) (((x)>(y))?(x):(y))
#endif

void MatrixProd(int* uu, int* cc, int** G, int dim, int len);

int sgn(double x);

void ProbClip(double *xx, int len_xx);

int min(int x, int y);

int max(int x, int y);

double Seqmax(double *x, int num);

//interleaver
void RandomIntl(int *pai, int period);
void N1N2RandomIntl(int *pai, int period, int N1, int N2);//|pai[i]-pai[j]| < N1 whenever 0<|i-j|>=N2
void ShiftIntl(int *pai, int shifted_pos, int period);

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

std::ostream & operator<<( std::ostream & os,const std::complex<double> & c);

#endif
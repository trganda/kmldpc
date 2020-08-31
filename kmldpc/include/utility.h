#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>

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

double FunctionQ(double x);
double erf_inv(double x);

int sgn(double x);

void ProbClip(double *xx, int len_xx);

void Dec2Bin(int d, int *b, int len_b);

void SeqDec2Bin(int *bin, int *dec, int len_dec, int len_symbol);

void SeqBin2Dec(int *bin, int *dec, int len_dec, int len_symbol);

int BitDotProd(int a, int b, int len);

int Systemizer(int num_row, int num_col, int **H, int **sysH, int *pai);

int min(int x, int y);

int max(int x, int y);

double Seqmax(double *x, int num);

double Seqmin(double *x, int num);

//interleaver
void RandomIntl(int *pai, int period);
void N1N2RandomIntl(int *pai, int period, int N1, int N2);//|pai[i]-pai[j]| < N1 whenever 0<|i-j|>=N2
void ShiftIntl(int *pai, int shifted_pos, int period);

void LLRClip(double *xx, int len_xx);//added by liangchulong @2013-01-30-11-59
double maxstar(double x, double y);

void DFT_gen(int size_t, double** F_re, double** F_im);//added by Leijun Wang @20151113

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

#endif
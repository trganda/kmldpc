#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>
#include <complex>
#include <iomanip>

const double m_PI = 3.14159265358979;
const double SMALLPROB = 1.0e-12;
const double DOUBLE_MAX = 1.0e+300;
const double SQRT2 = 1.4142135623730950488016;

const double SMALL_LLR = -28.0;
const double LARGE_LLR = 28.0;

void MatrixProd(int *uu, int *cc, int **G, int dim, int len);

void ProbClip(double *xx, int len_xx);

double Seqmax(double *x, int num);

template<class ForwardIterator>
void ProbClip(ForwardIterator first, ForwardIterator last);

std::ostream &operator<<(std::ostream &os, const std::complex<double> &c);

#endif
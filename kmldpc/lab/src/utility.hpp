#ifndef LAB_UTILITY_HPP
#define LAB_UTILITY_HPP

#include <string>
#include <vector>
#include <complex>
#include <iomanip>

namespace lab {

const double kPi = 3.14159265358979;
const double kSmallestProb = 1.0e-12;
const double kMaxDouble = 1.0e+300;
const double kSqrt2 = 1.4142135623730950488016;

const double kSmallestLLR = -28.0;
const double kLargestLLR = 28.0;

// Tanner Edge for LDPC
typedef struct Edge
{
    int m_row_no;
    int m_col_no;

    double m_alpha[2];
    double m_beta[2];
    double m_v2c[2];
    double m_c2v[2];

    struct Edge *left;
    struct Edge *right;
    struct Edge *up;
    struct Edge *down;
} Edge;

namespace utility {

// using inline to avoid multi definition of function
inline void MatrixProd(int *uu, int *cc, int **G, int dim, int len) {
    for (int n = 0; n < len; ++n) {
        cc[n] = 0;
    }
    for (int k = 0; k < dim; ++k) {
        for (int n = 0; n < len; ++n) {
            cc[n] ^= (uu[k] & G[k][n]);
        }
    }
}

inline double Seqmax(double *x, int num) {
    double temp = x[0];
    for (int i = 1; i < num; i++) {
        if (x[i] > temp)
            temp = x[i];
    }
    return temp;
}

inline void ProbClip(double *xx, int len_xx) {
    for (int i = 0; i < len_xx; i++) {
        if (xx[i] < kSmallestProb) {
            xx[i] = kSmallestProb;
        } else if (xx[i] > 1.0 - kSmallestProb) {
            xx[i] = 1.0 - kSmallestProb;
        }
    }
}

//
//template<class ForwardIterator>
//void ProbClip(ForwardIterator first, ForwardIterator last) {
//    while (first < last) {
//        if (*first < kSmallestProb)
//            *first = kSmallestProb;
//        else if (*first > (1.0 - kSmallestProb))
//            *first = 1.0 - kSmallestProb;
//        first++;
//    }
//}
//

inline std::ostream &operator<<(std::ostream &os, const std::complex<double> &c) {
    os << std::fixed << std::setprecision(14)
       << std::resetiosflags(std::ios::showpos)
       << '('
       << std::showpos << c.real() << std::showpos << c.imag() << 'i'
       << ')'
       << std::resetiosflags(std::ios::showpos);
    return os;
}

}

}

#endif
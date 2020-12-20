#ifndef LAB_UTILITY_H
#define LAB_UTILITY_H

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
    typedef struct Edge {
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
        void MatrixProd(int *uu, int *cc, int **G, int dim, int len);

        double Seqmax(double *x, int num);

        void ProbClip(double *xx, int len_xx);

        std::ostream &operator<<(std::ostream &os, const std::complex<double> &c);

    }

}

#endif
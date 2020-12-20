#ifndef LAB_LINEAR_SYSTEM_H
#define LAB_LINEAR_SYSTEM_H

#include <complex>

#include "modem.h"
#include "randnum.h"
#include "utility.h"

namespace lab {

    class LinearSystem {
    public:
        explicit LinearSystem();

        virtual ~LinearSystem();

        void Malloc(Modem *modem, int len_xx);

        void SoftAWGNDemodulation(const double *yy, double *sym_prob,
                                  std::complex<double> &theta) const;

        void AWGNLinearSystem(const double *xx, double *sym_prob) const;

        void PartitionHAWGNSystem(const double *xx,
                                  std::vector<std::complex<double>> &selectH) const;

        // Setter
        void SetSigma(double sigma);
        void SetVar(double var);

        // Getter
        double *GetYy() const;
        double *GetNyy() const;
        Modem *GetMModem() const;

    private:
        double sigma_;
        double var_;

        double *noise_;
        double *yy_;
        double *nyy_;
        double *symbol_prob_;

        int xx_len_;
        int num_symbol_;

        Modem *modem_;
    };

}

#endif
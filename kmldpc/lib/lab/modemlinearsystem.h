#ifndef LAB_MODEM_LINEAR_SYSTEM_HPP
#define LAB_MODEM_LINEAR_SYSTEM_HPP

#include <complex>
#include <cstring>

#include "modem.h"
#include "utility.h"
#include "randnum.h"
#include "linearsystem.h"
#include "toml.hpp"

namespace lab {

    class ModemLinearSystem {
    public:
        explicit ModemLinearSystem();
        virtual ~ModemLinearSystem();

        void Malloc(int len_cc, const toml::value& arguments);
        void MLSystem(int *cc, double *sym_prob);
        void MLSystemPartition(int *cc,
                               std::vector<std::complex<double>> &selectH) const;
        void SoftDemodulation(std::vector<std::pair<int, std::complex<double>>> &thetaList) const;

        std::vector<std::complex<double>> GetRSymbol() const;
        double *GetSymProb() const;
        const Modem &GetModem() const;
        LinearSystem &GetLinearSystem();

    private:
        int cc_len_;
        int xx_len_;
        double *xx_;
        double *sym_prob_;
        Modem modem_;
        LinearSystem linsym_;
    };

}

#endif
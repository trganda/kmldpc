#ifndef LAB_MODEM_H
#define LAB_MODEM_H

#include <iostream>
#include <complex>
#include <vector>

#include "log.h"
#include "utility.h"

namespace lab {

    class Modem {
    public:
        explicit Modem();
        virtual ~Modem();
    public:
        void Malloc(int code_no, const std::string& modem_file);
        void Mapping(const int *bin_cc, double *xx, int num_bit_in_blk) const;
        void DeMapping(double *bitLin, double *symRin,
                       double *bitLout, int num_bit_in_blk) const;

        std::vector<std::complex<double>> GetConstellations() const;
        int GetNumSymbol() const;
        int GetInputLen() const;
        int GetOutputLen() const;
        double **GetOutputSymbol() const;
    private:
        // Number of constellation points
        int num_symbol_;
        int input_len_;
        int output_len_;

        int **input_symbol_;
        double **output_symbol_;
        double *symbol_prob_;

        double energy_;
    };

}

#endif

#ifndef LAB_MODEM_H
#define LAB_MODEM_H

#include <iostream>
#include <complex>
#include <vector>

#include "toml.hpp"

#include "log.h"
#include "utility.h"

namespace lab {

    class Modem {
    public:
        explicit Modem(const toml::value& arguments);
        virtual ~Modem();
    public:
        void Malloc(const std::string& modem_file);
        void Mapping(const int *bin_cc, double *xx, int num_bit_in_blk) const;
        std::vector<std::complex<double>> Mapping(const int* bin_cc, int bin_cc_len);
        void DeMapping(double *bitLin, double *symRin,
                       double *bitLout, int num_bit_in_blk) const;

        std::vector<std::complex<double>> GetConstellations() const;
        int GetNumSymbol() const;
        int GetInputLen() const;
        int GetOutputLen() const;
        double **GetOutputSymbol() const;

    private:
        void init(const std::string& modem_file);
    private:
        std::vector<std::vector<int>> symbol_in_;
        std::vector<std::complex<double>> symbol_out_;
        // Number of constellation points
        int symbol_num_;
        int input_len_;
        int output_len_;

        int **input_symbol_;
        double **output_symbol_;
        double *symbol_prob_;

        double energy_;
    };

}

#endif

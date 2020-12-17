#ifndef LAB_MODEM_LINEAR_SYSTEM_HPP
#define LAB_MODEM_LINEAR_SYSTEM_HPP

#include <complex>
#include <cstring>

#include "modem.hpp"
#include "utility.hpp"
#include "randnum.hpp"
#include "linearsystem.hpp"

namespace lab {

    class ModemLinearSystem {
    public:
        explicit ModemLinearSystem()
                : cc_len_(0), xx_len_(0),
                  xx_(nullptr), sym_prob_(nullptr),
                  modem_(Modem()), linsym_(LinearSystem()) {}

        virtual ~ModemLinearSystem() {
            delete[]xx_;
            delete[]sym_prob_;
        }

        void Malloc(int len_cc, int code_no, char *file_name) {
            int temp0, temp;
            char temp_str[80] = {' '};
            char mark[80];
            char constellation_file_[255];
            FILE *fp;

            cc_len_ = len_cc;

            if ((fp = fopen(file_name, "r")) == nullptr) {
                fprintf(stderr, "\nCannot Open %s", file_name);
                exit(3);
            }

            sprintf(mark, "ModemLinearSystem***%d***PARAMETERS", code_no);
            while (strcmp(temp_str, mark) != 0)
                fscanf(fp, "%s", temp_str);

            fscanf(fp, "%s", temp_str);
            fscanf(fp, "%s", constellation_file_);
            fclose(fp);

            modem_.Malloc(0, constellation_file_);

            temp0 = cc_len_ / modem_.GetInputLen();
            if (cc_len_ % modem_.GetInputLen() != 0) {
                fprintf(stderr, "\n(cc_len_ = %d) %% (input_len_ = %d) != 0 !\n", cc_len_, modem_.GetInputLen());
                system("pause");
                exit(3);
            }

            xx_len_ = temp0 * modem_.GetOutputLen();
            xx_ = new double[xx_len_];

            temp = temp0 * modem_.GetNumSymbol();
            sym_prob_ = new double[temp];

            linsym_.Malloc(&modem_, xx_len_, 0, file_name);
        }

        void MLSystem(int *cc, double *sym_prob) {
            modem_.Mapping(cc, xx_, cc_len_);

            linsym_.AWGNLinearSystem(xx_, sym_prob);
        }

        void MLSystemPartition(int *cc,
                               std::vector<std::complex<double>> &selectH) const {
            modem_.Mapping(cc, xx_, cc_len_);

            linsym_.PartitionHAWGNSystem(xx_, selectH);
        }

        void SoftDemodulation(std::vector<std::pair<int, std::complex<double>>> &thetaList) const {
            std::vector<std::complex<double>> um = GetRSymbol();

            int symbolPerPart = um.size() / thetaList.size();

            std::vector<std::complex<double>> tum(um.size());
            for (size_t i = 0; i < thetaList.size(); i++) {
                for (int j = 0; j < symbolPerPart; j++) {
                    tum[j + i * symbolPerPart] = um[j + i * symbolPerPart];
                    linsym_.GetYy()[0] = tum[j + i * symbolPerPart].real();
                    linsym_.GetYy()[1] = tum[j + i * symbolPerPart].imag();
                    int temp = (j + i * symbolPerPart) * linsym_.GetMModem()->GetNumSymbol();
                    linsym_.SoftAWGNDemodulation(linsym_.GetYy(), (sym_prob_ + temp), thetaList[i].second);
                }
            }
        }

        std::vector<std::complex<double>> GetRSymbol() const {
            std::vector<std::complex<double>> um(xx_len_ / 2);
            for (size_t i = 0; i < um.size(); i++) {
                um[i] = std::complex<double>(linsym_.GetNyy()[i * 2], linsym_.GetNyy()[i * 2 + 1]);
            }
            return um;
        }

        double *GetSymProb() const {
            return sym_prob_;
        }

        const Modem &GetModem() const {
            return modem_;
        }

        LinearSystem &GetLinearSystem() {
            return linsym_;
        }

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
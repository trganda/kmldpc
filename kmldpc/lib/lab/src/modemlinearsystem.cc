#include "modemlinearsystem.h"

namespace lab {

    ModemLinearSystem::ModemLinearSystem(const toml::value& arguments, int cc_len)
            : cc_len_(cc_len), xx_len_(0),
              xx_(nullptr), sym_prob_(nullptr),
              modem_(Modem()), linsym_(LinearSystem()) {
        Malloc(cc_len, arguments);
    }

    ModemLinearSystem::~ModemLinearSystem() {
        delete[]xx_;
        delete[]sym_prob_;
    }

    void ModemLinearSystem::Malloc(int len_cc, const toml::value &arguments) {
        int temp0, temp;

        cc_len_ = len_cc;

        const auto modem = toml::find(arguments, "modem");
        std::string modem_file = toml::find<std::string>(modem, "modem_file");

        modem_.Malloc(0, modem_file);

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

        linsym_.Malloc(&modem_, xx_len_);
    }

    void ModemLinearSystem::MLSystem(int *cc, double *sym_prob) {
        modem_.Mapping(cc, xx_, cc_len_);

        linsym_.AWGNLinearSystem(xx_, sym_prob);
    }

    void ModemLinearSystem::MLSystemPartition(int *cc,
                                              std::vector<std::complex<double>> &selectH) const {
        modem_.Mapping(cc, xx_, cc_len_);

        linsym_.PartitionHAWGNSystem(xx_, selectH);
    }

    void ModemLinearSystem::SoftDemodulation(std::vector<std::pair<int, std::complex<double>>> &thetaList) const {
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

    std::vector<std::complex<double>> ModemLinearSystem::GetRSymbol() const {
        std::vector<std::complex<double>> um(xx_len_ / 2);
        for (size_t i = 0; i < um.size(); i++) {
            um[i] = std::complex<double>(linsym_.GetNyy()[i * 2], linsym_.GetNyy()[i * 2 + 1]);
        }
        return um;
    }

    double *ModemLinearSystem::GetSymProb() const {
        return sym_prob_;
    }

    const Modem &ModemLinearSystem::GetModem() const {
        return modem_;
    }

    LinearSystem &ModemLinearSystem::GetLinearSystem() {
        return linsym_;
    }

}
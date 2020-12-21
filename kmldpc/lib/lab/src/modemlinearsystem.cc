#include "modemlinearsystem.h"

namespace lab {

    ModemLinearSystem::ModemLinearSystem(const toml::value& arguments, int cc_len)
            : Modem(arguments), cc_len_(cc_len), xx_len_(0),
              xx_(nullptr), sym_prob_(nullptr),
              modem_(Modem(arguments)), linsym_(LinearSystem()) {
        int temp0, temp;

        temp0 = cc_len_ / modem_.GetInputLen();
        if (cc_len_ % modem_.GetInputLen() != 0) {
            LOG(logger::Error, true) << "(cc_len_ = " << cc_len_ << " %% (input_len_ = "
                                               << GetInputLen() << " ) != 0 !" << std::endl;
            exit(-1);
        }

        xx_len_ = temp0 * modem_.GetOutputLen();
        xx_ = new double[xx_len_];



        temp = temp0 * modem_.GetNumSymbol();
        sym_prob_ = new double[temp];

        linsym_.Malloc(&modem_, xx_len_);
    }

    ModemLinearSystem::~ModemLinearSystem() {
        delete[]xx_;
        delete[]sym_prob_;
    }

    void ModemLinearSystem::MLSystem(int *cc, double *sym_prob) {
        modem_.Mapping(cc, xx_, cc_len_);

        linsym_.AWGNLinearSystem(xx_, sym_prob);
    }

    void ModemLinearSystem::MLSystemPartition(const int *cc,
                                              std::vector<std::complex<double>> &selectH) {
//        modem_.Mapping(cc, xx_, cc_len_);
        auto ret = modem_.Mapping(cc, cc_len_);
        linsym_.PartitionHAWGNSystem(ret, selectH);
//        linsym_.PartitionHAWGNSystem(xx_, selectH);
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
//                linsym_.SoftAWGNDemodulation(um[j+i*symbolPerPart], sym_prob_+temp, thetaList[i].second);
            }
        }
    }

    std::vector<std::complex<double>> ModemLinearSystem::GetRSymbol() const {
//        std::vector<std::complex<double>> um(xx_len_ / 2);
//        for (size_t i = 0; i < um.size(); i++) {
//            um[i] = std::complex<double>(linsym_.GetNyy()[i * 2], linsym_.GetNyy()[i * 2 + 1]);
//        }
//        return um;
        return linsym_.GetTyy();
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
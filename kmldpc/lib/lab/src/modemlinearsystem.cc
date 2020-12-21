#include "modemlinearsystem.h"

namespace lab {

    ModemLinearSystem::ModemLinearSystem(const toml::value &arguments, int cc_len)
            : Modem(arguments), cc_len_(cc_len), sigma_(0.0),
              var_(0.0), sym_prob_(nullptr) {
        if (cc_len_ % input_len_ != 0) {
            LOG(logger::Error, true) << "(cc_len_ = " << cc_len_ << " %% (input_len_ = "
                                     << input_len_ << " ) != 0 !" << std::endl;
            exit(-1);
        }
        yy_ = std::vector<std::complex<double>>(cc_len_ / input_len_);
        sym_prob_ = new double[(cc_len_ / input_len_) * symbol_num_];
    }

    ModemLinearSystem::~ModemLinearSystem() {
        delete[]sym_prob_;
    }

    void ModemLinearSystem::PartitionModemLSystem(const int *cc,
                                                  std::vector<std::complex<double>> &select_h) {
        auto ret = Mapping(cc, cc_len_);
        PartitionHAWGNSystem(ret, select_h);
    }

    void ModemLinearSystem::PartitionHAWGNSystem(
            std::vector<std::complex<double>> &xx,
            std::vector<std::complex<double>> &selected_h) {

        std::vector<std::complex<double>> noise(xx.size());
        CLCRandNum::Get().Normal(noise);

        for (size_t i = 0; i < selected_h.size(); i++) {
            auto num_of_part = xx.size() / selected_h.size();
            for (size_t j = i * num_of_part; j < num_of_part; j++) {
                std::complex<double> temp = xx[j] * selected_h[i];
                yy_[j] = temp + noise[j] * std::complex<double>(sigma_ / kSqrt2, 0);
            }
        }
    }

    void ModemLinearSystem::SoftAWGNDemodulation(const std::complex<double> &yy, double *sym_prob,
                                                 std::complex<double> &theta_h) const {
        std::vector<double> symbol_prob(symbol_num_);
        double sqr_norm = 0.0;

        for (size_t i = 0; i < symbol_prob.size(); i++) {
            auto symbol = symbol_out_[i];
            symbol *= theta_h;
            symbol -= yy;
            sqr_norm = (symbol.real() * symbol.real() + symbol.imag() * symbol.imag()) / var_;
            symbol_prob[i] = -sqr_norm;
        }

        auto max_prob = std::max_element(symbol_prob.begin(), symbol_prob.end());
        for (int i = 0; i < symbol_num_; i++) {
            symbol_prob[i] = exp(symbol_prob[i] - *max_prob);
        }
        // normalization
        double sum = 0.0;
        for (int i = 0; i < symbol_num_; i++) {
            sum += symbol_prob[i];
        }

        for (int i = 0; i < symbol_num_; i++) {
            symbol_prob[i] /= sum;
        }

        for (int i = 0; i < symbol_num_; i++) {
            sym_prob[i] = symbol_prob[i];
        }

        utility::ProbClip(sym_prob, symbol_num_);
    }

    void ModemLinearSystem::SoftDemodulation(std::vector<std::pair<int, std::complex<double>>> &thetaList) const {
        auto symbolPerPart = yy_.size() / thetaList.size();

        for (size_t i = 0; i < thetaList.size(); i++) {
            for (size_t j = 0; j < symbolPerPart; j++) {
                auto temp = (j + i * symbolPerPart) * symbol_num_;
                SoftAWGNDemodulation(yy_[j + i * symbolPerPart], sym_prob_ + temp, thetaList[i].second);
            }
        }
    }

    void ModemLinearSystem::DeMapping(std::vector<std::pair<int, std::complex<double>>> &thetaList, double *bitLin,
                                      double *bitLout) {
        SoftDemodulation(thetaList);
        Modem::DeMapping(bitLin, sym_prob_, bitLout, yy_.size());
    }

    std::vector<std::complex<double>> ModemLinearSystem::GetRecvSymbol() const {
        return yy_;
    }

    void ModemLinearSystem::SetSigma(double sigma) {
        this->sigma_ = sigma;
    }

    void ModemLinearSystem::SetVar(double var) {
        this->var_ = var;
    }

}
#include "linearsystem.h"

namespace lab {
    LinearSystem::LinearSystem()
            : sigma_(0.0), var_(0.0),
              noise_(nullptr), yy_(nullptr), nyy_(nullptr), symbol_prob_(nullptr),
              xx_len_(0), num_symbol_(0), modem_(nullptr) {}

    LinearSystem::~LinearSystem() {
        delete[]noise_;
        delete[]yy_;
        delete[]nyy_;
        delete[]symbol_prob_;
    }

    void LinearSystem::Malloc(Modem *modem, int len_xx) {
        modem_ = modem;
        xx_len_ = len_xx;

        noise_ = new double[xx_len_];
        yy_ = new double[2];
        nyy_ = new double[xx_len_];

        symbol_prob_ = new double[modem_->GetNumSymbol()];

        num_symbol_ = modem_->GetNumSymbol();
    }

    void LinearSystem::SoftAWGNDemodulation(const double *yy, double *sym_prob,
                                            std::complex<double> &theta) const {
        int i, q;
        double sqr_norm, sum;
        double temp;

        temp = 0.0;
        for (i = 0; i < num_symbol_; i++) {
            sqr_norm = 0.0;

            double symbol_x = modem_->GetOutputSymbol()[i][0];
            double symbol_y = modem_->GetOutputSymbol()[i][1];
            double x_real = symbol_x * theta.real() - symbol_y * theta.imag();
            double x_imag = symbol_y * theta.real() + symbol_x * theta.imag();


            sqr_norm += ((yy[0] - x_real) * (yy[0] - x_real)
                         + (yy[1] - x_imag) * (yy[1] - x_imag)) / var_;

            symbol_prob_[i] = -sqr_norm;

        }

        temp = utility::Seqmax(symbol_prob_, num_symbol_);

        for (q = 0; q < num_symbol_; q++) {
            symbol_prob_[q] = exp(symbol_prob_[q] - temp);
        }

        //normalization
        sum = 0.0;
        for (q = 0; q < num_symbol_; q++) {
            sum += symbol_prob_[q];
        }

        for (q = 0; q < num_symbol_; q++) {
            symbol_prob_[q] /= sum;
        }

        for (q = 0; q < num_symbol_; q++) {
            sym_prob[q] = symbol_prob_[q];
        }

        utility::ProbClip(sym_prob, num_symbol_);
    }

    void LinearSystem::AWGNLinearSystem(const double *xx, double *sym_prob) const {
        int i;
        int xxidx, nnidx;
        int temp;

        int num_of_symbol_blk = xx_len_ / 2;

        CLCRandNum::Get().Normal(noise_, xx_len_);

        xxidx = 0;
        nnidx = 0;
        for (i = 0; i < num_of_symbol_blk; i++) {
            yy_[0] = xx[xxidx] + (sigma_ / kSqrt2) * noise_[nnidx];
            yy_[1] = xx[xxidx + 1] + (sigma_ / kSqrt2) * noise_[nnidx + 1];
            nyy_[xxidx] = yy_[0];
            nyy_[xxidx + 1] = yy_[1];
            xxidx += 2;
            nnidx += 2;

            temp = i * modem_->GetNumSymbol();
            std::complex<double> h(1, 0);
            SoftAWGNDemodulation(yy_, (sym_prob + temp), h);
        } //end of for (i = 0; i < num_sym_in_blk; i++)
    }

    void LinearSystem::PartitionHAWGNSystem(const double *xx,
                                            std::vector<std::complex<double>> &selectH) const {
        int xxidx, nnidx;

        int num_of_symbol_blk = xx_len_ / 2;
        int partLen = selectH.size();

        CLCRandNum::Get().Normal(noise_, xx_len_);

        for (int i = 0; i < partLen; i++) {
            int num_of_symbol_part_blk = num_of_symbol_blk / partLen;
            xxidx = 0 + i * num_of_symbol_part_blk * 2;
            nnidx = 0 + i * num_of_symbol_part_blk * 2;
            for (int j = 0; j < num_of_symbol_part_blk; j++) {
                double treal = xx[xxidx] * selectH[i].real() - xx[xxidx + 1] * selectH[i].imag();
                double timag = xx[xxidx + 1] * selectH[i].real() + xx[xxidx] * selectH[i].imag();
                nyy_[xxidx] = treal + (sigma_ / kSqrt2) * noise_[nnidx];
                nyy_[xxidx + 1] = timag + (sigma_ / kSqrt2) * noise_[nnidx + 1];

                xxidx += 2;
                nnidx += 2;
            }
        }
    }

    // Setter
    void LinearSystem::SetSigma(double sigma) {
        this->sigma_ = sigma;
    }

    void LinearSystem::SetVar(double var) {
        this->var_ = var;
    }

    // Getter
    double *LinearSystem::GetYy() const {
        return yy_;
    }

    double *LinearSystem::GetNyy() const {
        return nyy_;
    }

    Modem *LinearSystem::GetMModem() const {
        return modem_;
    }
}
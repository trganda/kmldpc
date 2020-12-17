#ifndef LAB_SOURCE_SINK_HPP
#define LAB_SOURCE_SINK_HPP

#include <complex>
#include <vector>
#include <iomanip>

#include "log.hpp"
#include "randnum.hpp"

namespace lab {

    class CSourceSink {
    public:
        void GetBitStr(int *uu, int len) {
            for (int t = 0; t < len; t++) {
                uu[t] = (CLCRandNum::Get().Uniform() < 0.5 ? 0 : 1);
            }
        }

        void GetSymStr(int *uu, int qary, int len) {
            for (int t = 0; t < len; t++) {
                uu[t] = qary;
                while (uu[t] == qary)
                    uu[t] = (int) (qary * CLCRandNum::Get().Uniform());
            }
        }

        void ClrCnt() {
            num_tot_blk_ = 0;
            num_tot_bit_ = 0;
            num_err_blk_ = 0;
            num_err_bit_ = 0;
        }

        void CntErr(const int *uu, const int *uu_hat,
                    int len, int accumulator) {
            temp_err_ = 0;
            for (int t = 0; t < len; t++) {
                if (uu_hat[t] != uu[t])
                    temp_err_++;
            }

            if (accumulator == 1) {
                if (temp_err_ > 0) {
                    num_err_bit_ += temp_err_;
                    num_err_blk_ += 1;
                }

                num_tot_blk_ += 1.0;
                num_tot_bit_ += len;

                ber_ = num_err_bit_ / num_tot_bit_;
                fer_ = num_err_blk_ / num_tot_blk_;
            }
        }

        void PrintResult(double snr) const {
            LOG(logger::Info, true) << std::fixed << std::setprecision(3) << std::setfill('0')
                                    << "SNR = "
                                    << std::setw(3) << std::right << snr << ' '
                                    << "Total blk = "
                                    << std::setw(7) << std::right << std::setprecision(0) << num_tot_blk_ << ' '
                                    << "Error blk = "
                                    << std::setw(7) << std::right << num_err_blk_ << ' '
                                    << "Error bit = "
                                    << std::setw(7) << std::right << num_err_bit_ << ' '
                                    << std::fixed << std::setprecision(14)
                                    << "BER = " << ber_ << ' '
                                    << "FER = " << fer_
                                    << std::endl;
        }

        void PrintResult(FILE *fp) const {
            fprintf(fp, "tot_blk = %d: err_blk = %d: err_bit = %d: ber = %12.10lf: fer = %12.10lf\n",
                    (int) num_tot_blk_, num_err_blk_, num_err_bit_, ber_, fer_);
        }

        // Getter
        double GetNumTotBlk() const {
            return num_tot_blk_;
        }

        int GetNumErrBlk() const {
            return num_err_blk_;
        }

        double GetBer() const {
            return ber_;
        }

        double GetFer() const {
            return fer_;
        }

    private:
        double num_tot_blk_;
        double num_tot_bit_;
        int num_err_blk_;
        int num_err_bit_;
        int temp_err_;
        double ber_;
        double fer_;
    };

}

#endif
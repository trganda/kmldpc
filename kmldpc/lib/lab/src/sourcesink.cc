#include "sourcesink.h"

namespace lab {
void SourceSink::GetBitStr(int *uu, int len) {
  for (int t = 0; t < len; t++) {
    uu[t] = (CLCRandNum::Get().Uniform() < 0.5 ? 0 : 1);
  }
}

void SourceSink::GetSymStr(int *uu, int qary, int len) {
  for (int t = 0; t < len; t++) {
    uu[t] = qary;
    while (uu[t] == qary)
      uu[t] = (int) (qary * CLCRandNum::Get().Uniform());
  }
}

void SourceSink::ClrCnt() {
  num_tot_blk_ = 0;
  num_tot_bit_ = 0;
  num_err_blk_ = 0;
  num_err_bit_ = 0;
}

void SourceSink::CntErr(const int *uu, const int *uu_hat,
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

void SourceSink::PrintResult(double snr) const {
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

double SourceSink::num_tot_blk() const {
  return num_tot_blk_;
}

int SourceSink::num_err_blk() const {
  return num_err_blk_;
}

double SourceSink::ber() const {
  return ber_;
}

double SourceSink::fer() const {
  return fer_;
}

}

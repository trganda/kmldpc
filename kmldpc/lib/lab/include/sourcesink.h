#ifndef LAB_SOURCE_SINK_H
#define LAB_SOURCE_SINK_H

#include <complex>
#include <vector>
#include <iomanip>

#include "log.h"
#include "randnum.h"

namespace lab {

class SourceSink {
 public:
  void GetBitStr(int *uu, int len);
  void GetSymStr(int *uu, int qary, int len);
  void ClrCnt();
  void CntErr(const int *uu, const int *uu_hat,
              int len, int accumulator);

  void PrintResult(double snr) const;
  // Getter
  int num_err_blk() const;
  double num_tot_blk() const;
  double ber() const;
  double fer() const;

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
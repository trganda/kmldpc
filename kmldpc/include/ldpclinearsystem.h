#ifndef KMLDPC_LDPC_LINEAR_SYSTEM_H
#define KMLDPC_LDPC_LINEAR_SYSTEM_H

#include <complex>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <thread>

#include "sourcesink.h"
#include "xorsegcodec.h"
#include "modemlinearsystem.h"
#include "randnum.h"
#include "log.h"
#include "toml.hpp"

#include "kmeans.h"

typedef struct CodecData {
  CodecData(const CodecData &codec_data) {
    uu_len_ = codec_data.uu_len_;
    cc_len_ = codec_data.cc_len_;
    uu_ = new int[uu_len_];
    uu_hat_ = new int[uu_len_];
    cc_ = new int[cc_len_];
    cc_hat_ = new int[cc_len_];
  }

  CodecData(int uu_len, int cc_len) {
    uu_len_ = uu_len;
    cc_len_ = cc_len;
    uu_ = new int[uu_len_];
    uu_hat_ = new int[uu_len_];
    cc_ = new int[cc_len_];
    cc_hat_ = new int[cc_len_];
  }

  ~CodecData() {
    delete[] uu_;
    delete[] uu_hat_;
    delete[] cc_;
    delete[] cc_hat_;
  }

  // Uncoded codeword
  int *uu_;
  int *uu_hat_;
  int uu_len_;
  // Encoded codeword
  int *cc_;
  int *cc_hat_;
  int cc_len_;
} CodecData;

class LDPCLinearSystem {
 public:
  explicit LDPCLinearSystem(toml::value arguments);
  virtual ~LDPCLinearSystem() = default;

  void Simulator();
 private:
  void Run(lab::XORSegCodec &codec, lab::ModemLinearSystem &mls, lab::SourceSink &ssink, CodecData &cdata,
           double snr, bool histogram_enable,
           std::pair<double, double> &ber, std::pair<double, double> &fer) const;

 private:
  const toml::value arguments_;
  // Simulation range of snr
  double min_snr_;
  double max_snr_;
  // Step size for increase snr
  double step_snr_;
  // Maximum error blocks, stop simulation of one snr while
  // errors block is equal to it
  int max_err_blk_;
  // Maximum blocks for simulation
  int max_num_blk_;

  lab::XORSegCodec codec_;
  CodecData codec_data_;
  lab::ModemLinearSystem modem_linear_system_;
};

#endif
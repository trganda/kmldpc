#ifndef LAB_XOR_SEG_CODEC_H
#define LAB_XOR_SEG_CODEC_H

#include "binary5gldpccodec.h"
#include "binaryldpccodec.h"
#include "log.h"
#include "mat.h"
#include "modemlinearsystem.h"
#include "randnum.h"
#include "utility.h"
#include <algorithm>
#include <complex>
#include <iomanip>
#include <vector>

class KmCodec {
 public:
  KmCodec() = default;
  KmCodec(const KmCodec &codec);
  explicit KmCodec(toml::value arguments);
  virtual ~KmCodec();
  void Encoder(int *uu, int *cc);
  void Decoder(
      lab::ModemLinearSystem &modem_linear_system, const std::vector<std::complex<double>> &h_hats,
      int *uu_hat);
  std::vector<double> GetHistogramData(
      lab::ModemLinearSystem &mlsystem,
      const std::vector<std::complex<double>> &h_hats, int *uu_hat);
  int uu_len() const;
  int cc_len() const;

 private:
  void DeMapping(
      lab::ModemLinearSystem &modem_linear_system,
      std::vector<std::pair<int, std::complex<double>>> &thetaList) const;
  int GetParityCheck() const;
  std::vector<double> GetMetrics(
      lab::ModemLinearSystem &modem_linear_system,
      const std::vector<std::complex<double>> &h_hats, int *uu_hat);
  double Metric(int *uu_hat);

 private:
  const toml::value arguments_;
  std::unique_ptr<lab::BinaryLDPCCodec> ldpc_codec_;
  int iter_cnt_;// iteration times while using LDPC on 5G
  bool using_ldpc_5g_;
  bool using_syndrom_metric_;
  int uu_len_;       // length of uu for LDPC
  int cc_len_;       // length of cc for LDPC
  int *rr_;          // hard decision
  double *bit_l_in_; // bit input probability
  double *bit_l_out_;// bit out probability
};
#endif
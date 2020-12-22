#ifndef LAB_MODEM_H
#define LAB_MODEM_H

#include <iostream>
#include <complex>
#include <vector>

#include "toml.hpp"

#include "log.h"
#include "utility.h"

namespace lab {

class Modem {
 public:
  explicit Modem(const toml::value &arguments);
  virtual ~Modem() = default;
 public:
  void Mapping(const int *bin_cc, int bin_cc_len, std::vector<std::complex<double>> &xx);

  void DeMapping(double *bitLin, double *symRin,
                 double *bitLout, int yy_len) const;

  std::vector<std::complex<double>> GetConstellations() const;
 private:
  void init(const std::string &modem_file);
 protected:
  // Number of constellation points
  int symbol_num_;
  int input_len_;
  int output_len_;
  std::vector<std::complex<double>> symbol_out_;
 private:
  std::vector<std::vector<int>> symbol_in_;
};

}

#endif

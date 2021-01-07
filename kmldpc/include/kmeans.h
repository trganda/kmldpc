#ifndef KMLDPC_KMEANS_H
#define KMLDPC_KMEANS_H

#include "mat.h"
#include <algorithm>
#include <complex>
#include <iostream>
#include <utility>
#include <vector>

namespace kmldpc {
class KMeans {
 public:
  explicit KMeans(
      std::vector<std::complex<double>> &data,
      std::vector<std::complex<double>> &constellations,
      int iter);
  ~KMeans();

 public:
  std::vector<std::complex<double>> clusters();
  std::vector<int> idx();
  void Run();
#ifdef USE_MATIO
  void DumpToMat(std::string &filename, std::vector<std::complex<double>> &append);
#endif
 private:
  const std::vector<std::complex<double>> data_;
  const std::vector<std::complex<double>> constellations_;
  std::vector<std::complex<double>> clusters_;
  std::vector<int> idx_;
  unsigned int iter_;
};
}// namespace kmldpc

#endif//KMLDPC_KMEANS_H

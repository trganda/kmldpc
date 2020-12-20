#ifndef KMLDPC_KMEANS_H
#define KMLDPC_KMEANS_H

#include <vector>
#include <complex>
#include <algorithm>
#include <utility>
#include <iostream>

#include "mat.h"

namespace kmldpc {

    class KMeans {
    public:
        explicit KMeans(std::vector<std::complex<double>> &data,
                        std::vector<std::complex<double>> &constellations,
                        int iter);

        ~KMeans();

    public:
        std::vector<std::complex<double>> GetClusters();

        std::vector<int> GetIdx();

        void Run();

        void DumpToMat(std::string &filename, std::vector<std::complex<double>> &append);

    private:
        const std::vector<std::complex<double>> data_;
        const std::vector<std::complex<double>> constellations_;
        std::vector<std::complex<double>> clusters_;
        std::vector<int> idx_;
        unsigned int iter_;
    };

} // namespace

#endif //KMLDPC_KMEANS_H

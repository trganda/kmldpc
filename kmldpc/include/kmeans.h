#ifndef KMLDPC_KMEANS_H
#define KMLDPC_KMEANS_H

#include <vector>
#include <complex>
#include <algorithm>
#include "Mat.h"

namespace kmldpc {
    class KMeans {
    public:
        explicit KMeans(std::vector<std::complex<double>> &data,
                        std::vector<std::complex<double>> &constellations,
                        int iter);
        ~KMeans();

    public:
        std::vector<std::complex<double>> getClusters();

        std::vector<int> getIdx();

        void run();

        void dumpToMat(std::string& filename, std::vector<std::complex<double>>& append);

    private:
        const std::vector<std::complex<double>> _data;
        const std::vector<std::complex<double>> _constellations;
        std::vector<std::complex<double>> _clusters;
        std::vector<int> _idx;
        unsigned int _iter;
    };
}

#endif //KMLDPC_KMEANS_H

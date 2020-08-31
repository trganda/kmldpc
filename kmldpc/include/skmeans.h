#ifndef KMLDPC_SKMEANS_H
#define KMLDPC_SKMEANS_H

#include <vector>
#include <complex>
#include <algorithm>

class SKMeans {
public:
    explicit SKMeans(std::vector<std::complex<double>> data, std::vector<std::complex<double>> constellations, int iter);
    ~SKMeans() = default;

public:
    std::vector<std::complex<double>> getClusters();
    std::vector<int> getIdx();
    void run();

private:
    const std::vector<std::complex<double>> _data;
    const std::vector<std::complex<double>> _constellations;
    std::vector<std::complex<double>> _clusters;
    std::vector<int> _idx;
    unsigned int _iter;
};


#endif //KMLDPC_SKMEANS_H

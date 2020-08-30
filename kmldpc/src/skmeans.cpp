//
// Created by Jong on 2020/8/29.
//

#include "skmeans.h"

#include <utility>

SKMeans::SKMeans(std::vector<std::complex<double>> data,
                 std::vector<std::complex<double>> constellations, int iter)
                 : _data(std::move(data)), _constellations(std::move(constellations)), _iter(iter) {
    _clusters = std::vector<std::complex<double>>(constellations.size());
    _idx = std::vector<int>(data.size());
}

void SKMeans::run() {
    // Find the point far from the origin to be the outlier
    std::vector<double> absValues(_data.size());
    for (int i=0; i<absValues.size(); i++) {
        absValues[i] = abs(_data[i]);
    }
    auto maxIndex = max_element(absValues.begin(), absValues.end()) - absValues.begin();
    auto maxValue = _data[maxIndex];

    // Initial clusters
    auto hatH = maxValue / _constellations[0];
    for (int i = 0; i < _clusters.size(); i++) {
        _clusters[i] = _constellations[i] * hatH;
    }

    std::vector<std::complex<double>> tempClusters(_clusters.size());
    for (int i = 0; i < _iter; i++) {
        std::vector<int> idxCount(_clusters.size());
        std::vector<std::complex<double>> idxSum(_clusters.size());
        for (int j = 0; j < _data.size(); j++) {
            std::vector<double> tempAbsValues(_clusters.size());
            for (int k = 0; k < _clusters.size(); k++) {
                tempAbsValues[k] = abs(_clusters[k] - _data[j]);
            }
            auto minIndex = min_element(tempAbsValues.begin(), tempAbsValues.end()) - tempAbsValues.begin();
            _idx[j] = minIndex;
            idxCount[_idx[j]]++;
            idxSum[_idx[j]] += _data[j];
        }

        bool flag = true;
        for (int j = 0; j < _clusters.size(); j++) {
            if (_clusters[j] != tempClusters[j]) {
                flag &= false;
                break;
            }
        }

        if (flag) {
            break;
        }

        // Updating the cluster
        tempClusters = _clusters;
        for (int j = 0; j < _clusters.size(); j++) {
            _clusters[j] = idxSum[j] / std::complex<double>(idxCount[j], 0);
        }

        // Form to the constellation schema
        absValues.clear();
        for (int j = 0; j < absValues.size(); j++) {
            absValues[j] = abs(_clusters[j]);
        }
        maxIndex = max_element(absValues.begin(), absValues.end()) - absValues.begin();
        maxValue = _clusters[maxIndex];

        hatH = maxValue / _constellations[0];
        for (int j = 0; j < _clusters.size(); j++) {
            _clusters[j] = _constellations[j] * hatH;
        }
    }

    for (int i = 0; i < _data.size(); i++) {
        std::vector<double> tempAbsValues(_clusters.size());
        for (int j = 0; j < _clusters.size(); j++) {
            tempAbsValues[j] = abs(_clusters[j] - _data[i]);
        }
        auto minIndex = min_element(tempAbsValues.begin(), tempAbsValues.end()) - tempAbsValues.begin();
        _idx[i] = minIndex;
    }
}

std::vector<std::complex<double>> SKMeans::getClusters() {
    return this->_clusters;
}

std::vector<int> SKMeans::getIdx() {
    return this->_idx;
}

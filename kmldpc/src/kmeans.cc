#include "kmeans.h"

namespace kmldpc {
KMeans::KMeans(
    std::vector<std::complex<double>> &data,
    std::vector<std::complex<double>> &constellations, int iter
)
    : data_(data), constellations_(constellations), iter_(iter) {
    clusters_ = std::vector<std::complex<double>>(constellations.size());
    idx_ = std::vector<int>(data.size());
}

KMeans::~KMeans() = default;

void KMeans::Run() {
    // Find the point far from the origin to be the outlier
    std::vector<double> absValues(data_.size());
    for (size_t i = 0; i < absValues.size(); i++) {
        absValues[i] = abs(data_[i]);
    }
    auto maxIndex = std::distance(absValues.begin(), max_element(absValues.begin(), absValues.end()));
    auto maxValue = data_[maxIndex];

    // Initial clusters
    auto hatH = maxValue / constellations_[0];
    for (size_t i = 0; i < clusters_.size(); i++) {
        clusters_[i] = constellations_[i] * hatH;
    }
    std::vector<std::complex<double>> tempClusters(clusters_.size());
    std::vector<int> idxCount(clusters_.size());
    std::vector<std::complex<double>> idxSum(clusters_.size());
    for (size_t i = 0; i < iter_; i++) {
        idxCount.clear();
        idxSum.clear();
        std::vector<double> tempAbsValues(clusters_.size());
        for (size_t j = 0; j < data_.size(); j++) {
            for (size_t k = 0; k < clusters_.size(); k++) {
                tempAbsValues[k] = abs(clusters_[k] - data_[j]);
            }
            auto minIndex = std::distance(
                tempAbsValues.begin(),
                min_element(tempAbsValues.begin(), tempAbsValues.end()));
            idx_[j] = minIndex;
            idxCount[idx_[j]]++;
            idxSum[idx_[j]] += data_[j];
        }
        bool flag = true;
        for (size_t j = 0; j < clusters_.size(); j++) {
            if (clusters_[j] != tempClusters[j]) {
                flag &= false;
                break;
            }
        }
        if (flag) {
            break;
        }

        // Updating the cluster
        tempClusters = clusters_;
        for (size_t j = 0; j < clusters_.size(); j++) {
            clusters_[j] = idxSum[j] / std::complex<double>(idxCount[j], 0);
        }

        // Form to the constellation schema
        absValues.clear();
        for (size_t j = 0; j < absValues.size(); j++) {
            absValues[j] = abs(clusters_[j]);
        }
        maxIndex = std::distance(absValues.begin(), max_element(absValues.begin(), absValues.end()));
        maxValue = clusters_[maxIndex];
        hatH = maxValue / constellations_[0];
        for (size_t j = 0; j < clusters_.size(); j++) {
            clusters_[j] = constellations_[j] * hatH;
        }
    }
    for (size_t i = 0; i < data_.size(); i++) {
        std::vector<double> tempAbsValues(clusters_.size());
        for (size_t j = 0; j < clusters_.size(); j++) {
            tempAbsValues[j] = abs(clusters_[j] - data_[i]);
        }
        auto minIndex = min_element(tempAbsValues.begin(), tempAbsValues.end()) - tempAbsValues.begin();
        idx_[i] = minIndex;
    }
}

std::vector<std::complex<double>> KMeans::clusters() {
    return this->clusters_;
}

std::vector<int> KMeans::idx() {
    return this->idx_;
}

#ifndef NO_MATIO

void KMeans::DumpToMat(std::string &filename, std::vector<std::complex<double>> &append) {
    lab::Mat mat = lab::Mat(filename);
    mat.Open();
    mat.WriteVector("data", data_);
    mat.WriteVector("cluster", clusters_);
    mat.WriteVector("idx", idx_);
    mat.WriteVector("constellations", constellations_);
    mat.WriteVector("hHats", std::vector<std::complex<double>>(append.begin(), append.begin() + 4));
    mat.WriteComplex("realH", append[4]);
    mat.Close();
}

#endif
}   // namespace
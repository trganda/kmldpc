#ifndef LAB_XOR_SEG_CODEC_HPP
#define LAB_XOR_SEG_CODEC_HPP

#include <vector>
#include <complex>
#include <algorithm>
#include <iomanip>

#include "binary5gldpccodec.hpp"
#include "binaryldpccodec.hpp"
#include "modemlinearsystem.hpp"
#include "log.hpp"
#include "mat.hpp"
#include "randnum.hpp"
#include "utility.hpp"

namespace lab {

    class XORSegCodec {
    public:
        explicit XORSegCodec()
                : ldpc_codec_5g_(Binary5GLDPCCodec()), ldpc_codec_(BinaryLDPCCodec()),
                  iter_cnt_(0), using_ldpc_5g_(0), using_syndrom_metric_(0),
                  uu_len_(0), cc_len_(0),
                  rr_(nullptr), bit_l_in_(nullptr), bit_l_out_(nullptr) {}

        virtual ~XORSegCodec() {
            delete[] rr_;
            delete[] bit_l_in_;
            delete[] bit_l_out_;
        }

        void Malloc(int code_no, const char *file_name) {
            // setup from "Setup_of_Codec.txt"
            char temp_str[80] = {' '};
            char mark[80];

            FILE *fp = fopen(file_name, "r");
            if (nullptr == fp) {
                fprintf(stderr, "\nCannot Open %s", file_name);
                exit(-1);
            }

            sprintf(mark, "RanXORLDPC***%d***PARAMETERS", code_no);
            while (strcmp(temp_str, mark)) {
                fscanf(fp, "%s", temp_str);
            }

            fscanf(fp, "%s", temp_str);
            fscanf(fp, "%u", &using_ldpc_5g_);

            fscanf(fp, "%s", temp_str);
            fscanf(fp, "%d", &iter_cnt_);

            fscanf(fp, "%s", temp_str);
            fscanf(fp, "%u", &using_syndrom_metric_);

            fscanf(fp, "%s", temp_str);
            fscanf(fp, "%s", temp_str);

            fclose(fp);

            if (using_ldpc_5g_ == 1) {
                LOG(logger::Info, true) << "Using 5G LDPC." << std::endl;
                ldpc_codec_5g_.Malloc(code_no, temp_str);
                uu_len_ = ldpc_codec_5g_.GetCodeDim();
                cc_len_ = ldpc_codec_5g_.GetCodeenPuncture();
            } else {
                LOG(logger::Info, true) << "Using traditional LDPC." << std::endl;
                ldpc_codec_.Malloc(code_no, temp_str);
                uu_len_ = ldpc_codec_.GetCodeDim();
                cc_len_ = ldpc_codec_.GetCodeLen();
            }

            rr_ = new int[cc_len_];
            bit_l_in_ = new double[cc_len_];
            bit_l_out_ = new double[cc_len_];
        }

        void Encoder(int *uu, int *cc) {
            if (using_ldpc_5g_ == 1) {
                ldpc_codec_5g_.Encoder(uu, cc);
            } else {
                ldpc_codec_.Encoder(uu, cc);
            }
        }

        void Decoder(ModemLinearSystem &modem_linear_system,
                     const std::vector<std::complex<double>> &hHats, int *uu_hat) {
            std::vector<double> metric_results;
            std::vector<std::pair<int, std::complex<double>>> temp;

            metric_results = GetMetrics(modem_linear_system, hHats, uu_hat);

            auto minIndex = std::distance(metric_results.begin(),
                                          min_element(metric_results.begin(), metric_results.end()));
            LOG(logger::Info, false) << "hatIndex = " << minIndex << std::endl;
            temp = {std::pair<int, std::complex<double>>(0, hHats[minIndex])};
            DeMapping(modem_linear_system, temp);

            if (using_ldpc_5g_ == 1) {
                ldpc_codec_5g_.Decoder(bit_l_out_, uu_hat, ldpc_codec_5g_.GetMaxIter());
            } else {
                ldpc_codec_.Decoder(bit_l_out_, uu_hat, ldpc_codec_.GetMaxIter());
            }
        }

        std::vector<double> GetHistogramData(ModemLinearSystem &modem_linear_system,
                                             const std::vector<std::complex<double>> &hHats, int *uu_hat) {
            return GetMetrics(modem_linear_system, hHats, uu_hat);
        }

        // Getter
        int GetUuLen() const {
            return uu_len_;
        }

        int GetCcLen() const {
            return cc_len_;
        }

    private:
        void DeMapping(ModemLinearSystem &modem_linear_system,
                       std::vector<std::pair<int, std::complex<double>>> &thetaList) const {
            modem_linear_system.SoftDemodulation(thetaList);

            // demapping to Get soft information
            for (int i = 0; i < cc_len_; i++) {
                bit_l_in_[i] = 0.5;
            }

            modem_linear_system.GetModem().DeMapping(
                    bit_l_in_, modem_linear_system.GetSymProb(),
                    bit_l_out_, cc_len_);
        }

        int GetParityCheck() const {
            if (using_ldpc_5g_ == 1) {
                return ldpc_codec_5g_.ParityCheck(ldpc_codec_5g_.GetCcHat());
            } else {
                // hard decision without decoding
                for (int i = 0; i < cc_len_; i++) {
                    if (bit_l_out_[i] > 0.5) {
                        rr_[i] = 1;
                    } else {
                        rr_[i] = 0;
                    }
                }
                // For Debug
                std::vector<int> rr_debug(rr_, rr_ + cc_len_);

                return ldpc_codec_.ParityCheck(rr_);
            }
        }

        std::vector<double> GetMetrics(ModemLinearSystem &modem_linear_system,
                                       const std::vector<std::complex<double>> &hHats, int *uu_hat) {
            std::vector<double> metric_results(hHats.size(), 0);
            std::vector<std::pair<int, std::complex<double>>> temp;
            for (auto i = 0; i < metric_results.size(); i++) {
                temp = {std::pair<int, std::complex<double>>(0, hHats[i])};
                DeMapping(modem_linear_system, temp);
                if (using_ldpc_5g_ == 1) {
                    metric_results[i] = Metric(ldpc_codec_5g_, uu_hat);
                } else {
                    metric_results[i] = Metric(ldpc_codec_, uu_hat);
                }

                LOG(lab::logger::Info, false) << std::fixed << std::setprecision(14)
                                              << "Hhat = " << hHats[i]
                                              << " Metric = "
                                              << std::setw(5) << std::right
                                              << metric_results[i] << std::endl;

                metric_results[i] = abs(metric_results[i]);
            }
            return metric_results;
        }

        double Metric(BinaryLDPCCodec &codec, int *uu_hat) {
            double metric_result;

            if (using_syndrom_metric_ == 1) {
                codec.Decoder(bit_l_out_, uu_hat, iter_cnt_);
                std::vector<double> soft_syndroms;
                soft_syndroms = std::vector<double>(codec.GetSyndromSoft(),
                                                    codec.GetSyndromSoft() + codec.GetNumRow());
                for (auto j = 0; j < codec.GetNumRow(); j++) {
                    metric_result += log(soft_syndroms[j]);
                }
            } else {
                if (using_ldpc_5g_ == 1) {
                    codec.Decoder(bit_l_out_, uu_hat, iter_cnt_);
                }
                metric_result = GetParityCheck();
            }
            return metric_result;
        }

    private:
        Binary5GLDPCCodec ldpc_codec_5g_;
        BinaryLDPCCodec ldpc_codec_;

        int iter_cnt_;       // iteration times while using LDPC on 5G
        unsigned int using_ldpc_5g_;
        unsigned int using_syndrom_metric_;
        int uu_len_;         // length of uu for LDPC
        int cc_len_;         // length of cc for LDPC
        int *rr_;             // hard decision
        double *bit_l_in_;     // bit input probability
        double *bit_l_out_;  // bit out probability
    };

}

#endif
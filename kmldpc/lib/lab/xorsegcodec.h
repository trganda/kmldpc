#ifndef LAB_XOR_SEG_CODEC_H
#define LAB_XOR_SEG_CODEC_H

#include <vector>
#include <complex>
#include <algorithm>
#include <iomanip>

#include "binary5gldpccodec.h"
#include "binaryldpccodec.h"
#include "modemlinearsystem.h"
#include "log.h"
#include "mat.h"
#include "randnum.h"
#include "utility.h"

namespace lab {

    class XORSegCodec {
    public:
        explicit XORSegCodec();
        virtual ~XORSegCodec();

        void Malloc(const toml::value& arguments);
        void Encoder(int *uu, int *cc);
        void Decoder(ModemLinearSystem &modem_linear_system,
                     const std::vector<std::complex<double>> &hHats, int *uu_hat);
        std::vector<double> GetHistogramData(ModemLinearSystem &modem_linear_system,
                                             const std::vector<std::complex<double>> &hHats, int *uu_hat);
        // Getter
        int GetUuLen() const;
        int GetCcLen() const;
    private:
        void DeMapping(ModemLinearSystem &modem_linear_system,
                       std::vector<std::pair<int, std::complex<double>>> &thetaList) const;
        int GetParityCheck() const;
        std::vector<double> GetMetrics(ModemLinearSystem &modem_linear_system,
                                       const std::vector<std::complex<double>> &hHats, int *uu_hat);
        double Metric(BinaryLDPCCodec &codec, int *uu_hat);

    private:
        Binary5GLDPCCodec ldpc_codec_5g_;
        BinaryLDPCCodec ldpc_codec_;

        int iter_cnt_;       // iteration times while using LDPC on 5G
        bool using_ldpc_5g_;
        bool using_syndrom_metric_;
        int uu_len_;         // length of uu for LDPC
        int cc_len_;         // length of cc for LDPC
        int *rr_;             // hard decision
        double *bit_l_in_;     // bit input probability
        double *bit_l_out_;  // bit out probability
    };

}

#endif
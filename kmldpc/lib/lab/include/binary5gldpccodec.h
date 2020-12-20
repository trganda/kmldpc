#ifndef LAB_BINARY_5GLDPC_CODEC_H
#define LAB_BINARY_5GLDPC_CODEC_H

#include <cstdio>
#include <cstring>

#include "utility.h"
#include "binaryldpccodec.h"

namespace lab {

    class Binary5GLDPCCodec : public BinaryLDPCCodec {
    public:
        explicit Binary5GLDPCCodec(const toml::value& arguments);
        ~Binary5GLDPCCodec() override;

        void init(const toml::value& arguments) override;

        void Encoder(int *uu, int *cc) const override;

        int Decoder(const double *M2V, int *uu_hat, int iter_cnt) override;

        int GetCodeenPuncture() const;

    private:
        void SystH_5G();

    private:
        int code_len_no_puncture_;  // Code length before puncture
        int code_len_puncture_;     // Code length after puncture
        int lifting_factor_;
        int *cc_no_puncture_;
        double *cc_soft_no_puncture_;
    };

}

#endif //KMLDPC_BINARY5GLDPCCODEC_H

#ifndef KMLDPC_LDPC_LINEAR_SYSTEM_H
#define KMLDPC_LDPC_LINEAR_SYSTEM_H

#include <complex>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>

#include "sourcesink.h"
#include "xorsegcodec.h"
#include "modemlinearsystem.h"
#include "randnum.h"
#include "log.h"
#include "toml.hpp"

#include "kmeans.h"

class LDPCLinearSystem {
public:
    explicit LDPCLinearSystem(toml::value arguments);

    virtual ~LDPCLinearSystem();

    void Simulator();

private:
    const toml::value arguments_;
    // Simulation range of snr
    double min_snr_;
    double max_snr_;
    // Step size for increase snr
    double step_snr_;
    // Maximum error blocks, stop simulation of one snr while
    // errors block is equal to it
    int max_err_blk_;
    // Maximum blocks for simulation
    int max_num_blk_;
    // Uncoded codeword
    int *uu_;
    int *uu_hat_;
    int uu_len_;
    // Encoded codeword
    int *cc_;
    int *cc_hat_;
    int cc_len_;

    double *sym_prob_;

    lab::CSourceSink source_sink_;
    lab::XORSegCodec codec_;
    lab::ModemLinearSystem modem_linear_system_;
};

#endif
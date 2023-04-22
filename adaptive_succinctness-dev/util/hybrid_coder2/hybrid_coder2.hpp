#ifndef _HYBRID_CODER2_HYBRID_CODER2_HPP_
#define _HYBRID_CODER2_HYBRID_CODER2_HPP_

#include <inttypes.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>
#include <sdsl/hyb_vector.hpp>

#include "huffman_coder2.hpp"
#include "tunstallCoder2.hpp"
#include "utils.hpp"

struct Hybrid_coder_test {
    //tunstall_coder tunstall_seq;
    //huffman_coder huffman_seq;
    Tunstall_coder_test tunstall_seq;
    sdsl::sd_vector<> B;
    sdsl::hyb_vector<> hyb_B;
    bool eliasfano_B = true;
    uint64_t top_most_freq;
    uint64_t n;  // original sequence length
};

namespace m_hybrid_coder {
    uint64_t size_bits(Hybrid_coder_test& h);
    uint64_t size_bytes(Hybrid_coder_test& h);
    uint64_t bytes_tunstall_seq(Hybrid_coder_test& h);
    uint64_t bytes_huffman_seq(Hybrid_coder_test& h);
    uint64_t bytes_bit_vec_b(Hybrid_coder_test& h);
    uint64_t bits_tunstall_seq(Hybrid_coder_test& h);
    uint64_t bits_huffman_seq(Hybrid_coder_test& h);
    uint64_t bits_bit_vec_b(Hybrid_coder_test& h);
    void print_info(Hybrid_coder_test& h);

}

#endif
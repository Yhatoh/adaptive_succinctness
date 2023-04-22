#ifndef _HYBRID_CODER_HPP_
#define _HYBRID_CODER_HPP_

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

#include "huffman_coder.hpp"
#include "tunstallCoder.hpp"
#include "utils.hpp"

class hybrid_coder {
    tunstall_coder tunstall_seq;
    huffman_coder huffman_seq;
    sdsl::sd_vector<> B;
    sdsl::hyb_vector<> hyb_B;
    bool eliasfano_B = true;
    uint64_t top_most_freq;
    uint64_t n;  // original sequence length

   public:
    hybrid_coder();
    hybrid_coder(const std::vector<uint32_t> &seq, const uint32_t block_size, 
        const uint64_t top_k, float &symbol_p, bool eliasfano);

    uint64_t size_bits();
    uint64_t size_bytes();
    uint64_t bytes_tunstall_seq();
    uint64_t bytes_huffman_seq();
    uint64_t bytes_bit_vec_b();
    uint64_t bits_tunstall_seq();
    uint64_t bits_huffman_seq();
    uint64_t bits_bit_vec_b();
    void print_info();

    // saving size results of hybrid coder
    void fetch_results(std::vector<uint64_t> &t_seq_bytes,
        std::vector<float> &t_seq_p, 
        std::vector<uint64_t> &h_seq_bytes,
        std::vector<float> &h_seq_p, 
        std::vector<uint64_t> &bit_vec_bytes,
        std::vector<float> &bit_vec_p,
        std::vector<uint64_t> &total_bytes
        );
};

#endif

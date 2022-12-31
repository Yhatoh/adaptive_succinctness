#ifndef _HYBRID_CODER2_HUFFMAN_CODER_HPP_
#define _HYBRID_CODER2_HUFFMAN_CODER_HPP_

#include <vector>
#include <map>
#include <inttypes.h>
#include <math.h>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>

#include "block_element.hpp"
#include "nqsort.h"
#include "mysort.h"


constexpr uint32_t L = 31;
constexpr uint32_t LUT_BITS = 8;
constexpr uint32_t LUT_SIZE = (1 << LUT_BITS);
constexpr uint32_t MAX_ULONG = 0xffffffff;
constexpr uint32_t MAX_IT = 0x00ffffff;
constexpr int32_t BUFF_BITS = sizeof(uint32_t) << 3; //# bits in a buffer element


struct Huffman_coder_test {
    uint32_t cw_lens[L + 1];  // OJO cuidado con el valor de L
    /* Canonical coding arrays */
    uint32_t min_code[L];
    uint32_t lj_base[L]; 
    uint32_t offset[L];
    uint32_t max_cw_len;
    uint32_t min_cw_len; /* used to start the linear search if lut==NULL */
    uint32_t* freq_table;
    uint32_t max_symbol;    
    uint32_t n_symbols;
    uint32_t* syms;  // symbol map
    int32_t* lens;
    uint32_t* lut[LUT_SIZE]; /* canonical decode array */
    uint64_t length;  // length of the original sequence
    std::vector<uint32_t> compressed_seq;  // array with the compressed sequence
    uint32_t buff;  // current uint32_t being decompressed/compressed
    uint32_t buff_btg;
    uint32_t block_size;
    sdsl::enc_vector<> block_info_prefix_sum;
    sdsl::enc_vector<> block_info_starting_position;
};

namespace my_huffman_coder {
    
} // namespace my_huffman_coder


#endif

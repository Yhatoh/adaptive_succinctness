#ifndef _MY_TYPES_HPP_
#define _MY_TYPES_HPP_

#include <sdsl/bit_vectors.hpp>
#include <sdsl/hyb_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <block_element.hpp>

// struct Tunstall_coder_test {
//     std::vector<std::vector<uint64_t>> D;  // Dictionary
//     std::vector<uint32_t> map_table;
//     std::vector<uint16_t> compressed_seq;
//     std::vector<blockElement> block;
//     //enc_vector<> block_info_prefix_sum;
//     //enc_vector<> block_info_starting_position;
//     uint32_t bSize;
//     uint32_t sigma;
//     //uint64_t D_size=65536 by default
//     uint64_t D_size;
// };

// struct Huffman_coder_test {
//     uint32_t cw_lens[L + 1];  // OJO cuidado con el valor de L
//     /* Canonical coding arrays */
//     uint32_t min_code[L];
//     uint32_t lj_base[L]; 
//     uint32_t offset[L];
//     uint32_t max_cw_len;
//     uint32_t min_cw_len; /* used to start the linear search if lut==NULL */
//     uint32_t* freq_table;
//     uint32_t max_symbol;    
//     uint32_t n_symbols;
//     uint32_t* syms;  // symbol map
//     int32_t* lens;
//     uint32_t* lut[LUT_SIZE]; /* canonical decode array */
//     uint64_t length;  // length of the original sequence
//     std::vector<uint32_t> compressed_seq;  // array with the compressed sequence
//     uint32_t buff;  // current uint32_t being decompressed/compressed
//     uint32_t buff_btg;
//     uint32_t block_size;
//     sdsl::enc_vector<> block_info_prefix_sum;
//     sdsl::enc_vector<> block_info_starting_position;
// };



// struct Hybrid_coder_test {
//     //tunstall_coder tunstall_seq;
//     //huffman_coder huffman_seq;
//     sdsl::sd_vector<> B;
//     sdsl::hyb_vector<> hyb_B;
//     bool eliasfano_B = true;
//     uint64_t top_most_freq;
//     uint64_t n;  // original sequence length
// };


#endif
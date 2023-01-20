#ifndef _RUNENCODER_HPP_
#define _RUNENCODER_HPP_

#include <sdsl/bit_vectors.hpp>

#include <vector>
#include <iostream>

#include "tunstallCoder.hpp"
#include "huffman_coder.hpp"

template< uint16_t w >
class RunEncoder {
public:
  //tunstall_coder<w> tc_r0;
  // r0 encoding
  tunstall_coder<w> tc_r0_top_k;
  huffman_coder huffman_r0;
  // r1 encoding
  tunstall_coder<w> tc_r1_top_k;
  huffman_coder huffman_r1;

  sdsl::bit_vector tc_or_huffman;
  //tunstall_coder<w> tc_r1;
  //tunstall_coder<w> tc;
  uint64_t u; // universe
  uint64_t n; // ones
  uint64_t top_most_freq;
  uint64_t symbol_tc_p;
  RunEncoder(sdsl::bit_vector &bv, uint64_t top_k);
  RunEncoder(std::vector<uint64_t> &pb, uint64_t top_k, uint64_t block_size);
  uint64_t bits_tunstall_seq();
private:
  uint64_t get_runs(std::vector<uint64_t> &pb, std::vector<uint32_t> &res);
  void top_k_encoding(std::vector< uint32_t > &pb, int block_size, bool type);
};

#endif

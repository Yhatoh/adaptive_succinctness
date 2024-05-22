#ifndef _RUNENCODERSDARRAY2_HPP_
#define _RUNENCODERSDARRAY2_HPP_

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>

#include "../../s9_vector/include/sdsl/s9_vector.hpp"

#include <algorithm>
#include <vector>
#include <iostream>

#include "../util/tunstallCoder.hpp"
//#include "../util/huffman_coder.hpp"

#include <chrono>

extern std::chrono::high_resolution_clock::time_point block_start, block_stop;
extern std::chrono::high_resolution_clock::time_point tunst_start, tunst_stop;
extern std::chrono::high_resolution_clock::time_point huff_start, huff_stop;
extern std::chrono::high_resolution_clock::time_point select_start, select_stop;
extern std::chrono::duration< double > block_time, huff_time, tunst_time, select_time;
extern double block_total_time;
extern double tunst_total_time;
extern double huff_total_time;
extern double select_total_time;

template< uint16_t w, uint64_t bs, uint64_t br, class gap_class, class _select_gap, class _bv, class _select, class _rank>
class RunEncoderSelect {
public:
  //tunstall_coder<w> tc_r0;
  // r0 encoding
  tunstall_coder<w> tc_r0_top_k;
  gap_class s9_r0;
  _select_gap select_s9_r0;
  vector< uint64_t > dict_s9_r0;

  // r1 encoding
  tunstall_coder<w> tc_r1_top_k;
  gap_class s9_r1;
  _select_gap select_s9_r1;
  vector< uint64_t > dict_s9_r1;

  _bv tc_or_huffman_r0;
  uint64_t n_tchuff_r0;
  _rank rank_tchuff_r0;
  _select select_tchuff_r0;

  _bv tc_or_huffman_r1;
  uint64_t n_tchuff_r1;
  _rank rank_tchuff_r1;
  _select select_tchuff_r1;

  _bv block_r1;
  uint64_t n_block_r1;
  _select select_block_r1;
  _rank rank_block_r1;

  _bv block_r0;
  uint64_t n_block_r0;
  _select select_block_r0;
  _rank rank_block_r0;

  uint64_t u; // universe
  uint64_t n; // ones
  uint64_t n_r0; // ones
  uint64_t n_r1; // ones
  uint64_t top_most_freq;
  double symbol_tc_p;

  RunEncoderSelect(sdsl::bit_vector &bv, uint64_t top_k);
  RunEncoderSelect(std::vector<uint64_t> &pb, uint64_t top_k);
  uint64_t bits_tunstall_seq();
  uint64_t size_block_r1();
  uint64_t size_block_r0();
  uint64_t select(uint64_t k);
  uint64_t rank(uint64_t i);
private:
  void top_k_encoding(std::vector< uint32_t > &pb, bool type);
};

#endif

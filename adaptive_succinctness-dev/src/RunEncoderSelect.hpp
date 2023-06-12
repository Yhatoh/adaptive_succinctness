#ifndef _RUNENCODERSDARRAY_HPP_
#define _RUNENCODERSDARRAY_HPP_

#include <sdsl/bit_vectors.hpp>

#include <algorithm>
#include <vector>
#include <iostream>

#include "../util/tunstallCoder.hpp"
#include "../util/huffman_coder.hpp"

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
class RunEncoderSDArray {
public:
  //tunstall_coder<w> tc_r0;
  // r0 encoding
  tunstall_coder<w> tc_r0_top_k;
  huffman_coder huffman_r0;
  // r1 encoding
  tunstall_coder<w> tc_r1_top_k;
  huffman_coder huffman_r1;

  _bv tc_or_huffman_r0;
  _rank rank_tchuff_r0;
  _select select_tchuff_r0;

  _bv tc_or_huffman_r1;
  _rank rank_tchuff_r1;
  _select select_tchuff_r1;

  _bv block_r1;
  _select select_block_r1;
  _rank rank_block_r1;

  _bv block_r0;
  _select select_block_r0;
  _rank rank_block_r0;

  uint64_t u; // universe
  uint64_t n; // ones
  uint64_t n_r0; // ones
  uint64_t n_r1; // ones
  uint64_t top_most_freq;
  uint64_t symbol_tc_p;

  RunEncoderSDArray(sdsl::bit_vector &bv, uint64_t top_k);
  RunEncoderSDArray(std::vector<uint64_t> &pb, uint64_t top_k);
  uint64_t bits_tunstall_seq();
  uint64_t size_block_r1();
  uint64_t size_block_r0();
  uint64_t select(uint64_t k);
  uint64_t rank(uint64_t i);
private:
  uint64_t get_runs(std::vector<uint64_t> &pb, std::vector<uint32_t> &res);
  void top_k_encoding(std::vector< uint32_t > &pb, bool type);
};

#endif

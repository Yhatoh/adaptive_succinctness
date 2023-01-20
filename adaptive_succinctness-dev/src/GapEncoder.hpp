#ifndef _GAPENCODER_HPP_
#define _GAPENCODER_HPP_

#include <sdsl/bit_vectors.hpp>

#include <vector>
#include <iostream>

#include "tunstallCoder.hpp"

template< uint16_t w >
class GapEncoder {
public:
  //tunstall_coder tc_r0;
  //tunstall_coder tc_r1;
  tunstall_coder<w> tc;
  uint64_t u; // universe
  uint64_t n; // ones
  GapEncoder(sdsl::bit_vector &bv);
  GapEncoder(std::vector<uint64_t> &pb);
  uint64_t bits_tunstall_seq();
private:
  uint64_t get_gaps(std::vector<uint64_t> &pb, std::vector<uint32_t> &res);
};

#endif

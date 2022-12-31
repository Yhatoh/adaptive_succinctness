#include "GapEncoder.hpp"

GapEncoder::GapEncoder(sdsl::bit_vector &bv) {

}

GapEncoder::GapEncoder(std::vector<uint64_t> &pb) {
  std::vector<uint32_t> gaps(pb.size(), 0);
  get_gaps(pb, gaps);
  tc = tunstall_coder(gaps, 512, 65536);
}

void GapEncoder::get_gaps(std::vector<uint64_t> &pb, std::vector<uint32_t> &res) {
  res[0] = pb[0] - 1;
  for(uint64_t i = 1; i < pb.size(); i++) {
    res[i] = pb[i] - pb[i - 1] - 1;
  }
}

uint64_t GapEncoder::bits_tunstall_seq() {
  return 8 * tc.size();
}

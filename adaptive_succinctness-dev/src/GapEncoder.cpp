#include "GapEncoder.hpp"

template< uint16_t w >
GapEncoder<w>::GapEncoder(sdsl::bit_vector &bv) {

}

template< uint16_t w >
GapEncoder<w>::GapEncoder(std::vector<uint64_t> &pb) {
  n = pb.size();
  std::vector<uint32_t> gaps(pb.size(), 0);
  std::cout << "Calculating Gaps...\n";
  std::cout << "Receiving a vector of " << pb.size() << " elements...\n";
  uint64_t sum = get_gaps(pb, gaps);

  std::cout << "Creating Tunstall...\n";
  tc = tunstall_coder<w>(gaps, 1024, 1 << w);
  std::cout << "Done...\n";
}

template< uint16_t w >
uint64_t GapEncoder<w>::get_gaps(std::vector<uint64_t> &pb, std::vector<uint32_t> &res) {
  uint64_t sum = 0;
  res[0] = pb[0] - 1;
  for(uint64_t i = 1; i < pb.size(); i++) {
    res[i] = pb[i] - pb[i - 1] - 1;
    sum += pb[i] - pb[i - 1];
    if(res[i] < 0) std::cout << res[i] << "\n";
  }
  return sum;
}

template< uint16_t w >
uint64_t GapEncoder<w>::bits_tunstall_seq() {
  std::cout << "Dict size: " << tc.dict_size() * 8 << " " << (double) tc.dict_size() * 8 / n << std::endl;
  std::cout << "Compressed seq size: " << tc.compressed_seq_size() * 8 << " " << (double) tc.compressed_seq_size() * 8 / n << std::endl;
  std::cout << "Block size: " << tc.block_vec_size() * 8 << " " << (double) tc.block_vec_size() * 8 / n << std::endl;
  return 8 * tc.size();
}

template class GapEncoder<16>;
template class GapEncoder<18>;
template class GapEncoder<20>;
template class GapEncoder<22>;
template class GapEncoder<24>;

#include "RunEncoder.hpp"

template< uint16_t w, uint64_t bs >
RunEncoder<w,bs>::RunEncoder(sdsl::bit_vector &bv, uint64_t top_k) {

}

template< uint16_t w, uint64_t bs >
RunEncoder<w,bs>::RunEncoder(std::vector<uint64_t> &pb, uint64_t top_k, uint64_t block_size) {
  top_most_freq = top_k;
  std::cerr << "Receiving a vector of " << pb.size() << " elements..." << std::endl;

  std::cerr << "Calculating Runs..." << std::endl;
  
  n = pb.size();
  u = pb[pb.size() - 1] + 1;

  std::cerr << "Creating sets R0 and R1..." << std::endl;
  std::vector< uint32_t > R0;
  std::vector< uint32_t > R1;
  uint64_t last = pb[0];
  if(pb[0] > 0) {
    R0.push_back(pb[0] - 1);
  }

  for(uint64_t i = 1; i < pb.size(); i++) {
    if(pb[i] > pb[i - 1] + 1) {
      R1.push_back(pb[i - 1] - last + 1);
      R0.push_back(pb[i] - pb[i - 1] - 1);
      last = pb[i];
    }
  }

  R1.push_back(pb[pb.size() - 1] - last + 1);

  std::cerr << "Creating bit_vectors..." << std::endl;
  std::vector< uint32_t > PB_R0(R0.size(), 0);

  PB_R0[0] = R0[0];
  for(uint64_t i = 1; i < R0.size(); i++) {
    PB_R0[i] = R0[i] + PB_R0[i - 1];
  }

  R0.clear();

  std::vector< uint32_t > PB_R1(R1.size(), 0);

  PB_R1[0] = R1[0];
  for(uint64_t i = 1; i < R1.size(); i++) {
    PB_R1[i] = R1[i] + PB_R1[i - 1];
  }

  R1.clear();

  n_r0 = PB_R0.size();
  n_r1 = PB_R1.size();

  std::cerr << "Creating gaps representation..." << std::endl; 
  std::vector< uint32_t > GapsR0(PB_R0.size(), 0);
  std::vector< uint32_t > GapsR1(PB_R1.size(), 0);

  std::cerr << "Calculating Gaps R0..." << std::endl;
  //GapsR0[0] = PB_R0[0] - 1;
  GapsR0[0] = PB_R0[0];
  for(uint64_t i = 1; i < GapsR0.size(); i++) {
    //GapsR0[i] = PB_R0[i] - PB_R0[i - 1] - 1;
    GapsR0[i] = PB_R0[i] - PB_R0[i - 1];
  }
  /*
  std::cerr << "Creating IDX_BITS AND IDX_ONES R0..." << std::endl;
  idx_bits_r0 = sdsl::int_vector<>(n_r0 / bs + 2, 0);
  idx_ones_r0 = sdsl::int_vector<>(n_r0 / bs + 2, 0);
  uint64_t s_idbr0 = 1;
  uint64_t s_idor0 = 1;
  uint64_t index = 0;
  while(index < n_r0) {
    uint64_t amount = (index + bs < n_r0 ? bs : n_r0 - index + 1);
    uint64_t next = (index + bs < n_r0 ? index + bs : n_r0);
    idx_bits_r0[s_idbr0++] = PB_R0[next - 1] + 1;
    idx_ones_r0[s_idor0++] = idx_ones_r0[s_idor0 - 1] + amount;
    index = next;
  }
  */
  PB_R0.clear();

  /*
  idx_bits_r0.resize(s_idbr0);
  idx_ones_r0.resize(s_idor0);

  std::cerr << "Creating L2_BITS AND L2_ONES R0..." << std::endl;

  uint64_t s_l2_r0 = s_idbr0;
  l2_bits_r0.resize(s_l2_r0);
  l2_ones_r0.resize(s_l2_r0);

  l2_bits_r0_div = u / s_l2_r0 + (u % s_l2_r0 != 0);
  for(uint64_t i = 0; i < s_l2_r0; i++) {
    auto it = std::upper_bound(idx_bits_r0.begin(), idx_bits_r0.end(), i * l2_bits_r0_div);
    l2_bits_r0[i] - std::distance(idx_bits_r0.begin(), it);
  }

  l2_ones_r0_div = (n + 1) / s_l2_r0 + ((n + 1) % s_l2_r0 != 0);
  for(uint64_t i = 0; i < s_l2_r0; i++) {
    auto it = std::upper_bound(idx_ones_r0.begin(), idx_ones_r0.end(), i * l2_ones_r0_div);
    l2_ones_r0[i] - std::distance(idx_ones_r0.begin(), it);
  }
  sdsl::util::bit_compress(idx_bits_r0);
  sdsl::util::bit_compress(idx_ones_r0);
  sdsl::util::bit_compress(l2_bits_r0);
  sdsl::util::bit_compress(l2_ones_r0);
  */

  std::cerr << "Calculating Gaps R1..." << std::endl;
  //GapsR1[0] = PB_R1[0] - 1;
  GapsR1[0] = PB_R1[0];
  for(uint64_t i = 1; i < GapsR1.size(); i++) {
    //GapsR1[i] = PB_R1[i] - PB_R1[i - 1] - 1;
    GapsR1[i] = PB_R1[i] - PB_R1[i - 1];
  }
  /*
  std::cerr << "Creating IDX_BITS AND IDX_ONES R1..." << std::endl;
  idx_bits_r1 = sdsl::int_vector<>(n_r1 / bs + 2, 0);
  idx_ones_r1 = sdsl::int_vector<>(n_r1 / bs + 2, 0);

  uint64_t s_idbr1 = 1;
  uint64_t s_idor1 = 1;
  index = 0;
  while(index < n_r1) {
    uint64_t amount = (index + bs < n_r1 ? bs : n_r1 - index + 1);
    uint64_t next = (index + bs < n_r1 ? index + bs : n_r1);
    idx_bits_r1[s_idbr1++] = PB_R1[next - 1] + 1;
    idx_ones_r1[s_idor1++] = idx_ones_r1[s_idor1 - 1] + amount;
    index = next;
  }

  idx_bits_r1.resize(s_idbr1);
  idx_ones_r1.resize(s_idor1);
  */
  /*
  std::cerr << "Creating L2_BITS AND L2_ONES R0..." << std::endl;
  uint64_t s_l2_r1 = s_idbr1;
  l2_bits_r1.resize(s_l2_r1);
  l2_ones_r1.resize(s_l2_r1);

  l2_bits_r1_div = u / s_l2_r1 + (u % s_l2_r1 != 0);
  for(uint64_t i = 0; i < s_l2_r1; i++) {
    auto it = std::upper_bound(idx_bits_r1.begin(), idx_bits_r1.end(), i * l2_bits_r1_div);
    l2_bits_r1[i] - std::distance(idx_bits_r1.begin(), it);
  }

  l2_ones_r1_div = (n + 1) / s_l2_r1 + ((n + 1) % s_l2_r1 != 0);
  for(uint64_t i = 0; i < s_l2_r1; i++) {
    auto it = std::upper_bound(idx_ones_r1.begin(), idx_ones_r1.end(), i * l2_ones_r1_div);
    l2_ones_r1[i] - std::distance(idx_ones_r1.begin(), it);
  }

  sdsl::util::bit_compress(idx_bits_r1);
  sdsl::util::bit_compress(idx_ones_r1);
  sdsl::util::bit_compress(l2_bits_r1);
  sdsl::util::bit_compress(l2_ones_r1);
  */

  PB_R1.clear();

  //std::cerr << "Creating Tunstall..." << std::endl;
  //tc_r1 = tunstall_coder<w>(GapsR0, 512, 1 << w);
  std::cerr << "Creating Top-k Encoding..." << std::endl;
  top_k_encoding(GapsR0, block_size, false);
  top_k_encoding(GapsR1, block_size, true);

  //tc_r0 = tunstall_coder<w>(GapsR1, 512, 1 << w);
  //tc = tunstall_coder<w>(gaps, 512, 1 << w);
  //std::cerr << "Done..." << std::endl;
}

template< uint16_t w, uint64_t bs >
uint64_t RunEncoder<w,bs>::get_runs(std::vector<uint64_t> &pb, std::vector<uint32_t> &res) {
  uint64_t sum = 0;
  res[0] = pb[0] - 1;
  for(uint64_t i = 1; i < pb.size(); i++) {
    res[i] = pb[i] - pb[i - 1] - 1;
    sum += pb[i] - pb[i - 1];
    if(res[i] < 0) std::cout << res[i] << "\n";
  }
  return sum;
}

template< uint16_t w, uint64_t bs >
void RunEncoder<w,bs>::top_k_encoding(std::vector< uint32_t > &seq, int block_size, bool type) {
  //std::cerr << "Creating frequency map..." << std::endl;
  std::map< uint32_t, uint64_t > freq_map;
  for(uint64_t i = 0; i < seq.size(); i++) {
    if(freq_map.count(seq[i]) == 0) 
      freq_map[seq[i]] = 1;
    else
      freq_map[seq[i]] += 1;
  }
  //std::cerr << "Done..." << std::endl;
  //std::cerr << "Creating priority queue..." << std::endl;

  std::priority_queue< std::pair< uint32_t, uint64_t > > pq;
  for(std::map< uint32_t, uint64_t >::iterator it = freq_map.begin(); it != freq_map.end(); it++) {
    pq.push(std::pair< uint64_t, uint32_t >(it->second, it->first));
  }

  //std::cerr << "Done..." << std::endl;
  //std::cerr << "Tunstall alphabet..." << std::endl;

  std::set< uint32_t > tc_alphabet;
  for(uint64_t i = 0; i < top_most_freq; i++) {
    tc_alphabet.insert(pq.top().second);
    pq.pop();
  }

  //std::cerr << "Done..." << std::endl;
  //std::cerr << "Filling bitvector..." << std::endl;

  tc_or_huffman = sdsl::bit_vector(seq.size(), 0);
  std::vector< uint32_t > tc_seq;
  std::vector< uint32_t > hf_seq;

  for(uint64_t i = 0; i < seq.size(); i++) {
    if(tc_alphabet.count(seq[i]) == 1) {
      tc_seq.push_back(seq[i]);
      tc_or_huffman[i] = 1;
    } else {
      hf_seq.push_back(seq[i]);
    }
  }
 
  //std::cerr << "Done..." << std::endl;
  //std::cerr << "Creating tunstall and huffman..." << std::endl;
  std::cout << "Percentage of symbols in Tunstall sequence: " << (double)tc_seq.size() / seq.size() << std::endl;
  symbol_tc_p = (double)tc_seq.size() / seq.size() ;

  // r1
  if(type) {
    tc_r1_top_k = tunstall_coder<w>(tc_seq, block_size, 1 << w);
    huffman_r1.encode(hf_seq, block_size);
  // r0
  } else {
    tc_r0_top_k = tunstall_coder<w>(tc_seq, block_size, 1 << w); 
    huffman_r0.encode(hf_seq, block_size);
  }
  //std::cerr << "Done..." << std::endl;
}

template< uint16_t w, uint64_t bs >
uint64_t RunEncoder<w,bs>::bits_tunstall_seq() {
  /*std::cout << "--- Gaps R0 ---" << std::endl;
  std::cout << "Dict size: " << tc_r0.dict_size() * 8 << " " << (double) tc_r0.dict_size() * 8 / n << std::endl;
  std::cout << "Compressed seq size: " << tc_r0.compressed_seq_size() * 8 << " " << (double) tc_r0.compressed_seq_size() * 8 / n << std::endl;
  std::cout << "Block size: " << tc_r0.block_vec_size() * 8 << " " << (double) tc_r0.block_vec_size() * 8 / n << std::endl;
  std::cout << "--- Gaps R1 ---" << std::endl;
  std::cout << "Dict size: " << tc_r1.dict_size() * 8 << " " << (double) tc_r1.dict_size() * 8 / n << std::endl;
  std::cout << "Compressed seq size: " << tc_r1.compressed_seq_size() * 8 << " " << (double) tc_r1.compressed_seq_size() * 8 / n << std::endl;
  std::cout << "Block size: " << tc_r1.block_vec_size() * 8 << " " << (double) tc_r1.block_vec_size() * 8 / n << std::endl;
*/
  std::cout << "TC_R0: " << tc_r0_top_k.size() * 8 << std::endl;
  std::cout << "HF_R0: " << huffman_r0.size() * 8 << std::endl;
  std::cout << "TC_R1: " << tc_r1_top_k.size() * 8 << std::endl;
  std::cout << "HF_R1: " << huffman_r1.size() * 8 << std::endl;
  std::cout << "IDX: "<< bits_idx() << std::endl;
  std::cout << "L2: "<< bits_l2() << std::endl;
  return 8 * sdsl::size_in_bytes(tc_or_huffman) + 
         8 * tc_r0_top_k.size() + 8 * tc_r1_top_k.size() +
         8 * huffman_r0.size() + 8 * huffman_r1.size() +
         bits_idx() + bits_l2();
  //return 8 * tc.size();
}

template< uint16_t w, uint64_t bs >
uint64_t RunEncoder<w,bs>::bits_idx() {
  return 8 * sdsl::size_in_bytes(idx_bits_r0) + 8 * sdsl::size_in_bytes(idx_ones_r0) +
         8 * sdsl::size_in_bytes(idx_bits_r1) + 8 * sdsl::size_in_bytes(idx_ones_r1);
}

template< uint16_t w, uint64_t bs >
uint64_t RunEncoder<w,bs>::bits_l2() {
  return 8 * sdsl::size_in_bytes(l2_bits_r0) + 8 * sdsl::size_in_bytes(l2_ones_r0) +
         8 * sdsl::size_in_bytes(l2_bits_r1) + 8 * sdsl::size_in_bytes(l2_ones_r1);
}

template< uint16_t w, uint64_t bs >
uint64_t RunEncoder<w,bs>::select(uint64_t k) {
  uint64_t pos =
}

template< uint16_t w, uint64_t bs >
uint64_t RunEncoder<w,bs>::rank(uint64_t i) {

}

template class RunEncoder<16, 256>;
template class RunEncoder<16, 512>;
template class RunEncoder<16, 1024>;
template class RunEncoder<18, 1024>;
template class RunEncoder<20, 1024>;
template class RunEncoder<22, 1024>;
template class RunEncoder<24, 1024>;

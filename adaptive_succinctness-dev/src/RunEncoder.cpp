#include "RunEncoder.hpp"

template< uint16_t w >
RunEncoder<w>::RunEncoder(sdsl::bit_vector &bv, uint64_t top_k) {

}

template< uint16_t w >
RunEncoder<w>::RunEncoder(std::vector<uint64_t> &pb, uint64_t top_k, uint64_t block_size) {
  top_most_freq = top_k;
  //std::cerr << "Receiving a vector of " << pb.size() << " elements..." << std::endl;

  //std::cerr << "Calculating Runs..." << std::endl;
  
  n = pb.size();
  //std::cerr << "Creating sets R0 and R1..." << std::endl;
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
  //std::cerr << R0.size() << std::endl;
  //std::cerr << R1.size() << std::endl;

  //std::cerr << "Done..." << std::endl;

  //std::cerr << "Creating bit_vectors..." << std::endl;
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

  //std::cerr << "Done..." << std::endl;

  //std::cerr << "Creating gaps representation..." << std::endl; 

  std::vector< uint32_t > GapsR0(PB_R0.size(), 0);
  std::vector< uint32_t > GapsR1(PB_R1.size(), 0);

  GapsR0[0] = PB_R0[0] - 1;
  for(uint64_t i = 1; i < GapsR0.size(); i++) {
    GapsR0[i] = PB_R0[i] - PB_R0[i - 1] - 1;
  }
  
  PB_R0.clear();

  GapsR1[0] = PB_R1[0] - 1;
  for(uint64_t i = 1; i < GapsR1.size(); i++) {
    GapsR1[i] = PB_R1[i] - PB_R1[i - 1] - 1;
  }

  PB_R1.clear();

  //std::cerr << "Done..." << std::endl;
  //std::cerr << "Creating Tunstall..." << std::endl;
  //tc_r1 = tunstall_coder<w>(GapsR0, 512, 1 << w);
  top_k_encoding(GapsR0, block_size, false);
  top_k_encoding(GapsR1, block_size, true);

  //tc_r0 = tunstall_coder<w>(GapsR1, 512, 1 << w);
  //tc = tunstall_coder<w>(gaps, 512, 1 << w);
  //std::cerr << "Done..." << std::endl;
}

template< uint16_t w >
uint64_t RunEncoder<w>::get_runs(std::vector<uint64_t> &pb, std::vector<uint32_t> &res) {
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
void RunEncoder<w>::top_k_encoding(std::vector< uint32_t > &seq, int block_size, bool type) {
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

template< uint16_t w >
uint64_t RunEncoder<w>::bits_tunstall_seq() {
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
  return 8 * sdsl::size_in_bytes(tc_or_huffman) + 8 * tc_r0_top_k.size() + 8 * tc_r1_top_k.size() + huffman_r0.size() * 8 + huffman_r1.size() * 8;
  //return 8 * tc.size();
}

template class RunEncoder<16>;
template class RunEncoder<18>;
template class RunEncoder<20>;
template class RunEncoder<22>;
template class RunEncoder<24>;

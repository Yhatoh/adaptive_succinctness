#include "RunEncoderSDArray.hpp"

#define debug(x) std::cout << #x << ": " << x << std::endl

template< uint16_t w, uint64_t bs, uint64_t br>
RunEncoderSDArray<w,bs,br>::RunEncoderSDArray(sdsl::bit_vector &bv, uint64_t top_k) {

}

template< uint16_t w, uint64_t bs, uint64_t br >
RunEncoderSDArray<w,bs,br>::RunEncoderSDArray(std::vector<uint64_t> &pb, uint64_t top_k) {
  top_most_freq = top_k;
  std::cerr << "Receiving a vector of " << pb.size() << " elements..." << std::endl;

  std::cerr << "Adding a 0 at the beginning" << std::endl;
  
  for(uint64_t i = 0; i < pb.size(); i++) pb[i]++;

  std::cerr << "Calculating Runs..." << std::endl;
  
  n = pb.size();
  u = pb[pb.size() - 1] + 1;

  std::cerr << "Creating sets R0 and R1..." << std::endl;
  std::vector< uint32_t > R0;
  std::vector< uint32_t > R1;

  uint64_t last = pb[0];
  R0.push_back(pb[0]);

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
  //std::cerr << "R0: " << R0[0] << " ";
  for(uint64_t i = 1; i < R0.size(); i++) {
    //std::cerr << R0[i] << " ";
    PB_R0[i] = R0[i] + PB_R0[i - 1];
  }

  R0.clear();

  std::vector< uint32_t > PB_R1(R1.size(), 0);

  PB_R1[0] = R1[0];
  //std::cerr << "R1: " << R1[0] << " ";
  for(uint64_t i = 1; i < R1.size(); i++) {
    //std::cerr << R1[i] << " ";
    PB_R1[i] = R1[i] + PB_R1[i - 1];
  }

  R1.clear();

  n_r0 = PB_R0.size();
  n_r1 = PB_R1.size();

  std::cerr << "N R0: " << n_r0 << std::endl;
  std::cerr << "N R1: " << n_r1 << std::endl;

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

  std::cerr << "Creatings blocks of R0..." << std::endl;
  std::vector< uint64_t > br_R0;
  uint64_t i; 
  for(i = br; i <= PB_R0.size(); i += br) {
    br_R0.push_back(PB_R0[i - 1]);
  }
  
  if(PB_R0.size() % br != 0)
    br_R0.push_back(PB_R0[PB_R0.size() - 1]);
 
  //block_r0.resize(br_R0[br_R0.size() - 1] + 1);
  //sdsl::util::_set_zero_bits(block_r0);
  //for(uint64_t j = 0; j < br_R0.size(); j++) {
  //  block_r0[br_R0[j]] = 1;
  //}
  block_r0 = sdsl::sd_vector<>(br_R0.begin(), br_R0.end());

  br_R0.clear();
  PB_R0.clear();

  std::cerr << "Calculating Gaps R1..." << std::endl;
  //GapsR1[0] = PB_R1[0] - 1;
  GapsR1[0] = PB_R1[0];
  for(uint64_t i = 1; i < GapsR1.size(); i++) {
    //GapsR1[i] = PB_R1[i] - PB_R1[i - 1] - 1;
    GapsR1[i] = PB_R1[i] - PB_R1[i - 1];
  } 


  std::cerr << "Creatings blocks of R1..." << std::endl;
  std::vector< uint64_t > br_R1;
  for(i = br; i <= PB_R1.size(); i += br) {
    br_R1.push_back(PB_R1[i - 1]);
  }
  
  if(PB_R1.size() % br != 0)
    br_R1.push_back(PB_R1[PB_R1.size() - 1]);
 
  //block_r1.resize(br_R1[br_R1.size() - 1] + 1);
  //sdsl::util::_set_zero_bits(block_r1);
  //for(auto j : br_R1) {
  //  block_r1[j] = 1;
  //}
  
  block_r1 = sdsl::sd_vector<>(br_R1.begin(), br_R1.end());

  br_R1.clear();
  PB_R1.clear();

  sdsl::util::init_support(select_block_r0, &block_r0);
  sdsl::util::init_support(rank_block_r0, &block_r0);
  sdsl::util::init_support(select_block_r1, &block_r1);
  sdsl::util::init_support(rank_block_r1, &block_r1);

  //std::cerr << "Creating Tunstall..." << std::endl;
  std::cerr << "Creating Top-k Encoding..." << std::endl;
  top_k_encoding(GapsR0, false);
  top_k_encoding(GapsR1, true);

  sdsl::util::init_support(rank_tchuff_r0, &tc_or_huffman_r0);
  sdsl::util::init_support(rank_tchuff_r1, &tc_or_huffman_r1);
  sdsl::util::init_support(select_tchuff_r0, &tc_or_huffman_r0);
  sdsl::util::init_support(select_tchuff_r1, &tc_or_huffman_r1);

  std::cerr << "Done..." << std::endl;
}

template< uint16_t w, uint64_t bs, uint64_t br >
uint64_t RunEncoderSDArray<w,bs,br>::get_runs(std::vector<uint64_t> &pb, std::vector<uint32_t> &res) {
  uint64_t sum = 0;
  res[0] = pb[0] - 1;
  for(uint64_t i = 1; i < pb.size(); i++) {
    res[i] = pb[i] - pb[i - 1] - 1;
    sum += pb[i] - pb[i - 1];
    if(res[i] < 0) std::cout << res[i] << "\n";
  }
  return sum;
}

template< uint16_t w, uint64_t bs, uint64_t br >
void RunEncoderSDArray<w,bs,br>::top_k_encoding(std::vector< uint32_t > &seq, bool type) {
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

  std::vector< uint64_t > tc_or_huff_seq;

  std::vector< uint32_t > tc_seq;
  std::vector< uint32_t > hf_seq;

  for(uint64_t i = 0; i < seq.size(); i++) {
    if(tc_alphabet.count(seq[i]) == 1) {
      tc_seq.push_back(seq[i]);
      tc_or_huff_seq.push_back(i);
    } else {
      hf_seq.push_back(seq[i]);
    }
  }
 
  if(type) {
    tc_or_huffman_r1 = sdsl::sd_vector<>(tc_or_huff_seq.begin(), tc_or_huff_seq.end());
  } else {
    tc_or_huffman_r0 = sdsl::sd_vector<>(tc_or_huff_seq.begin(), tc_or_huff_seq.end());
  }

  //std::cerr << "Done..." << std::endl;
  //std::cerr << "Creating tunstall and huffman..." << std::endl;
  std::cerr << "Percentage of symbols in Tunstall sequence: " << (double)tc_seq.size() / seq.size() << std::endl;
  std::cerr << tc_seq.size() << " " << hf_seq.size() << "\n";
  symbol_tc_p = (double)tc_seq.size() / seq.size() ;

  // r1
  if(type) {
    tc_r1_top_k = tunstall_coder<w>(tc_seq, bs, 1 << w);
    huffman_r1.encode(hf_seq, bs);
  // r0
  } else {
    tc_r0_top_k = tunstall_coder<w>(tc_seq, bs, 1 << w); 
    huffman_r0.encode(hf_seq, bs);
  }
  //std::cerr << "Done..." << std::endl;
}

template< uint16_t w, uint64_t bs, uint64_t br >
uint64_t RunEncoderSDArray<w,bs,br>::bits_tunstall_seq() {
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
  return 8 * sdsl::size_in_bytes(tc_or_huffman_r1) + 
         8 * sdsl::size_in_bytes(tc_or_huffman_r0) +
         8 * tc_r0_top_k.size() + 8 * tc_r1_top_k.size() +
         8 * huffman_r0.size() + 8 * huffman_r1.size() +
         8 * size_block_r0() + 8 * size_block_r1();
}

template< uint16_t w, uint64_t bs, uint64_t br >
uint64_t RunEncoderSDArray<w,bs,br>::size_block_r1() {
  return sdsl::size_in_bytes(block_r1) + sdsl::size_in_bytes(select_block_r1) + sdsl::size_in_bytes(rank_block_r1);
}

template< uint16_t w, uint64_t bs, uint64_t br >
uint64_t RunEncoderSDArray<w,bs,br>::size_block_r0() {
  return sdsl::size_in_bytes(block_r0) + sdsl::size_in_bytes(select_block_r0) + sdsl::size_in_bytes(rank_block_r0);
}

template< uint16_t w, uint64_t bs, uint64_t br >
uint64_t RunEncoderSDArray<w,bs,br>::select(uint64_t k) {
  assert(k <= n);

  uint64_t block = rank_block_r1(k);

  uint64_t pos = 0;
  bool take_gr0 = true;

  // ones and zeros until block
  uint64_t ones = 0;
  uint64_t zeros = 0;

  // actual gap
  uint64_t gap_pos = 0;

  // ith one R0
  uint64_t one_r0 = 1;
  // ith one R1
  uint64_t one_r1 = 1;

  // current gap in tunstall R0
  uint64_t gap_tc_r0 = 0;
  // current gap in huffman R0
  uint64_t gap_huff_r0 = 0;
  // current gap in tunstall R1
  uint64_t gap_tc_r1 = 0;
  // current gap in huffman R1
  uint64_t gap_huff_r1 = 0;
  
  if(block > 0) {
    // ones and zeros until block
    ones = select_block_r1(block);
    zeros = select_block_r0(block);

    // curr gap pos
    gap_pos = block * br;

    // actual gaps in each structure
    gap_tc_r0 = rank_tchuff_r0(gap_pos);
    one_r0 = gap_tc_r0;

    gap_huff_r0 = gap_pos - gap_tc_r0;
    gap_tc_r1 = rank_tchuff_r1(gap_pos);
    one_r1 = gap_tc_r1++;
    gap_huff_r1 = gap_pos  - gap_tc_r1;

    // pos before block
    pos += zeros;
    pos += ones;
  }

  uint64_t prev_one_r0 = 0;
  bool flag_acum_r0 = false;
  uint64_t prev_one_r1 = 0;
  bool flag_acum_r1 = false;

  uint64_t res_select_r1 = select_tchuff_r1(one_r1);
  //std::cout << "BEFORE WHILE" << std::endl;
  //debug(one_r1);
  //debug(gap_pos);
  //debug(res_select_r1);
  if(one_r1 > 1) {
    prev_one_r1 = res_select_r1;
    if(res_select_r1 > gap_pos) {
      prev_one_r1 = gap_pos;
      flag_acum_r0 = true;
    }
  }
  //debug(prev_one_r1);
  //std::cout << "IN WHILE" << std::endl;
  while(ones < k) {
    if(take_gr0) {
      // read from gap of run 0
      uint64_t act_zeros = 0;
      if(tc_or_huffman_r0[gap_pos]) {
        // is in tunstall
        if(gap_tc_r0 == 0) act_zeros = tc_r0_top_k.decode(gap_tc_r0);
        else act_zeros = tc_r0_top_k.decode(gap_tc_r0) - tc_r0_top_k.decode(gap_tc_r0 - 1);
        gap_tc_r0++;
      } else {
        // is in huffman
        if(gap_huff_r0 == 0) act_zeros = huffman_r0.decode(gap_huff_r0);
        else act_zeros = huffman_r0.decode(gap_huff_r0) - huffman_r0.decode(gap_huff_r0 - 1);
        gap_huff_r0++;
      }
      pos += act_zeros;
    } else {
      // read from gap of run 1
      uint64_t act_ones = 0;
      //if(tc_or_huffman_r1[gap_pos]) {
      if(!flag_acum_r1) res_select_r1 = select_tchuff_r1(one_r1);
      //debug(tc_or_huffman_r1[gap_pos]);
      //debug(res_select_r1);
      //debug(one_r1);
      //debug(prev_one_r1);
      //debug(flag_acum_r1);
      //std::cout << "CYCLE " << std::endl;
      if(res_select_r1 <= prev_one_r1 + 1 && !flag_acum_r1) {
        // is in tunstall
        if(gap_tc_r1 == 0) act_ones = tc_r1_top_k.decode(gap_tc_r1);
        else act_ones = tc_r1_top_k.decode(gap_tc_r1) - tc_r1_top_k.decode(gap_tc_r1 - 1);
        gap_tc_r1++;
        one_r1++;
        prev_one_r1 = res_select_r1;
      } else {
        flag_acum_r1 = true;
        // is in huffman
        if(gap_huff_r1 == 0) act_ones = huffman_r1.decode(gap_huff_r1);
        else act_ones = huffman_r1.decode(gap_huff_r1) - huffman_r1.decode(gap_huff_r1 - 1);
        gap_huff_r1++;
        prev_one_r1++;
        if(res_select_r1 <= prev_one_r1 + 1) {
          flag_acum_r1 = false;
        }
      }
      pos += act_ones;
      ones += act_ones;
      gap_pos++;
    }
    take_gr0 = !take_gr0;
  }
  while(ones > k) {
    pos--;
    ones--;
  }
  return pos - 2;
}

template< uint16_t w, uint64_t bs, uint64_t br >
uint64_t RunEncoderSDArray<w,bs,br>::rank(uint64_t i) {
  i++;
  if(i >= u) return n;
  uint64_t l = 1;
  uint64_t r = rank_block_r1(block_r1.size());
 
  uint64_t block = 0;
  uint64_t _ones = 0;
  uint64_t _zeros = 0;
  uint64_t pos = 0;

  // binary search, find block
  while(l < r) {
    uint64_t mid = (l + r) / 2;

    uint64_t ones_mid = select_block_r1(mid);
    uint64_t zeros_mid = select_block_r0(mid);
    uint64_t pos_mid = ones_mid + zeros_mid;
    if(pos_mid >= i) {
      r = mid - 1;
    } else {
      block = mid;
      _ones = ones_mid;
      _zeros = zeros_mid;
      pos = _ones + _zeros;
      l = mid + 1;
    }
  }

  // actual gap
  uint64_t gap_pos = 0;

  // current gap in tunstall R0
  uint64_t gap_tc_r0 = 0;
  // current gap in huffman R0
  uint64_t gap_huff_r0 = 0;
  // current gap in tunstall R1
  uint64_t gap_tc_r1 = 0;
  // current gap in huffman R1
  uint64_t gap_huff_r1 = 0;

  bool take_gr0 = true;
  if(block > 0 && pos < i) {
    // curr gap pos
    gap_pos = block * br;

    // actual gaps in each structure
    gap_tc_r0 = rank_tchuff_r0(gap_pos);
    gap_huff_r0 = gap_pos - gap_tc_r0;
    gap_tc_r1 = rank_tchuff_r1(gap_pos);
    gap_huff_r1 = gap_pos - gap_tc_r1;
  } else if(block == 1 && pos >= i) {
    _ones = 0;
    _zeros = 0;
    pos = 0;
  }

  while(pos <= i) {
    if(take_gr0) {
      // read from gap of run 0
      uint64_t act_zeros = 0;
      if(tc_or_huffman_r0[gap_pos]) {
        // is in tunstall
        if(gap_tc_r0 == 0) act_zeros = tc_r0_top_k.decode(gap_tc_r0);
        else act_zeros = tc_r0_top_k.decode(gap_tc_r0) - tc_r0_top_k.decode(gap_tc_r0 - 1);
        gap_tc_r0++;
      } else {
        // is in huffman
        if(gap_huff_r0 == 0) act_zeros = huffman_r0.decode(gap_huff_r0);
        else act_zeros = huffman_r0.decode(gap_huff_r0) - huffman_r0.decode(gap_huff_r0 - 1);
        gap_huff_r0++;
      }
      pos += act_zeros;
    } else {
      // read from gap of run 1
      uint64_t act_ones = 0;
      if(tc_or_huffman_r1[gap_pos]) {
        // is in tunstall
        if(gap_tc_r1 == 0) act_ones = tc_r1_top_k.decode(gap_tc_r1);
        else act_ones = tc_r1_top_k.decode(gap_tc_r1) - tc_r1_top_k.decode(gap_tc_r1 - 1);
        gap_tc_r1++;
      } else {
        // is in huffman
        if(gap_huff_r1 == 0) act_ones = huffman_r1.decode(gap_huff_r1);
        else act_ones = huffman_r1.decode(gap_huff_r1) - huffman_r1.decode(gap_huff_r1 - 1);
        gap_huff_r1++;
      }
      pos += act_ones;
      _ones += act_ones;
      gap_pos++;
      if(pos >= i) {
        _ones -= act_ones;
        pos -= act_ones;
        uint64_t to_sum = i - pos;
        pos = i + 1;
        _ones += to_sum;
      }
    }
    take_gr0 = !take_gr0;
  }
  return _ones;
}

template class RunEncoderSDArray<16, 256, 512>;
template class RunEncoderSDArray<16, 512, 512>;
template class RunEncoderSDArray<16, 1024, 512>;
//template class RunEncoderSDArray<18, 1024>;
//template class RunEncoderSDArray<20, 1024>;
//template class RunEncoderSDArray<22, 1024>;
//template class RunEncoderSDArray<24, 1024>;

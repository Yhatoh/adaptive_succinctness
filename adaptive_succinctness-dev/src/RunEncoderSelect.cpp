#include "RunEncoderSelect.hpp"

#include <chrono>
using namespace std;

chrono::high_resolution_clock::time_point block_start, block_stop;
chrono::high_resolution_clock::time_point tunst_start, tunst_stop;
chrono::high_resolution_clock::time_point huff_start, huff_stop;
chrono::high_resolution_clock::time_point select_start, select_stop;
chrono::duration< double > block_time, huff_time, tunst_time, select_time;
double block_total_time = 0;
double tunst_total_time = 0;
double huff_total_time = 0;
double select_total_time = 0;

#define debug(x) std::cout << #x << ": " << x << std::endl

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
RunEncoderSelect<w,bs,br,_bv,_select,_rank>::RunEncoderSelect(sdsl::bit_vector &bv, uint64_t top_k) {

}

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
RunEncoderSelect<w,bs,br,_bv,_select,_rank>::RunEncoderSelect(std::vector<uint64_t> &pb, uint64_t top_k) {
  top_most_freq = top_k;
  std::cerr << "Receiving a vector of " << pb.size() << " elements..." << std::endl;

  std::cerr << "Calculating Runs..." << std::endl;
  
  n = pb.size();
  u = pb[pb.size() - 1] + 1;

  std::cerr << "Creating sets R0 and R1..." << std::endl;
  std::vector< uint64_t > R0;
  std::vector< uint64_t > R1;

  uint64_t last = pb[0] + 1;
  R0.push_back(pb[0] + 1);

  for(uint64_t i = 1; i < pb.size(); i++) {
    if(pb[i] > pb[i - 1] + 1) {
      R1.push_back(pb[i - 1] + 1 - last + 1);
      R0.push_back(pb[i] - pb[i - 1] - 1);
      last = pb[i] + 1;
    }
  }

  R1.push_back(pb[pb.size() - 1] + 1 - last + 1);

  std::cerr << "Creating bit_vectors..." << std::endl;
  std::vector< uint64_t > PB_R0(R0.size(), 0);

  PB_R0[0] = R0[0];
  //std::cerr << "R0: " << R0[0] << " ";
  for(uint64_t i = 1; i < R0.size(); i++) {
    //std::cerr << R0[i] << " ";
    PB_R0[i] = R0[i] + PB_R0[i - 1];
  }

  R0.clear();

  std::vector< uint64_t > PB_R1(R1.size(), 0);

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
  uint32_t maxi = PB_R0[0];
  for(uint64_t i = 1; i < GapsR0.size(); i++) {
    //GapsR0[i] = PB_R0[i] - PB_R0[i - 1] - 1;
    GapsR0[i] = PB_R0[i] - PB_R0[i - 1];
    maxi = max(GapsR0[i], maxi);
  }

  std::cerr << "Max value gap R0 " << maxi << std::endl;

  std::cerr << "Creatings blocks of R0..." << std::endl;
  std::vector< uint64_t > br_R0;
  uint64_t i; 
  for(i = br; i <= PB_R0.size(); i += br) {
    br_R0.push_back(PB_R0[i]);
  }
  
  if(PB_R0.size() % br != 0)
    br_R0.push_back(PB_R0[PB_R0.size() - 1]);

  block_r0 = _bv(br_R0.begin(), br_R0.end());

  n_block_r0 = br_R0.size();
  br_R0.clear();
  PB_R0.clear();

  std::cerr << "Calculating Gaps R1..." << std::endl;
  //GapsR1[0] = PB_R1[0] - 1;
  GapsR1[0] = PB_R1[0];
  uint32_t maxi_ = PB_R1[0];
  for(uint64_t i = 1; i < GapsR1.size(); i++) {
    //GapsR1[i] = PB_R1[i] - PB_R1[i - 1] - 1;
    GapsR1[i] = PB_R1[i] - PB_R1[i - 1];
    maxi_ = max(maxi_, GapsR1[i]);
  } 

  std::cerr << "Max value gap R0 " << maxi << std::endl;
  std::cerr << "Creatings blocks of R1..." << std::endl;
  std::vector< uint64_t > br_R1;
  for(i = br; i <= PB_R1.size(); i += br) {
    br_R1.push_back(PB_R1[i]);
  }
  
  if(PB_R1.size() % br != 0)
    br_R1.push_back(PB_R1[PB_R1.size() - 1]);
 
  block_r1 = _bv(br_R1.begin(), br_R1.end());

  n_block_r1 = br_R1.size();
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

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,_bv,_select,_rank>::get_runs(std::vector<uint64_t> &pb, std::vector<uint32_t> &res) {
  uint64_t sum = 0;
  res[0] = pb[0] - 1;
  for(uint64_t i = 1; i < pb.size(); i++) {
    res[i] = pb[i] - pb[i - 1] - 1;
    sum += pb[i] - pb[i - 1];
    if(res[i] < 0) std::cout << res[i] << "\n";
  }
  return sum;
}

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
void RunEncoderSelect<w,bs,br,_bv,_select,_rank>::top_k_encoding(std::vector< uint32_t > &seq, bool type) {
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
    } else {
      tc_or_huff_seq.push_back(i);
      hf_seq.push_back(seq[i]);
    }
  }
 
  if(type) {
    n_tchuff_r1 = tc_or_huff_seq.size();
    tc_or_huffman_r1 = _bv(tc_or_huff_seq.begin(), tc_or_huff_seq.end());
  } else {
    n_tchuff_r0 = tc_or_huff_seq.size();
    tc_or_huffman_r0 = _bv(tc_or_huff_seq.begin(), tc_or_huff_seq.end());
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

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,_bv,_select,_rank>::bits_tunstall_seq() {
  /*std::cout << "--- Gaps R0 ---" << std::endl;
  std::cout << "Dict size: " << tc_r0.dict_size() * 8 << " " << (double) tc_r0.dict_size() * 8 / n << std::endl;
  std::cout << "Compressed seq size: " << tc_r0.compressed_seq_size() * 8 << " " << (double) tc_r0.compressed_seq_size() * 8 / n << std::endl;
  std::cout << "Block size: " << tc_r0.block_vec_size() * 8 << " " << (double) tc_r0.block_vec_size() * 8 / n << std::endl;
  std::cout << "--- Gaps R1 ---" << std::endl;
  std::cout << "Dict size: " << tc_r1.dict_size() * 8 << " " << (double) tc_r1.dict_size() * 8 / n << std::endl;
  std::cout << "Compressed seq size: " << tc_r1.compressed_seq_size() * 8 << " " << (double) tc_r1.compressed_seq_size() * 8 / n << std::endl;
  std::cout << "Block size: " << tc_r1.block_vec_size() * 8 << " " << (double) tc_r1.block_vec_size() * 8 / n << std::endl;
  std::cout << "TC_R0: " << tc_r0_top_k.size() * 8 << std::endl;
  std::cout << "HF_R0: " << huffman_r0.size() * 8 << std::endl;
  std::cout << "TC_R1: " << tc_r1_top_k.size() * 8 << std::endl;
  std::cout << "HF_R1: " << huffman_r1.size() * 8 << std::endl;
*/
  return 8 * sdsl::size_in_bytes(tc_or_huffman_r1) + 
         8 * sdsl::size_in_bytes(rank_tchuff_r1) +
         8 * sdsl::size_in_bytes(select_tchuff_r1) +
         8 * sdsl::size_in_bytes(tc_or_huffman_r0) +
         8 * sdsl::size_in_bytes(rank_tchuff_r0) +
         8 * sdsl::size_in_bytes(select_tchuff_r0) +
         8 * tc_r0_top_k.size() + 8 * tc_r1_top_k.size() +
         8 * huffman_r0.size() + 8 * huffman_r1.size() +
         8 * size_block_r0() + 8 * size_block_r1();
}

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,_bv,_select,_rank>::size_block_r1() {
  return sdsl::size_in_bytes(block_r1) + sdsl::size_in_bytes(select_block_r1) + sdsl::size_in_bytes(rank_block_r1);
}

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,_bv,_select,_rank>::size_block_r0() {
  return sdsl::size_in_bytes(block_r0) + sdsl::size_in_bytes(select_block_r0) + sdsl::size_in_bytes(rank_block_r0);
}

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,_bv,_select,_rank>::select(uint64_t k) {
  assert(k <= n);

#ifdef SPLIT_TIME
  block_start = chrono::high_resolution_clock::now();
#endif
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

    // curr gap pos
    gap_pos = block * br + 1;

    // actual gaps in each structure
    /*
    if(gap_pos <= rank_tchuff_r0.size()) gap_huff_r0 = rank_tchuff_r0(gap_pos);
    else gap_huff_r0 = n_tchuff_r0;
    gap_tc_r0 = gap_pos - gap_huff_r0;
    one_r0 = gap_huff_r0; // are the same value in this point
    */

    if(gap_pos <= rank_tchuff_r1.size()) gap_huff_r1 = rank_tchuff_r1(gap_pos);
    else gap_huff_r1 = n_tchuff_r1;
    gap_tc_r1 = gap_pos - gap_huff_r1;
    one_r1 = gap_huff_r1; // are the same value in this point

    // pos before block
    pos += zeros;
    pos += ones;
  }

  bool flag_acum_r0 = false;
  bool flag_acum_r1 = false;
  uint64_t prev_r1 = -1;
  uint64_t res_select_r1 = -1;
  uint64_t prev_r0 = -1;
  uint64_t res_select_r0 = -1;

  if(block > 0) {
    if(one_r1 >= 1) {
      prev_r1 = select_tchuff_r1(one_r1);
      if(one_r1 + 1 <= n_tchuff_r1) res_select_r1 = select_tchuff_r1(one_r1 + 1);
      if(one_r1 + 1 > n_tchuff_r1 || prev_r1 + 1 != res_select_r1) {
        flag_acum_r1 = true;
        prev_r1 = gap_pos - 1;
        if(prev_r1 + 1 == res_select_r1) flag_acum_r1 = false;
      }
    } else {
      if(one_r1 + 1 <= n_tchuff_r1) res_select_r1 = select_tchuff_r1(one_r1 + 1);
      flag_acum_r1 = true;
      prev_r1 = gap_pos - 1;
    }
    one_r1++;

    /*
    if(one_r0 >= 1) {
      prev_r0 = select_tchuff_r0(one_r0);
      if(one_r0 + 1 <= n_tchuff_r0) res_select_r0 = select_tchuff_r0(one_r0 + 1);
      if(one_r0 + 1 > n_tchuff_r0 || prev_r0 + 1 != res_select_r0) {
        flag_acum_r0 = true;
        prev_r0 = gap_pos - 1;
        if(prev_r0 + 1 == res_select_r0) flag_acum_r0 = false;
      }
    } else {
      if(one_r0 + 1 <= n_tchuff_r0) res_select_r0 = select_tchuff_r0(one_r0 + 1);
      flag_acum_r0 = true;
      prev_r0 = gap_pos - 1;
    }
    one_r0++;
    */
  }

  /*
  uint32_t n_decoded_r0, cur_int_r0;
  uint32_t code_r0 = 0;
  uint32_t bits_needed_r0 = sizeof(uint32_t) << 3;
  uint32_t currcode_r0;
  uint32_t currlen_r0 = sizeof(uint32_t) << 3;
  uint32_t* lj_r0;
  uint32_t* start_linear_search_r0 = huffman_r0.lj_base + max(LUT_BITS, huffman_r0.min_cw_len) - 1;

  uint64_t symb_r0_huff;
  uint64_t i_r0 = gap_huff_r0;
  if(gap_huff_r0 == 0) {
    symb_r0_huff = 0;
    symb_r0_huff = huffman_r0.block_info_prefix_sum[gap_huff_r0/huffman_r0.block_size];
    cur_int_r0 = huffman_r0.block_info_starting_position[gap_huff_r0/huffman_r0.block_size];
    i_r0 = gap_huff_r0 - (gap_huff_r0/huffman_r0.block_size)*huffman_r0.block_size;

    huffman_r0.buff = huffman_r0.compressed_seq[cur_int_r0];
    huffman_r0.buff_btg = BUFF_BITS;

  } else {
    i_r0--;
    symb_r0_huff = huffman_r0.block_info_prefix_sum[i_r0/huffman_r0.block_size];
    cur_int_r0 = huffman_r0.block_info_starting_position[i_r0/huffman_r0.block_size];
    i_r0 = i_r0 - (i_r0/huffman_r0.block_size)*huffman_r0.block_size;

    huffman_r0.buff = huffman_r0.compressed_seq[cur_int_r0];
    huffman_r0.buff_btg = BUFF_BITS;

    for (n_decoded_r0 = 0; n_decoded_r0 <= i_r0; ++n_decoded_r0) {
      code_r0 |= huffman_r0.INPUT_ULONG(cur_int_r0, bits_needed_r0);

      lj_r0 = huffman_r0.lut[code_r0 >> ((sizeof(uint32_t) << 3) - LUT_BITS)];
      if (lj_r0 == NULL)
        for (lj_r0 = start_linear_search_r0; code_r0 < *lj_r0; lj_r0++)
          ;
      currlen_r0 = lj_r0 - huffman_r0.lj_base + 1;

      // calculate symbol number
      currcode_r0 = code_r0 >> ((sizeof(uint32_t) << 3) - currlen_r0);
      currcode_r0 -= huffman_r0.min_code[currlen_r0 - 1];
      currcode_r0 += huffman_r0.offset[currlen_r0 - 1];

      symb_r0_huff += huffman_r0.syms[currcode_r0]-1;

      code_r0 <<= currlen_r0;
      bits_needed_r0 = currlen_r0;
    }
    i_r0++;
  }
  */
  uint32_t n_decoded_r1, cur_int_r1;
  uint32_t code_r1 = 0;
  uint32_t bits_needed_r1 = sizeof(uint32_t) << 3;
  uint32_t currcode_r1;
  uint32_t currlen_r1 = sizeof(uint32_t) << 3;
  uint32_t* lj_r1;
  uint32_t* start_linear_search_r1 = huffman_r1.lj_base + max(LUT_BITS, huffman_r1.min_cw_len) - 1;

  uint64_t symb_r1_huff;
  uint64_t i_r1 = gap_huff_r1;
  if(gap_huff_r1 == 0) {
    symb_r1_huff = 0;
    symb_r1_huff = huffman_r1.block_info_prefix_sum[gap_huff_r1/huffman_r1.block_size];
    cur_int_r1 = huffman_r1.block_info_starting_position[gap_huff_r1/huffman_r1.block_size];
    i_r1 = gap_huff_r1 - (gap_huff_r1/huffman_r1.block_size)*huffman_r1.block_size;

    huffman_r1.buff = huffman_r1.compressed_seq[cur_int_r1];
    huffman_r1.buff_btg = BUFF_BITS;

  } else {
    i_r1--;
    symb_r1_huff = huffman_r1.block_info_prefix_sum[i_r1/huffman_r1.block_size];
    cur_int_r1 = huffman_r1.block_info_starting_position[i_r1/huffman_r1.block_size];
    i_r1 = i_r1 - (i_r1/huffman_r1.block_size)*huffman_r1.block_size;

    huffman_r1.buff = huffman_r1.compressed_seq[cur_int_r1];
    huffman_r1.buff_btg = BUFF_BITS;

    for (n_decoded_r1 = 0; n_decoded_r1 <= i_r1; ++n_decoded_r1) {
      code_r1 |= huffman_r1.INPUT_ULONG(cur_int_r1, bits_needed_r1);

      lj_r1 = huffman_r1.lut[code_r1 >> ((sizeof(uint32_t) << 3) - LUT_BITS)];
      if (lj_r1 == NULL)
        for (lj_r1 = start_linear_search_r1; code_r1 < *lj_r1; lj_r1++)
          ;
      currlen_r1 = lj_r1 - huffman_r1.lj_base + 1;

      // calculate symbol number
      currcode_r1 = code_r1 >> ((sizeof(uint32_t) << 3) - currlen_r1);
      currcode_r1 -= huffman_r1.min_code[currlen_r1 - 1];
      currcode_r1 += huffman_r1.offset[currlen_r1 - 1];

      symb_r1_huff += huffman_r1.syms[currcode_r1]-1;

      code_r1 <<= currlen_r1;
      bits_needed_r1 = currlen_r1;
    }
    i_r1++;
  }

  // cuidado cuando el gap_huff_r0 == 0
  //uint64_t prev_tunst_r0 = (gap_tc_r0 == 0 ? 0 : tc_r0_top_k.decode(gap_tc_r0 - 1));
  uint64_t prev_tunst_r1 = (gap_tc_r1 == 0 ? 0 : tc_r1_top_k.decode(gap_tc_r1 - 1));

  uint64_t symb_r1_tunst = prev_tunst_r1;
  //uint64_t symb_r0_tunst = prev_tunst_r0;
  while(symb_r1_tunst + symb_r1_huff < k) {
    /*
    if(take_gr0) {
      // read from gap of run 0
      uint64_t act_zeros = 0;
      if(!flag_acum_r0 && one_r0 <= n_tchuff_r0)
        res_select_r0 = select_tchuff_r0(one_r0);
      else if(one_r0 > n_tchuff_r0) flag_acum_r0 = true;

      if(!flag_acum_r0 && res_select_r0 == prev_r0 + 1) {
        if(i_r0 >= huffman_r0.block_size) {
          code_r0 = 0;
          bits_needed_r0 = sizeof(uint32_t) << 3;
          currlen_r0 = sizeof(uint32_t) << 3;
          start_linear_search_r0 = huffman_r0.lj_base + max(LUT_BITS, huffman_r0.min_cw_len) - 1;

          symb_r0_huff = huffman_r0.block_info_prefix_sum[gap_huff_r0/huffman_r0.block_size];
          cur_int_r0 = huffman_r0.block_info_starting_position[gap_huff_r0/huffman_r0.block_size];
          i_r0 = gap_huff_r0 - (gap_huff_r0/huffman_r0.block_size)*huffman_r0.block_size;

          huffman_r0.buff = huffman_r0.compressed_seq[cur_int_r0];
          huffman_r0.buff_btg = BUFF_BITS;
        }
        code_r0 |= huffman_r0.INPUT_ULONG(cur_int_r0, bits_needed_r0);

        lj_r0 = huffman_r0.lut[code_r0 >> ((sizeof(uint32_t) << 3) - LUT_BITS)];
        if (lj_r0 == NULL)
          for (lj_r0 = start_linear_search_r0; code_r0 < *lj_r0; lj_r0++)
            ;
        currlen_r0 = lj_r0 - huffman_r0.lj_base + 1;

        // calculate symbol number
        currcode_r0 = code_r0 >> ((sizeof(uint32_t) << 3) - currlen_r0);
        currcode_r0 -= huffman_r0.min_code[currlen_r0 - 1];
        currcode_r0 += huffman_r0.offset[currlen_r0 - 1];

        symb_r0_huff += huffman_r0.syms[currcode_r0]-1;

        code_r0 <<= currlen_r0;
        bits_needed_r0 = currlen_r0;

        i_r0++;
        gap_huff_r0++;
        one_r0++;
        prev_r0 = res_select_r0;

      } else {

        // is in tunstall
        flag_acum_r0 = true;

        uint64_t decode = tc_r0_top_k.decode(gap_tc_r0);
        act_zeros = decode - prev_tunst_r0;
        prev_tunst_r0 = decode;

        symb_r0_tunst += act_zeros;
        gap_tc_r0++;
        prev_r0++;
        if(prev_r0 + 1 == res_select_r0)
          flag_acum_r0 = false;
      }
    } else {
    */
      // read from gap of run 1
      uint64_t act_ones = 0;

      if(!flag_acum_r1 && one_r1 <= n_tchuff_r1)
        res_select_r1 = select_tchuff_r1(one_r1);
      else if(one_r1 > n_tchuff_r1) flag_acum_r1 = true;

      if(!flag_acum_r1 && res_select_r1 == prev_r1 + 1) {
        if(i_r1 >= huffman_r1.block_size) {
          code_r1 = 0;
          bits_needed_r1 = sizeof(uint32_t) << 3;
          currlen_r1 = sizeof(uint32_t) << 3;
          start_linear_search_r1 = huffman_r1.lj_base + max(LUT_BITS, huffman_r1.min_cw_len) - 1;

          symb_r1_huff = huffman_r1.block_info_prefix_sum[gap_huff_r1/huffman_r1.block_size];
          cur_int_r1 = huffman_r1.block_info_starting_position[gap_huff_r1/huffman_r1.block_size];
          i_r1 = gap_huff_r1 - (gap_huff_r1/huffman_r1.block_size)*huffman_r1.block_size;

          huffman_r1.buff = huffman_r1.compressed_seq[cur_int_r1];
          huffman_r1.buff_btg = BUFF_BITS;
        }
        code_r1 |= huffman_r1.INPUT_ULONG(cur_int_r1, bits_needed_r1);

        lj_r1 = huffman_r1.lut[code_r1 >> ((sizeof(uint32_t) << 3) - LUT_BITS)];
        if (lj_r1 == NULL)
          for (lj_r1 = start_linear_search_r1; code_r1 < *lj_r1; lj_r1++)
            ;
        currlen_r1 = lj_r1 - huffman_r1.lj_base + 1;

        // calculate symbol number
        currcode_r1 = code_r1 >> ((sizeof(uint32_t) << 3) - currlen_r1);
        currcode_r1 -= huffman_r1.min_code[currlen_r1 - 1];
        currcode_r1 += huffman_r1.offset[currlen_r1 - 1];

        symb_r1_huff += huffman_r1.syms[currcode_r1]-1;

        code_r1 <<= currlen_r1;
        bits_needed_r1 = currlen_r1;

        i_r1++;

        gap_huff_r1++;
        one_r1++;
        prev_r1 = res_select_r1;

      } else {

        // is in tunstall
        flag_acum_r1 = true;

        uint64_t decode = tc_r1_top_k.decode(gap_tc_r1);
        act_ones = decode - prev_tunst_r1;
        prev_tunst_r1 = decode;
        symb_r1_tunst += act_ones;

        gap_tc_r1++;
        prev_r1++;
        if(prev_r1 + 1 == res_select_r1)
          flag_acum_r1 = false;
      }
      gap_pos++;
      /*
    }
    take_gr0 = !take_gr0;
    */
  }

  if(gap_pos <= rank_tchuff_r0.size()) gap_huff_r0 = rank_tchuff_r0(gap_pos);
  else gap_huff_r0 = n_tchuff_r0;
  gap_tc_r0 = gap_pos - gap_huff_r0;

  uint64_t symb_r0_huff = (gap_huff_r0 == 0 ? 0 : huffman_r0.decode(gap_huff_r0 - 1));
  uint64_t symb_r0_tunst = (gap_tc_r0 == 0 ? 0 : tc_r0_top_k.decode(gap_tc_r0 - 1));
  pos = symb_r0_huff + symb_r0_tunst + symb_r1_huff + symb_r1_tunst;

  if(symb_r1_huff + symb_r1_tunst > k) {
    pos -= symb_r1_huff + symb_r1_tunst - k;
  }

  return pos - 2;
}

template< uint16_t w, uint64_t bs, uint64_t br, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,_bv,_select,_rank>::rank(uint64_t i) {
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

  bool take_gr0 = true;
  if(block > 0 && pos < i) {
    // curr gap pos
    gap_pos = block * br + 1;

    // actual gaps in each structure
    if(gap_pos <= rank_tchuff_r0.size()) gap_huff_r0 = rank_tchuff_r0(gap_pos);
    else gap_huff_r0 = n_tchuff_r0;
    gap_tc_r0 = gap_pos - gap_huff_r0;
    one_r0 = gap_huff_r0; // are the same value in this point

    if(gap_pos <= rank_tchuff_r1.size()) gap_huff_r1 = rank_tchuff_r1(gap_pos);
    else gap_huff_r1 = n_tchuff_r1;
    gap_tc_r1 = gap_pos - gap_huff_r1;
    one_r1 = gap_huff_r1; // are the same value in this point
  } else if(block == 1 && pos >= i) {
    _ones = 0;
    _zeros = 0;
    pos = 0;
  }

  bool flag_acum_r0 = false;
  bool flag_acum_r1 = false;
  uint64_t prev_r1 = -1;
  uint64_t res_select_r1 = -1;
  uint64_t prev_r0 = -1;
  uint64_t res_select_r0 = -1;

  if(block > 0 && pos < i) {
    if(one_r1 >= 1) {
      prev_r1 = select_tchuff_r1(one_r1);
      if(one_r1 + 1 <= n_tchuff_r1) res_select_r1 = select_tchuff_r1(one_r1 + 1);
      if(one_r1 + 1 > n_tchuff_r1 || prev_r1 + 1 != res_select_r1) {
        flag_acum_r1 = true;
        prev_r1 = gap_pos - 1;
        if(prev_r1 + 1 == res_select_r1) flag_acum_r1 = false;
      }
    } else {
      if(one_r1 + 1 <= n_tchuff_r1) res_select_r1 = select_tchuff_r1(one_r1 + 1);
      flag_acum_r1 = true;
      prev_r1 = gap_pos - 1;
    }
    one_r1++;

    if(one_r0 >= 1) {
      prev_r0 = select_tchuff_r0(one_r0);
      if(one_r0 + 1 <= n_tchuff_r0) res_select_r0 = select_tchuff_r0(one_r0 + 1);
      if(one_r0 + 1 > n_tchuff_r0 || prev_r0 + 1 != res_select_r0) {
        flag_acum_r0 = true;
        prev_r0 = gap_pos - 1;
        if(prev_r0 + 1 == res_select_r0) flag_acum_r0 = false;
      }
    } else {
      if(one_r0 + 1 <= n_tchuff_r0) res_select_r0 = select_tchuff_r0(one_r0 + 1);
      flag_acum_r0 = true;
      prev_r0 = gap_pos - 1;
    }
    one_r0++;
  }

  uint64_t prev_huff_r0 = (gap_huff_r0 == 0 ? 0 : huffman_r0.decode(gap_huff_r0 - 1));
  uint64_t prev_huff_r1 = (gap_huff_r1 == 0 ? 0 : huffman_r1.decode(gap_huff_r1 - 1));
  uint64_t prev_tunst_r0 = (gap_tc_r0 == 0 ? 0 : tc_r0_top_k.decode(gap_tc_r0 - 1));
  uint64_t prev_tunst_r1 = (gap_tc_r1 == 0 ? 0 : tc_r1_top_k.decode(gap_tc_r1 - 1));

  while(pos <= i) {
    if(take_gr0) {
      // read from gap of run 0
      uint64_t act_zeros = 0;
      if(!flag_acum_r0 && one_r0 <= n_tchuff_r0) res_select_r0 = select_tchuff_r0(one_r0);
      else if(one_r0 > n_tchuff_r0) flag_acum_r0 = true;

      if(!flag_acum_r0 && res_select_r0 == prev_r0 + 1) {
        // is in huffman

        uint64_t decode = huffman_r0.decode(gap_huff_r0);
        act_zeros = decode - prev_huff_r0;
        prev_huff_r0 = decode;

        gap_huff_r0++;
        one_r0++;
        prev_r0 = res_select_r0;
      } else {
        // is in tunstall
        flag_acum_r0 = true;

        uint64_t decode = tc_r0_top_k.decode(gap_tc_r0);
        act_zeros = decode - prev_tunst_r0;
        prev_tunst_r0 = decode;

        gap_tc_r0++;
        prev_r0++;
        if(prev_r0 + 1 == res_select_r0)
          flag_acum_r0 = false;
      }
      pos += act_zeros;
    } else {
      // read from gap of run 1
      uint64_t act_ones = 0;
      if(!flag_acum_r1 && one_r1 <= n_tchuff_r1) res_select_r1 = select_tchuff_r1(one_r1);
      else if(one_r1 > n_tchuff_r1) flag_acum_r1 = true;


      if(!flag_acum_r1 && res_select_r1 == prev_r1 + 1) {
        // is in huffman

        uint64_t decode = huffman_r1.decode(gap_huff_r1);
        act_ones = decode - prev_huff_r1;
        prev_huff_r1 = decode;

        gap_huff_r1++;
        one_r1++;
        prev_r1 = res_select_r1;
      } else {
        // is in tunstall
        flag_acum_r1 = true;

        uint64_t decode = tc_r1_top_k.decode(gap_tc_r1);
        act_ones = decode - prev_tunst_r1;
        prev_tunst_r1 = decode;

        gap_tc_r1++;
        prev_r1++;
        if(prev_r1 + 1 == res_select_r1)
          flag_acum_r1 = false;
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

template class RunEncoderSelect<16, 256, 32, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
template class RunEncoderSelect<16, 256, 64, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
template class RunEncoderSelect<16, 256, 128, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
template class RunEncoderSelect<16, 256, 256, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
template class RunEncoderSelect<16, 256, 512, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
template class RunEncoderSelect<16, 256, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 512, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 1024, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;

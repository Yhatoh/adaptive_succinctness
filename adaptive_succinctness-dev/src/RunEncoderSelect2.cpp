#include "RunEncoderSelect2.hpp"
using namespace std;

template< uint16_t w, uint64_t bs, uint64_t br, class gap_class, class _select_gap, class _bv, class _select, class _rank>
RunEncoderSelect<w,bs,br,gap_class,_select_gap,_bv,_select,_rank>::RunEncoderSelect(sdsl::bit_vector &bv, uint64_t top_k) {
  // to do later
}

template< uint16_t w, uint64_t bs, uint64_t br, class gap_class, class _select_gap, class _bv, class _select, class _rank>
RunEncoderSelect<w,bs,br,gap_class,_select_gap,_bv,_select,_rank>::RunEncoderSelect(std::vector<uint64_t> &pb, uint64_t top_k) {
#ifdef DEBUG
  cout << "Starting compression..." << endl;
#endif
  top_most_freq = top_k;
  n = pb.size();
  u = pb[pb.size() - 1] + 1;

  std::vector< uint64_t > R0, R1;

#ifdef DEBUG
  cout << "Getting runs of 0's and 1's..." << endl;
#endif
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

#ifdef DEBUG
  cout << "Prefix sum (runs of 0's)..." << endl;
#endif
  std::vector< uint64_t > PB_R0(R0.size(), 0);

  PB_R0[0] = R0[0];
  for(uint64_t i = 1; i < R0.size(); i++) {
    PB_R0[i] = R0[i] + PB_R0[i - 1];
  }

  R0.clear();

#ifdef DEBUG
  cout << "Prefix sum (runs of 1's)..." << endl;
#endif
  std::vector< uint64_t > PB_R1(R1.size(), 0);

  PB_R1[0] = R1[0];
  for(uint64_t i = 1; i < R1.size(); i++) {
    PB_R1[i] = R1[i] + PB_R1[i - 1];
  }

  R1.clear();

  n_r0 = PB_R0.size();
  n_r1 = PB_R1.size();

#ifdef DEBUG
  cout << "Getting Gaps from runs of 0's..." << endl;
#endif
  std::vector< uint32_t > GapsR0(PB_R0.size(), 0);
  std::vector< uint32_t > GapsR1(PB_R1.size(), 0);

  GapsR0[0] = PB_R0[0];
  for(uint64_t i = 1; i < GapsR0.size(); i++) {
    GapsR0[i] = PB_R0[i] - PB_R0[i - 1];
  }

#ifdef DEBUG
  cout << "Dividing gaps of 0's in blocks..." << endl;
#endif
  std::vector< uint64_t > br_R0;
  uint64_t i; 
  for(i = br; i < PB_R0.size(); i += br) {
    br_R0.push_back(PB_R0[i]);
  }
  
  // this is just if the final block is incomplete
  if(i - br != PB_R0.size() - 1) br_R0.push_back(PB_R0[PB_R0.size() - 1]);

#ifdef DEBUG
  cout << "Compressing blocks with a bit vector representation..." << endl;
#endif
  block_r0 = _bv(br_R0.begin(), br_R0.end());
  n_block_r0 = br_R0.size();

#ifdef DEBUG
  cout << "Nums of blocks: " << n_block_r0 << endl;
  cout << "Deleting auxiliar vectors..." << endl;
  cout << br_R0[15338 - 1] << endl;
  cout << br_R0[15338] << endl;
#endif
  br_R0.clear();
  PB_R0.clear();

#ifdef DEBUG
  cout << "Getting Gaps from runs of 1's..." << endl;
#endif
  GapsR1[0] = PB_R1[0];
  for(uint64_t i = 1; i < GapsR1.size(); i++) {
    GapsR1[i] = PB_R1[i] - PB_R1[i - 1];
  } 

  std::vector< uint64_t > br_R1;
  for(i = br; i < PB_R1.size(); i += br) {
    br_R1.push_back(PB_R1[i]);
  }
  
  // this is just if the final block is incomplete
  if(i - br != PB_R1.size() - 1) br_R1.push_back(PB_R1[PB_R1.size() - 1]);
 
  block_r1 = _bv(br_R1.begin(), br_R1.end());

  n_block_r1 = br_R1.size();

#ifdef DEBUG
  cout << "Nums of blocks: " << n_block_r1 << endl;
  cout << "Deleting auxiliar vectors..." << endl;
  cout << br_R1[15338 - 1] << endl;
  cout << br_R1[15338] << endl; // erase later
#endif

  br_R1.clear();
  PB_R1.clear();


#ifdef DEBUG
  cout << "Initiliaze support of rank and select operations on blocks..." << endl;
#endif
  sdsl::util::init_support(select_block_r0, &block_r0);
  sdsl::util::init_support(rank_block_r0, &block_r0);
  sdsl::util::init_support(select_block_r1, &block_r1);
  sdsl::util::init_support(rank_block_r1, &block_r1);

#ifdef DEBUG
  cout << "Applying top_k_encoding on gaps of 0's..." << endl;
#endif
  top_k_encoding(GapsR0, false);
#ifdef DEBUG
  cout << "Applying top_k_encoding on gaps of 1's..." << endl;
#endif
  top_k_encoding(GapsR1, true);

#ifdef DEBUG
  cout << "Initialize support of rank oper in tc || huffman bv..." << endl;
#endif
  sdsl::util::init_support(rank_tchuff_r0, &tc_or_huffman_r0);
  sdsl::util::init_support(rank_tchuff_r1, &tc_or_huffman_r1);
#ifdef DEBUG
  cout << "Initialize support of select oper in tc || huffman bv..." << endl;
#endif
  sdsl::util::init_support(select_tchuff_r0, &tc_or_huffman_r0);
  sdsl::util::init_support(select_tchuff_r1, &tc_or_huffman_r1);
}


template< uint16_t w, uint64_t bs, uint64_t br, class gap_class, class _select_gap, class _bv, class _select, class _rank>
void RunEncoderSelect<w,bs,br,gap_class,_select_gap,_bv,_select,_rank>::top_k_encoding(std::vector< uint32_t > &seq, bool type) {
#ifdef DEBUG
  cout << "Executing top_k_encoding..." << endl;
  cout << "Getting frequencies of gaps..." << endl;
#endif
  std::map< uint32_t, uint64_t > freq_map;
  for(uint64_t i = 0; i < seq.size(); i++) {
    if(freq_map.count(seq[i]) == 0) 
      freq_map[seq[i]] = 1;
    else
      freq_map[seq[i]] += 1;
  }
#ifdef DEBUG
  cout << "Sorting using a priority queue..." << endl;
#endif
  std::priority_queue< std::pair< uint32_t, uint64_t > > pq;
  for(std::map< uint32_t, uint64_t >::iterator it = freq_map.begin(); it != freq_map.end(); it++) {
    pq.push(std::pair< uint64_t, uint32_t >(it->second, it->first));
  }

#ifdef DEBUG
  cout << "Getting top k symbols..." << endl;
#endif
  std::set< uint32_t > tc_alphabet;
  for(uint64_t i = 0; i < top_most_freq; i++) {
    tc_alphabet.insert(pq.top().second);
    pq.pop();
  }

#ifdef DEBUG
  cout << "Building tunstall or huffman sequences..." << endl;
#endif 
  std::vector< uint64_t > tc_or_huff_seq;

  std::vector< uint32_t > tc_seq;
  std::vector< uint64_t > hf_seq;

  for(uint64_t i = 0; i < seq.size(); i++) {
    if(tc_alphabet.count(seq[i]) == 1) {
      tc_seq.push_back(seq[i]);
    } else {
      tc_or_huff_seq.push_back(i);
      hf_seq.push_back(seq[i]);
    }
  }
 
#ifdef DEBUG
  cout << "Building tunstall || huffman bv..." << endl;
#endif 
  if(type) {
    n_tchuff_r1 = tc_or_huff_seq.size();
    tc_or_huffman_r1 = _bv(tc_or_huff_seq.begin(), tc_or_huff_seq.end());
  } else {
    n_tchuff_r0 = tc_or_huff_seq.size();
    tc_or_huffman_r0 = _bv(tc_or_huff_seq.begin(), tc_or_huff_seq.end());
  }

  symbol_tc_p = (double)tc_seq.size() / seq.size() ;

#ifdef DEBUG
  cout << (type ? "RUNS OF 1's INFO" : "RUNS OF 0's INFO" ) << endl;
  cout << "  Symbols in tunstall: " << tc_seq.size() << endl;
  cout << "  Symbols in huffman: " << hf_seq.size() << endl;
  cout << "  % Symbols in tunstall: " << symbol_tc_p << endl;
  uint64_t sum = 0;
  for(uint64_t i = 0; i < (type ? 5343 : 25149); i++) {
    sum += hf_seq[i];
  }
  cout << sum << endl;
  sum = 0;
  for(uint64_t i = (type ? 5343 : 25149); i < (type ? 5346 : 25149); i++) {
    sum += hf_seq[i];
  }
  cout << sum << endl;

  sum = 0;
  for(uint64_t i = 0; i < (type ? 976290 : 956484); i++) {
    sum += tc_seq[i];
  }
  cout << sum << endl;
  sum = 0;
  for(uint64_t i = (type ? 976290 : 956484); i < (type ? 976494 : 956692); i++) {
    sum += tc_seq[i];
  }
  cout << sum << endl;
#endif

#ifdef DEBUG
  cout << "Encoding using tunstall and huffman..." << endl;
#endif 

  map< uint64_t, uint64_t > diff_symbol;
  for(auto &symbol : hf_seq) {
    diff_symbol[symbol]++;
  }

  map< uint64_t, uint64_t > encode_symbols;
  int encode = 1;
  for(auto &pos : diff_symbol) {
    encode_symbols[pos.first] = encode++;
  }

  vector< uint64_t > new_seq;
  for(auto &symbol : hf_seq) {
    new_seq.push_back(encode_symbols[symbol] + (new_seq.size() > 0 ? new_seq.back() : 0));
  }

  sdsl::bit_vector aux(new_seq.back() + 1, 0);
  for(auto &pos : new_seq) {
    aux[pos] = 1;
  }

#ifdef DEBUG
  cout << (type ? "RUNS OF 1's INFO" : "RUNS OF 0's INFO" ) << endl;
  cout << "Distinct Symbols in huff: " << encode << endl;
#endif

  if(type) { // r1
    tc_r1_top_k = tunstall_coder<w>(tc_seq, bs, 1 << w);

    dict_s9_r1.assign(encode, 0);
    for(auto &code : encode_symbols) {
      dict_s9_r1[code.second] = code.first;
    }

    s9_r1 = gap_class(aux);
    select_s9_r1 = _select_gap(&s9_r1);
  } else { // r0
    tc_r0_top_k = tunstall_coder<w>(tc_seq, bs, 1 << w);

    dict_s9_r0.assign(encode, 0);
    for(auto &code : encode_symbols) {
      dict_s9_r0[code.second] = code.first;
    }

    s9_r0 = gap_class(aux);
    select_s9_r0 = _select_gap(&s9_r0); 
  }
}

template< uint16_t w, uint64_t bs, uint64_t br, class gap_class, class _select_gap, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,gap_class,_select_gap,_bv,_select,_rank>::bits_tunstall_seq() {
  return 8 * sdsl::size_in_bytes(tc_or_huffman_r1) + 
         8 * sdsl::size_in_bytes(rank_tchuff_r1) +
         8 * sdsl::size_in_bytes(select_tchuff_r1) +
         8 * sdsl::size_in_bytes(tc_or_huffman_r0) +
         8 * sdsl::size_in_bytes(rank_tchuff_r0) +
         8 * sdsl::size_in_bytes(select_tchuff_r0) +
         8 * tc_r0_top_k.size() + 8 * tc_r1_top_k.size() +
         8 * sdsl::size_in_bytes(s9_r0) + 8 * sdsl::size_in_bytes(s9_r1) +
         8 * sizeof(uint64_t) * (dict_s9_r0.size() + dict_s9_r1.size()) + 
         8 * size_block_r0() + 8 * size_block_r1();
}

template< uint16_t w, uint64_t bs, uint64_t br, class gap_class, class _select_gap, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,gap_class,_select_gap,_bv,_select,_rank>::size_block_r1() {
  return sdsl::size_in_bytes(block_r1) +
         sdsl::size_in_bytes(select_block_r1) +
         sdsl::size_in_bytes(rank_block_r1);
}

template< uint16_t w, uint64_t bs, uint64_t br, class gap_class, class _select_gap, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,gap_class,_select_gap,_bv,_select,_rank>::size_block_r0() {
  return sdsl::size_in_bytes(block_r0) +
         sdsl::size_in_bytes(select_block_r0) +
         sdsl::size_in_bytes(rank_block_r0);
}

template< uint16_t w, uint64_t bs, uint64_t br, class gap_class, class _select_gap, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,gap_class,_select_gap,_bv,_select,_rank>::select(uint64_t k) {
  assert(k <= n);

  uint64_t block = rank_block_r1(k);
#ifdef DEBUG
  cout << "Block obtained: " << block << endl;
#endif

  uint64_t pos = 0;
  bool take_gr0 = true;

  // ones and zeros until block
  uint64_t ones = 0;
  uint64_t zeros = 0;

  // actual gap
  uint64_t gap_pos = 0;

  // ith one R0/R1
  uint64_t one_r0 = 1;
  uint64_t one_r1 = 1;

  // current gap in tunstall/huffman R0
  uint64_t gap_tc_r0 = 0;
  uint64_t gap_huff_r0 = 0;
  // current gap in tunstall/huffman R1
  uint64_t gap_tc_r1 = 0;
  uint64_t gap_huff_r1 = 0;
  
  if(block > 0) {
    // curr gap pos
    gap_pos = block * br + 1;

    // actual gaps in each structure
    if(gap_pos <= rank_tchuff_r1.size()) gap_huff_r1 = rank_tchuff_r1(gap_pos);
    else gap_huff_r1 = n_tchuff_r1;
    gap_tc_r1 = gap_pos - gap_huff_r1;
    one_r1 = gap_huff_r1; // are the same value in this point
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
  }

  uint64_t i_r1_tunst = gap_tc_r1;
  uint64_t b_, p, block_sum, symb_r1_tunst;
  uint64_t read_D;
  if(gap_tc_r1 == 0) {
    symb_r1_tunst = 0;
    block_sum = symb_r1_tunst = tc_r1_top_k.block[b_].prefix_sum;
    p = tc_r1_top_k.block[b_].starting_position;
    read_D = 0;
  } else {
    i_r1_tunst--;
    b_ = i_r1_tunst / tc_r1_top_k.bSize;

    symb_r1_tunst = tc_r1_top_k.block[b_].prefix_sum;
    block_sum = symb_r1_tunst;
    p = tc_r1_top_k.block[b_].starting_position;

    uint64_t j, size, nDecode = i_r1_tunst % tc_r1_top_k.bSize + 1;
    i_r1_tunst = i_r1_tunst % tc_r1_top_k.bSize;

    read_D = 0;
    for (j = 0; j <= nDecode; ++p) {
      size = tc_r1_top_k.D[tc_r1_top_k.compressed_seq[p]].size();
      if (j + size < nDecode) {
        symb_r1_tunst += tc_r1_top_k.D[tc_r1_top_k.compressed_seq[p]][size - 1];
        block_sum += tc_r1_top_k.D[tc_r1_top_k.compressed_seq[p]][size - 1];
        j += size;
      } else {
        symb_r1_tunst += tc_r1_top_k.D[tc_r1_top_k.compressed_seq[p]][nDecode - j - 1];
        read_D = nDecode - j - 1 + 1;
        break;
      }
    }
    i_r1_tunst++;
  }

  uint64_t symb_r1_huff = select_block_r1(block) - symb_r1_tunst;
  uint64_t prev_select = (gap_huff_r1 == 0 ? 0 : select_s9_r1(gap_huff_r1));

  uint64_t backup_gap_pos = gap_pos;

#ifdef DEBUG
  cout << "Nums of 1's until block " << block << endl;
  cout << "  Gap Huffman " << gap_huff_r1 << "\n";
  cout << "  Huffman     " << symb_r1_huff << endl;
  cout << "  Gap Tunstall " << gap_tc_r1 << "\n";
  cout << "  Tunstall    " << symb_r1_tunst << endl;
  cout << "  Total       " << symb_r1_huff + symb_r1_tunst << endl;
#endif
  while(symb_r1_tunst + symb_r1_huff < k) {
    // read from gap of run 1

    if(!flag_acum_r1 && one_r1 <= n_tchuff_r1)
      res_select_r1 = select_tchuff_r1(one_r1);
    else if(one_r1 > n_tchuff_r1) flag_acum_r1 = true;

    if(!flag_acum_r1 && res_select_r1 == prev_r1 + 1) {

      //symb_r1_huff = select_s9_r1(gap_huff_r1);
      uint64_t prefix_run = select_s9_r1(gap_huff_r1 + 1);
      uint64_t symbol = prefix_run - prev_select;
      prev_select = prefix_run;
      symb_r1_huff += dict_s9_r1[symbol];

      // sum 1's in huffman
      gap_huff_r1++;
      one_r1++;
      prev_r1 = res_select_r1;

    } else {

      // is in tunstall
      flag_acum_r1 = true;

      if(i_r1_tunst >= tc_r1_top_k.bSize) {
        b_ = gap_tc_r1 / tc_r1_top_k.bSize;
        symb_r1_tunst = block_sum = tc_r1_top_k.block[b_].prefix_sum;
        p = tc_r1_top_k.block[b_].starting_position;
        read_D = 0;
        i_r1_tunst = gap_tc_r1 % tc_r1_top_k.bSize;
      }
      if(read_D == tc_r1_top_k.D[tc_r1_top_k.compressed_seq[p]].size()) {
        block_sum += tc_r1_top_k.D[tc_r1_top_k.compressed_seq[p]][tc_r1_top_k.D[tc_r1_top_k.compressed_seq[p]].size() - 1];
        p++;
        read_D = 0;
      }
      symb_r1_tunst = block_sum + tc_r1_top_k.D[tc_r1_top_k.compressed_seq[p]][read_D++];
      i_r1_tunst++;

      gap_tc_r1++;
      prev_r1++;
      if(prev_r1 + 1 == res_select_r1)
        flag_acum_r1 = false;
    }
    gap_pos++;

  }

  if(gap_pos <= rank_tchuff_r0.size()) gap_huff_r0 = rank_tchuff_r0(gap_pos);
  else gap_huff_r0 = n_tchuff_r0;
  gap_tc_r0 = gap_pos - gap_huff_r0;

  uint64_t symb_r0_tunst = (gap_tc_r0 == 0 ? 0 : tc_r0_top_k.decode(gap_tc_r0 - 1));
 
  uint64_t gap_huff_r0_2 = 0;
  uint64_t gap_tc_r0_2 = 0;

  if(backup_gap_pos <= rank_tchuff_r0.size()) gap_huff_r0_2 = rank_tchuff_r0(backup_gap_pos);
  else gap_huff_r0_2 = n_tchuff_r0;
  gap_tc_r0_2 = backup_gap_pos - gap_huff_r0_2; 

  uint64_t symb_tunst_block = (gap_tc_r0_2 == 0 ? 0 : tc_r0_top_k.decode(gap_tc_r0_2 - 1));
  uint64_t symb_r0_huff = select_block_r0(block) - symb_tunst_block;
  
#ifdef DEBUG
  cout << "Nums of 0's until block " << block << endl;
  cout << "  Huffman   " << symb_r0_huff << endl;
  cout << "  Tunstall " << symb_tunst_block << endl;
  cout << "  Total    " << symb_r0_huff + symb_tunst_block << endl;
#endif

#ifdef DEBUG
  cout << "Gap Huffman R1 obj : " << gap_huff_r1 << "\n";
  cout << "Gap Tunstall R1 obj: " << gap_tc_r1 << "\n";
  cout << "Gap Huffman R0 now : " << gap_huff_r0_2 << "\n";
  cout << "Gap Huffman R0 obj : " << gap_huff_r0 << "\n";
  cout << "Gap Tunstall R0 now: " << gap_tc_r0_2 << "\n";
  cout << "Gap Tunstall R0 obj: " << gap_tc_r0 << "\n";
#endif
  prev_select = (gap_huff_r0_2 == 0 ? 0 : select_s9_r0(gap_huff_r0_2));

  for(gap_huff_r0_2 = gap_huff_r0_2; gap_huff_r0_2 < gap_huff_r0; gap_huff_r0_2++) {
    uint64_t prefix_run = select_s9_r0(gap_huff_r0_2 + 1);
    uint64_t symbol = prefix_run - prev_select;
    prev_select = prefix_run; 
    symb_r0_huff += dict_s9_r0[symbol];
  }

#ifdef DEBUG
  cout << "Symbols huff R0: " << symb_r0_huff << endl;
  cout << "Symbols tunst R0: " << symb_r0_tunst << endl;
  cout << "Total R0: " << symb_r0_tunst + symb_r0_huff << endl;
  cout << "Symbols huff R1: " << symb_r1_huff << endl;
  cout << "Symbols tunst R1: " << symb_r1_tunst << endl;
  cout << "Total R1: " << symb_r1_tunst + symb_r1_huff << endl;
#endif


  pos = symb_r0_huff + symb_r0_tunst + symb_r1_huff + symb_r1_tunst;

  if(symb_r1_huff + symb_r1_tunst > k) {
    pos -= symb_r1_huff + symb_r1_tunst - k;
  }

  return pos - 2;
}

template< uint16_t w, uint64_t bs, uint64_t br, class gap_class, class _select_gap, class _bv, class _select, class _rank>
uint64_t RunEncoderSelect<w,bs,br,gap_class,_select_gap,_bv,_select,_rank>::rank(uint64_t i) {
  /*
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
  */
  return 0;
}

//template class RunEncoderSelect<16, 64, 64, sdsl::s9_vector<128, sdsl::int_vector<32>>, sdsl::select_support_s9<1, 128, sdsl::int_vector<32>>, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 128, 64, sdsl::s9_vector<128, sdsl::int_vector<32>>, sdsl::select_support_s9<1, 128, sdsl::int_vector<32>>, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 256, 64, sdsl::s9_vector<128, sdsl::int_vector<32>>, sdsl::select_support_s9<1, 128, sdsl::int_vector<32>>, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
template class RunEncoderSelect<16, 64, 64, sdsl::s9_vector<>, sdsl::select_support_s9<>, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
template class RunEncoderSelect<16, 128, 64, sdsl::s9_vector<>, sdsl::select_support_s9<>, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
template class RunEncoderSelect<16, 256, 64, sdsl::s9_vector<>, sdsl::select_support_s9<>, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 256, 128, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 256, 256, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 256, 512, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 256, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 128, 32, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 128, 64, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 512, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;
//template class RunEncoderSelect<16, 1024, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>>;

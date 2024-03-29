#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/vectors.hpp>
#include <vector>
//#include "RunEncoderAccess.hpp"
#include "RunEncoderSelect.hpp"

int main(int argc, char **argv) {
  /*
     std::vector<uint64_t> seq;
     std::ifstream rf(argv[1], std::ios::binary);
     std::cerr << "Reading file..." << std::endl;
     uint64_t reader = 0;
     while(!rf.eof()) {
     uint64_t one;
     rf.read((char*) &one, sizeof(uint64_t));
     seq.push_back(one);
     reader++;
  //if(reader >= 1000000) break;
  }
  std::cerr << "Done reading file..." << std::endl;
  rf.close();
  */
  /*
  sdsl::bit_vector bv;
  bv.resize(100000);
  sdsl::util::set_random_bits(bv, 42);

  std::vector< uint64_t > seq;
  uint64_t count = 0;
  for(uint64_t i = 0; i < 100000; i++) {
    if(bv[i] == 1) {
      seq.push_back(i);
      count++;
    }
  }
  */
  
  bv.resize(61);
  std::vector< uint64_t > seq = {2,3,4,7,8,9,13,14,19,20,21,25,26,27,30,31,36,37,38,41,42,45,49,50,51,55,56,59,60};
  for(uint64_t i = 0; i < 61; i++) {
    bv[i] = 0;
  }

  for(uint64_t i = 0; i < seq.size(); i++) {
    bv[seq[i]] = 1;
  }
  sdsl::select_support_mcl<1> select(&bv);
  sdsl::rank_support_v5<1> rank(&bv);
  std::vector< uint64_t > ks = {8, 16, 32, 64};
  for(uint64_t k : ks) {
    std::cout << "TOP_K = " << k << std::endl;
    std::cerr << "Creating top_k = 16" << std::endl;
    RunEncoderSelect<16, 256, 512, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>> ge_16(seq, k);
    std::cout << "SELECT OP" << std::endl;
    for(uint64_t i = 1; i <= seq.size(); i++) {
      if(ge_16.select(i) != select(i)) {
        std::cout << "Query " << i << std::endl;
        std::cout << "SELECT RUNENCODER" << std::endl;
        std::cout << ge_16.select(i) << std::endl;
        std::cout << "SELECT REAL" << std::endl;
        std::cout << select(i) << std::endl;
        break;
      }
    }
    std::cout << "RANK OP" << std::endl;
    for(uint64_t i = 0; i <= seq[seq.size() - 1]; i++) {
      if(ge_16.rank(i) != rank(i)) {
        std::cout << "Query " << i << std::endl;
        std::cout << "RANK RUNENCODER" << std::endl;
        std::cout << ge_16.rank(i) << std::endl;
        std::cout << "rank REAL" << std::endl;
        std::cout << rank(i) << std::endl;
        break;
      }
    }
    std::cerr << "Done" << std::endl;
    std::cout << "INT VECTOR OF W = 16" << std::endl << (double) ge_16.bits_tunstall_seq() / seq.size() << std::endl;
  }
  return 0;
}

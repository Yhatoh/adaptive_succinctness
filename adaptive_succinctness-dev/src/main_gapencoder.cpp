#include <cmath>
#include <ctime>
#include <iostream>
#include <ratio>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/vectors.hpp>
#include <vector>
#include <string>

#include "../util/utils.hpp"
#include "RunEncoder.hpp"

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

  //std::vector< uint64_t > seq = {2,3,4,7,8,9,13,14,19,20,21,25,26,27,30,31,36,37,38,41,42,45,49,50,51,55,56,59,60};
  sdsl::select_support_mcl<1> select(&bv);
  std::vector< uint64_t > ks = {4};//, 8, 16, 32, 64};
  for(uint64_t k : ks) {
    std::cout << "TOP_K = " << k << std::endl;
    std::cerr << "Creating top_k = 16" << std::endl;
    RunEncoder<16, 1024, 512> ge_16(seq, k);
    std::cout << seq.size() << "\n";
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
    std::cerr << "Done" << std::endl;
    std::cout << "INT VECTOR OF W = 16" << std::endl << (double) ge_16.bits_tunstall_seq() / seq.size() << std::endl;
  }
  return 0;
}

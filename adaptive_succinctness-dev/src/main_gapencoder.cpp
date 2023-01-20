#include <cmath>
#include <ctime>
#include <iostream>
#include <ratio>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/vectors.hpp>
#include <vector>
#include <string>

#include "utils.hpp"
#include "GapEncoder.hpp"
#include "RunEncoder.hpp"

int main(int argc, char **argv) {
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
  /*
  GapEncoder<18> ge_18(seq);
  std::cout << "INT VECTOR OF W = 18 " << (double) ge_18.bits_tunstall_seq() / seq.size() << std::endl;
  GapEncoder<20> ge_20(seq);
  std::cout << "INT VECTOR OF W = 20 " << (double) ge_20.bits_tunstall_seq() / seq.size() << std::endl;
  GapEncoder<22> ge_22(seq);
  std::cout << "INT VECTOR OF W = 22 " << (double) ge_22.bits_tunstall_seq() / seq.size() << std::endl;
  GapEncoder<24> ge_24(seq);
  std::cout << "INT VECTOR OF W = 24 " << (double) ge_24.bits_tunstall_seq() / seq.size() << std::endl;
  RunEncoder<16> ge_16(seq);
  std::cout << "INT VECTOR OF W = 16 " << (double) ge_16.bits_tunstall_seq() / seq.size() << std::endl;
  */
  
  std::vector< uint64_t > ks = {4, 8, 16, 32, 64};
  for(uint64_t k : ks) {
    std::cout << "TOP_K = " << k << std::endl;
    std::cerr << "Creating top_k = 16" << std::endl;
    RunEncoder<16> ge_16(seq, k, 512);
    std::cerr << "Done" << std::endl;
    std::cout << "INT VECTOR OF W = 16" << std::endl << (double) ge_16.bits_tunstall_seq() / seq.size() << std::endl;
    std::cerr << "Creating top_k = 18" << std::endl;
    RunEncoder<18> ge_18(seq, k, 512);
    std::cerr << "Done" << std::endl;
    std::cout << "INT VECTOR OF W = 18 " << std::endl << (double) ge_18.bits_tunstall_seq() / seq.size() << std::endl;
    std::cerr << "Creating top_k = 20" << std::endl;
    RunEncoder<20> ge_20(seq, k, 512);
    std::cerr << "Done" << std::endl;
    std::cout << "INT VECTOR OF W = 20 " << std::endl << (double) ge_20.bits_tunstall_seq() / seq.size() << std::endl;
    std::cerr << "Creating top_k = 22" << std::endl;
    RunEncoder<22> ge_22(seq, k, 512);
    std::cerr << "Done" << std::endl;
    std::cout << "INT VECTOR OF W = 22 " << std::endl << (double) ge_22.bits_tunstall_seq() / seq.size() << std::endl;
    std::cerr << "Creating top_k = 24" << std::endl;
    RunEncoder<24> ge_24(seq, k, 512);
    std::cerr << "Done" << std::endl;
    std::cout << "INT VECTOR OF W = 24 " << std::endl << (double) ge_24.bits_tunstall_seq() / seq.size() << std::endl;
  }
  return 0;
}

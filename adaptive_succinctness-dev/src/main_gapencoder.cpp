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

int main(int argc, char **argv) {
  vector<uint64_t> seq;
  utils::read_input_file(argv[1], seq);
  GapEncoder ge(seq);
  std::cout << ge.bits_tunstall_seq() << "\n";
  return 0;
}
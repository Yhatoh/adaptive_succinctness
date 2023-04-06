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
#include "tunstallCoder.hpp"

int main() {
  std::vector< uint32_t > vec;
  vec.push_back(2);
  vec.push_back(4);
  vec.push_back(3);
  vec.push_back(5);
  vec.push_back(1);
  vec.push_back(2);
  vec.push_back(3);
  vec.push_back(5);

  // 0 1 2 3 4 5 6 7
  // 2 4 3 5 1 2 3 5
  // 3 
  tunstall_coder< 16 > t(vec,  512, 1);

  std::cout << t.decode(6) - t.decode(5) << "\n";
}

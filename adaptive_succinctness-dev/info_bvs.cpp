#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>


int main() {
  std::vector<uint64_t> seq;
  std::ifstream rf("/data/bitvectors/ii/gov2/url/gov2_ii_nofreq_url_dif.txt.dat.100000", std::ios::binary);
  //std::ifstream rf("/mnt/c/Users/gacar/Downloads/gov2_ii_nofreq_url_dif.txt.dat.100000", std::ios::binary);
  std::cout << "Reading file..." << std::endl;
  uint64_t reader = 0;
  while(!rf.eof()) {
    uint64_t one;
    rf.read((char*) &one, sizeof(uint64_t));
    seq.push_back(one);
    reader++;
    //if(reader >= 10000000) break;
  }
  seq.pop_back();
  std::cout << "Done reading file..." << std::endl;
  rf.close();

  uint64_t n = seq.size();
  uint64_t u = seq[seq.size() - 1] + 1;

  std::vector< uint64_t > R0, R1;

  uint64_t last = seq[0] + 1;
  R0.push_back(seq[0] + 1);

  for(uint64_t i = 1; i < seq.size(); i++) {
    if(seq[i] > seq[i - 1] + 1) {
      R1.push_back(seq[i - 1] + 1 - last + 1);
      R0.push_back(seq[i] - seq[i - 1] - 1);
      last = seq[i] + 1;
    }
  }

  R1.push_back(seq[seq.size() - 1] + 1 - last + 1);

  std::vector< uint64_t > PB_R0(R0.size(), 0);

  PB_R0[0] = R0[0];
  for(uint64_t i = 1; i < R0.size(); i++) {
    PB_R0[i] = R0[i] + PB_R0[i - 1];
  }

  R0.clear();

  std::vector< uint64_t > PB_R1(R1.size(), 0);

  PB_R1[0] = R1[0];
  for(uint64_t i = 1; i < R1.size(); i++) {
    PB_R1[i] = R1[i] + PB_R1[i - 1];
  }

  R1.clear();

  std::vector< uint32_t > GapsR0(PB_R0.size(), 0);
  std::vector< uint32_t > GapsR1(PB_R1.size(), 0);

  GapsR0[0] = PB_R0[0];
  for(uint64_t i = 1; i < GapsR0.size(); i++) {
    GapsR0[i] = PB_R0[i] - PB_R0[i - 1];
  }

  GapsR1[0] = PB_R1[0];
  for(uint64_t i = 1; i < GapsR1.size(); i++) {
    GapsR1[i] = PB_R1[i] - PB_R1[i - 1];
  } 

  PB_R1.clear();
  PB_R0.clear();

  std::map< uint64_t, uint64_t > countR1;

  for(auto gap : GapsR1) countR1[gap]++;

  for(auto &p : countR1) std::cout << p.first << " " << p.second << "\n";

  std::map< uint64_t, uint64_t > countR0;

  for(auto gap : GapsR0) countR0[gap]++;

  for(auto &p : countR0) std::cout << p.first << " " << p.second << "\n";
  return 0;
}

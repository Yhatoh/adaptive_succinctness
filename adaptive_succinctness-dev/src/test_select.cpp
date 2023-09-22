#include <cmath>
#include <ctime>
#include <iostream>
#include <ratio>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/vectors.hpp>
#include <vector>
#include <string>

#include <chrono>

#include "../util/utils.hpp"
#include "RunEncoderSelect.hpp"
#include "randomer.hpp"

#define n_queries 1000
#define U 196433254820
#define N 5055078462

using namespace std;

//typedef RunEncoderSelect<16, 256, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>> res_256_sd;
//typedef RunEncoderSelect<16, 256, 512, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>> res_256_sd;
typedef RunEncoderSelect<16, 256, 256, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>> res_256_sd;
typedef RunEncoderSelect<16, 512, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>> res_512_sd;
typedef RunEncoderSelect<16, 1024, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>> res_1024_sd;


int main() {
  std::vector<uint64_t> seq;
  //std::ifstream rf("/data/bitvectors/ii/gov2/url/gov2_ii_nofreq_url_dif.txt.dat.100000", std::ios::binary);
  std::ifstream rf("/mnt/c/Users/gacar/Downloads/gov2_ii_nofreq_url_dif.txt.dat.100000", std::ios::binary);
  std::cerr << "Reading file..." << std::endl;
  uint64_t reader = 0;
  while(!rf.eof()) {
    uint64_t one;
    rf.read((char*) &one, sizeof(uint64_t));
    seq.push_back(one);
    reader++;
    //if(seq.size() > 2 * 122850) break;
    if(reader >= 10000000) break;
  }
  std::cerr << "Done reading file..." << std::endl;
  rf.close(); 

  cout << seq.size() << "\n";
  sdsl::bit_vector bv(seq[seq.size() - 1] + 1, 0);
  for(auto bit : seq) bv[bit] = 1;

  Randomer r_randomer{0, seq[seq.size() - 1] - 1, 2048};
  Randomer s_randomer{1, seq.size() - 1, 2048};

  vector< uint64_t > vrank;
  vector< uint64_t > vselect;

  for(uint64_t q = 0; q < n_queries; q++){
    vrank.push_back(r_randomer());
    vselect.push_back(s_randomer());
  }
  /*
  {
    sdsl::sd_vector<> sd(bv);
    sdsl::rank_support_sd<1> rank_sd(&sd);
    sdsl::select_support_sd<1> select_sd(&sd);
    chrono::high_resolution_clock::time_point start, stop;
    double total_time = 0;
    chrono::duration< double > time_span;
    uint64_t Q = 0;
    start = chrono::high_resolution_clock::now();
    for(uint64_t q = 0; q < n_queries; q++){
      Q += rank_sd(vrank[q]);
    }
    stop = chrono::high_resolution_clock::now();
    time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
    total_time = time_span.count();
    cout << "SD RANK: " << total_time * 1000000 / n_queries << " " << total_time << " " << Q << endl;
    chrono::high_resolution_clock::time_point start_s, stop_s;
    double total_time_s = 0;
    chrono::duration< double > time_span_s;
    uint64_t Q_s = 0;
    start_s = chrono::high_resolution_clock::now();
    for(uint64_t q = 0; q < n_queries; q++){
      Q_s += select_sd(vselect[q]);
    }
    stop_s = chrono::high_resolution_clock::now();
    time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
    total_time_s = time_span_s.count();
    cout << "SD SELECT: " << Q_s << " " << total_time_s * 1000000 / n_queries << endl;
  }
  */
  std::vector< uint64_t > ks = {4, 8, 16, 32, 64};
  {
    for(uint64_t k : ks) {
      cout << "TEST H0RUN 256 " << k << "\n"; 
      res_256_sd run(seq, k);
      cout << "run sd 256 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << endl;
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t q = 0;
      start = chrono::high_resolution_clock::now();
      for(uint64_t i = 0; i < n_queries; i++){
        q += run.rank(vrank[i]);
      }
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();
      cout << "run sd 256 blocksize 512 top-" << k << " rank: " << q << " " << total_time*1000000/n_queries << endl;
      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t i = 0; i < vselect.size(); i++){
        q_s += run.select(vselect[i]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run sd 256 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s*1000000/n_queries << endl;
      cout << "BLOCK TIME: " << block_total_time * 1000000 / n_queries << "\n";
      cout << "TUNSTALL TIME " << tunst_total_time * 1000000 / n_queries << "\n";
      cout << "HUFFMAN TIME " << huff_total_time * 1000000 / n_queries << "\n";
      cout << "SELECT TIME " << select_total_time * 1000000 / n_queries << "\n";
      block_total_time = 0;
      tunst_total_time = 0;
      huff_total_time = 0;
      select_total_time = 0;
    }
  }
      /*
  {
    for(uint64_t k : ks) {
      cout << "TEST sd 512 " << k << "\n"; 
      res_512_sd run(seq, k);
      cout << "run sd 512 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t Q = 0;
      start = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        Q += run.rank(vrank[q]);
      }
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();
      cout << "run sd 512 blocksize 512 top-" << k << " rank: " << Q << " " << total_time*1000000/vrank.size() << "\n";
      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t Q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        Q_s += run.select(vselect[q]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run sd 512 blocksize 512 top-" << k << " select: " << Q_s << " " << total_time_s*1000000/vselect.size() << "\n";
    }
  }
  {
    for(uint64_t k : ks) {
      cout << "TEST sd 1024 " << k << "\n"; 
      res_1024_sd run(seq, k);
      cout << "run sd 1024 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t Q = 0;
      start = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        Q += run.rank(vrank[q]);
      }
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();
      cout << "run sd 1024 blocksize 512 top-" << k << " rank: " << Q << " " << total_time << "\n";
      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t Q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        Q_s += run.select(vselect[q]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run sd 1024 blocksize 512 top-" << k << " select: " << Q_s << " " << total_time_s << "\n";
    }
  }
    */
  return 0;
}

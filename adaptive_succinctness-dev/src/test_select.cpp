#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/io.hpp>
#include <sdsl/vectors.hpp>
#include <vector>
#include <string>

#include <chrono>

#include "RunEncoderSelect2.hpp"
//#include "RunEncoderAccess.hpp"
#include "randomer.hpp"

#define n_queries 100
#define U 196433254820
#define N 5055078462

using namespace std;

//typedef RunEncoderSelect<16, 256, 32,
//                         sdsl::sd_vector<>, sdsl::select_support_sd<1>,
//                         sdsl::rank_support_sd<1>> res_256_sd_32;'
//typedef RunEncoderSelect<16, 256, 64, sdsl::s9_vector<128, sdsl::int_vector<32>>, sdsl::select_support_s9<1, 128, sdsl::int_vector<32>>, sdsl::sd_vector<>,
//                         sdsl::select_support_sd<1>,
//                         sdsl::rank_support_sd<1>> res_256_sd_64;
//typedef RunEncoderSelect<16, 128, 64, sdsl::s9_vector<128, sdsl::int_vector<32>>, sdsl::select_support_s9<1, 128, sdsl::int_vector<32>>, sdsl::sd_vector<>,
//                         sdsl::select_support_sd<1>,
//                         sdsl::rank_support_sd<1>> res_128_sd_64;
//typedef RunEncoderSelect<16, 64, 64, sdsl::s9_vector<128, sdsl::int_vector<32>>, sdsl::select_support_s9<1, 128, sdsl::int_vector<32>>, sdsl::sd_vector<>,
//                         sdsl::select_support_sd<1>,
//                         sdsl::rank_support_sd<1>> res_64_sd_64;
typedef RunEncoderSelect<16, 256, 64, sdsl::s9_vector<>, sdsl::select_support_s9<>, sdsl::sd_vector<>,
                         sdsl::select_support_sd<1>,
                         sdsl::rank_support_sd<1>> res_256_sd_64;
typedef RunEncoderSelect<16, 128, 64, sdsl::s9_vector<>, sdsl::select_support_s9<>, sdsl::sd_vector<>,
                         sdsl::select_support_sd<1>,
                         sdsl::rank_support_sd<1>> res_128_sd_64;
typedef RunEncoderSelect<16, 64, 64, sdsl::s9_vector<>, sdsl::select_support_s9<>, sdsl::sd_vector<>,
                         sdsl::select_support_sd<1>,
                         sdsl::rank_support_sd<1>> res_64_sd_64;


//typedef RunEncoderAccess<16, 256, 64,
//                         sdsl::bit_vector,
//                         sdsl::rank_support_v5<1>,
//                         sdsl::sd_vector<>,
//                         sdsl::select_support_sd<1>,
//                         sdsl::rank_support_sd<1>> res_256_access_64;
/*
typedef RunEncoderSelect<16, 256, 128, sdsl::sd_vector<>,
                         sdsl::select_support_sd<1>,
                         sdsl::rank_support_sd<1>> res_256_sd_128;
typedef RunEncoderSelect<16, 256, 256, sdsl::sd_vector<>,
                         sdsl::select_support_sd<1>,
                         sdsl::rank_support_sd<1>> res_256_sd_256;
typedef RunEncoderSelect<16, 256, 512, sdsl::sd_vector<>,
                         sdsl::select_support_sd<1>,
                         sdsl::rank_support_sd<1>> res_64_sd_512;
                         */
//typedef RunEncoderSelect<16, 64, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>,
//                         sdsl::rank_support_sd<1>> res_64_sd_2048;
//typedef RunEncoderSelect<16, 512, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>> res_512_sd;
//typedef RunEncoderSelect<16, 1024, 2048, sdsl::sd_vector<>, sdsl::select_support_sd<1>, sdsl::rank_support_sd<1>> res_1024_sd;

void print_info(string type, uint64_t blockcode, uint64_t blockprefix,
                uint64_t k, double size, double size_per_one,
                double size_per_bit, uint64_t q_rank, double time_rank,
                uint64_t q_select, double time_select) {
  cout << "RUN " << type << endl;
  cout << "PARAMS: block code = " << blockcode << ", block prefix =  " << blockprefix << ", k = " << k << endl;
  cout << "  SIZE INFO: total = " << size << ", per one = " << size_per_one << ", per bit = " << size_per_bit << endl;
  cout << "  TIME INFO:" << endl;
  cout << "    RANK: " << q_rank << " time in microsecond = " << time_rank << endl;
  cout << "    SELECT: " << q_select << " time in microsecond = " << time_select << endl;
}

void print_info(string type, double size, double size_per_one,
                double size_per_bit, uint64_t q_rank, double time_rank,
                uint64_t q_select, double time_select) {
  cout << type << "\n";
  cout << "  SIZE INFO: total = " << size << ", per one = " << size_per_one << ", per bit = " << size_per_bit << endl;
  cout << "  TIME INFO:" << endl;
  cout << "    RANK: " << q_rank << " time in microsecond = " << time_rank << endl;
  cout << "    SELECT: " << q_select << " time in microsecond = " << time_select << endl;
}
int main() {
  std::vector<uint64_t> seq;
  std::ifstream rf("/data/bitvectors/ii/gov2/url/gov2_ii_nofreq_url_dif.txt.dat.100000", std::ios::binary);
  //std::ifstream rf("/mnt/c/Users/gacar/Downloads/gov2_ii_nofreq_url_dif.txt.dat.100000", std::ios::binary);
  std::cerr << "Reading file..." << std::endl;
  uint64_t reader = 0;
  while(!rf.eof()) {
    uint64_t one;
    rf.read((char*) &one, sizeof(uint64_t));
    seq.push_back(one);
    reader++;
    //if(reader >= 10000000) break;
  }
  std::cerr << "Done reading file..." << std::endl;
  rf.close();

  std::cerr << "Amount of 1's: " << seq.size() << "\n";
  std::cerr << "Universe: " << seq.back() + 1 << "\n";
  sdsl::bit_vector bv(seq[seq.size() - 1] + 1, 0);
  //sdsl::bit_vector bv(U, 0);
  for(auto bit : seq) bv[bit] = 1;

  //Randomer r_randomer{0, U - 1, 2048};
  //Randomer s_randomer{1, N - 1, 2048};

  Randomer r_randomer{0, seq[seq.size() - 1] - 1, 2048};
  Randomer s_randomer{1, seq.size() - 1, 2048};

  vector< uint64_t > vrank;
  vector< uint64_t > vselect;

  for(uint64_t q = 0; q < n_queries; q++){
    vrank.push_back(r_randomer());
    vselect.push_back(s_randomer());
  }
//  {
//    sdsl::sd_vector<> sd(bv);
//    sdsl::rank_support_sd<1> rank_sd(&sd);
//    sdsl::select_support_sd<1> select_sd(&sd);
//    chrono::high_resolution_clock::time_point start, stop;
//    double total_time = 0;
//    chrono::duration< double > time_span;
//    uint64_t Q = 0;
//    start = chrono::high_resolution_clock::now();
//    for(uint64_t q = 0; q < n_queries; q++){
//      Q += rank_sd(vrank[q]);
//    }
//    stop = chrono::high_resolution_clock::now();
//    time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
//    total_time = time_span.count();
//    chrono::high_resolution_clock::time_point start_s, stop_s;
//    double total_time_s = 0;
//    chrono::duration< double > time_span_s;
//    uint64_t Q_s = 0;
//    start_s = chrono::high_resolution_clock::now();
//    for(uint64_t q = 0; q < n_queries; q++){
//      Q_s += select_sd(vselect[q]);
//    }
//    stop_s = chrono::high_resolution_clock::now();
//    time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
//    total_time_s = time_span_s.count();
//    uint64_t total_size = (sdsl::size_in_bytes(sd) + sdsl::size_in_bytes(rank_sd) + sdsl::size_in_bytes(select_sd)) * 8;
//    print_info("SD VECTOR", total_size, (double) total_size / seq.size(), (double) total_size / (seq.back() + 1),
//               Q, (double) total_time * 1000000 / n_queries, Q_s, (double) total_time_s * 1000000 / n_queries);
//  }
  std::vector< uint64_t > ks = {64};//{4, 8, 16, 32, 64}; 
  {
    for(uint64_t k : ks) {
      res_64_sd_64 run(seq, k);

      cout << "FINISH" << endl;

      sdsl::sd_vector<> sd(bv);
      sdsl::rank_support_sd<1> rank_sd(&sd);
      sdsl::select_support_sd<1> select_sd(&sd);
      sdsl::rank_support_sd<0> rank_sd_0(&sd);
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t q = 0;
      start = chrono::high_resolution_clock::now();
      /*
      for(uint64_t i = 0; i < vrank.size(); i++){
        q += run.rank(vrank[i]);
      }
      */
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();


      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t i = 0; i < vselect.size(); i++) {
        auto res_run = run.select(vselect[i]);
        auto res_sd = select_sd(vselect[i]);
        if(res_run != res_sd) {
          cout << "TEST " << i << endl;
          cout << "TEST " << vselect[i] << endl;
          cout << "CORRECT " << res_sd << endl;
          cout << "FIND " << res_run << endl;
          cout << "AMOUNT OF ZEROS " << rank_sd_0(res_sd) << "\n";
          break;
        }
        q_s += run.select(vselect[i]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();

      print_info("SELECT", 64, 64, k, run.bits_tunstall_seq(),
                 (double) run.bits_tunstall_seq() / seq.size(),
                 (double) run.bits_tunstall_seq() / (seq.back() + 1),
                 q, total_time * 1000000 / n_queries, q_s, total_time_s * 1000000 / n_queries);
    }
  }
//  {
//    for(uint64_t k : ks) {
//      res_128_sd_64 run(seq, k);
//
//      sdsl::sd_vector<> sd(bv);
//      sdsl::rank_support_sd<1> rank_sd(&sd);
//      sdsl::select_support_sd<1> select_sd(&sd);
//      chrono::high_resolution_clock::time_point start, stop;
//      double total_time = 0;
//      chrono::duration< double > time_span;
//      uint64_t q = 0;
//      start = chrono::high_resolution_clock::now();
//      /*
//      for(uint64_t i = 0; i < vrank.size(); i++){
//        q += run.rank(vrank[i]);
//      }
//      */
//      stop = chrono::high_resolution_clock::now();
//      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
//      total_time = time_span.count();
//
//
//      chrono::high_resolution_clock::time_point start_s, stop_s;
//      double total_time_s = 0;
//      chrono::duration< double > time_span_s;
//      uint64_t q_s = 0;
//      start_s = chrono::high_resolution_clock::now();
//      for(uint64_t i = 0; i < vselect.size(); i++) {
//        q_s += run.select(vselect[i]);
//      }
//      stop_s = chrono::high_resolution_clock::now();
//      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
//      total_time_s = time_span_s.count();
//
//      print_info("SELECT", 128, 64, k, run.bits_tunstall_seq(),
//                 (double) run.bits_tunstall_seq() / seq.size(),
//                 (double) run.bits_tunstall_seq() / (seq.back() + 1),
//                 q, total_time * 1000000 / n_queries, q_s, total_time_s * 1000000 / n_queries);
//    }
//  }
//  {
//    for(uint64_t k : ks) {
//      res_256_sd_64 run(seq, k);
//
//      sdsl::sd_vector<> sd(bv);
//      sdsl::rank_support_sd<1> rank_sd(&sd);
//      sdsl::select_support_sd<1> select_sd(&sd);
//      chrono::high_resolution_clock::time_point start, stop;
//      double total_time = 0;
//      chrono::duration< double > time_span;
//      uint64_t q = 0;
//      start = chrono::high_resolution_clock::now();
//      /*
//      for(uint64_t i = 0; i < vrank.size(); i++){
//        q += run.rank(vrank[i]);
//      }
//      */
//      stop = chrono::high_resolution_clock::now();
//      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
//      total_time = time_span.count();
//
//
//      chrono::high_resolution_clock::time_point start_s, stop_s;
//      double total_time_s = 0;
//      chrono::duration< double > time_span_s;
//      uint64_t q_s = 0;
//      start_s = chrono::high_resolution_clock::now();
//      for(uint64_t i = 0; i < vselect.size(); i++) {
//        q_s += run.select(vselect[i]);
//      }
//      stop_s = chrono::high_resolution_clock::now();
//      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
//      total_time_s = time_span_s.count();
//
//      print_info("SELECT", 256, 64, k, run.bits_tunstall_seq(),
//                 (double) run.bits_tunstall_seq() / seq.size(),
//                 (double) run.bits_tunstall_seq() / (seq.back() + 1),
//                 q, total_time * 1000000 / n_queries, q_s, total_time_s * 1000000 / n_queries);
//    }
//  }
  return 0;
}

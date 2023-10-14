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
#include "RunEncoderAccess.hpp"
#include "randomer.hpp"

#define n_queries 1000000
#define U 196433254820
#define N 5055078462

using namespace std;

typedef RunEncoderAccess<16, 256, 523776, sdsl::bit_vector, sdsl::select_support_mcl<1>, sdsl::rank_support_v5<1>> rea_256_bv;
typedef RunEncoderAccess<16, 512, 523776, sdsl::bit_vector, sdsl::select_support_mcl<1>, sdsl::rank_support_v5<1>> rea_512_bv;
typedef RunEncoderAccess<16, 1024, 523776, sdsl::bit_vector, sdsl::select_support_mcl<1>, sdsl::rank_support_v5<1>> rea_1024_bv;
typedef RunEncoderAccess<16, 256, 523776, sdsl::rrr_vector<15>, sdsl::select_support_rrr<1,15>, sdsl::rank_support_rrr<1,15>> rea_256_rrr;
typedef RunEncoderAccess<16, 512, 523776, sdsl::rrr_vector<15>, sdsl::select_support_rrr<1,15>, sdsl::rank_support_rrr<1,15>> rea_512_rrr;
typedef RunEncoderAccess<16, 1024, 523776, sdsl::rrr_vector<15>, sdsl::select_support_rrr<1,15>, sdsl::rank_support_rrr<1,15>> rea_1024_rrr;
typedef RunEncoderAccess<16, 256, 523776, sdsl::rrr_vector<31>, sdsl::select_support_rrr<1,31>, sdsl::rank_support_rrr<1,31>> rea_256_rrr_31;
typedef RunEncoderAccess<16, 512, 523776, sdsl::rrr_vector<31>, sdsl::select_support_rrr<1,31>, sdsl::rank_support_rrr<1,31>> rea_512_rrr_31;
typedef RunEncoderAccess<16, 1024, 523776, sdsl::rrr_vector<31>, sdsl::select_support_rrr<1,31>, sdsl::rank_support_rrr<1,31>> rea_1024_rrr_31;


int main() {
  std::vector<uint64_t> seq;
  std::ifstream rf("/data/bitvectors/ii/gov2/url/gov2_ii_nofreq_url_dif.txt.dat.100000", std::ios::binary);
  std::cerr << "Reading file..." << std::endl;
  uint64_t reader = 0;
  while(!rf.eof()) {
    uint64_t one;
    rf.read((char*) &one, sizeof(uint64_t));
    seq.push_back(one);
    reader++;
  }
  std::cerr << "Done reading file..." << std::endl;
  rf.close();
  cout << seq.size() << "\n";
  sdsl::bit_vector bv(U, 0);
  for(auto bit : seq) bv[bit] = 1;

  Randomer r_randomer{0, U - 1};
  Randomer s_randomer{1, N - 1};

  vector< uint64_t > vrank;
  vector< uint64_t > vselect;

  for(uint64_t q = 0; q < n_queries; q++){
    vrank.push_back(r_randomer());
    vselect.push_back(s_randomer());
  }
  std::vector< uint64_t > ks = {4, 8, 16, 32, 64};

  {
    for(uint64_t k : ks) {
      cout << "TEST BV 256 " << k << "\n"; 
      rea_256_bv run(seq, k);
      cout << "run bv 256 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      /*
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
      cout << "run bv 256 blocksize 512 top-" << k << " rank: " << q << " " << total_time << "\n";

      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t i = 0; i < n_queries; i++){
        q_s += run.select(vselect[i]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run bv 256 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s << "\n";
      */
    }
  }
  {
    for(uint64_t k : ks) {
      cout << "TEST BV 512 " << k << "\n"; 
      rea_512_bv run(seq, k);
      /*
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t q = 0;
      start = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q += run.rank(vrank[q]);
      }
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();
      cout << "run bv 512 blocksize 512 top-" << k << " rank: " << q << " " << total_time << "\n";

      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q_s += run.select(vselect[q]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run bv 512 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s << "\n";
      cout << "run bv 512 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      */
    }
  }
  {
    for(uint64_t k : ks) {
      cout << "TEST bv 1024 " << k << "\n"; 
      rea_1024_bv run(seq, k);
      cout << "run bv 1024 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      /*
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t q = 0;
      start = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q += run.rank(vrank[q]);
      }
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();
      cout << "run bv 1024 blocksize 512 top-" << k << " rank: " << q << " " << total_time << "\n";

      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q_s += run.select(vselect[q]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run bv 1024 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s << "\n";
      */
    }
  }
  {
    for(uint64_t k : ks) {
      cout << "TEST RRR 256 " << k << "\n"; 
      rea_256_rrr run(seq, k);
      cout << "run rrr 15 256 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      /*
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
      cout << "run rrr 256 blocksize 512 top-" << k << " rank: " << q << " " << total_time << "\n";

      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t i = 0; i < n_queries; i++){
        q_s += run.select(vselect[i]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run rrr 256 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s << "\n";
      */
    }
  }
  {
    for(uint64_t k : ks) {
      cout << "TEST RRR 512 " << k << "\n"; 
      rea_512_rrr run(seq, k);
      cout << "run rrr 15 512 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      /*
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t q = 0;
      start = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q += run.rank(vrank[q]);
      }
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();
      cout << "run rrr 512 blocksize 512 top-" << k << " rank: " << q << " " << total_time << "\n";

      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q_s += run.select(vselect[q]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run rrr 512 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s << "\n";
      */
    }
  }
  {
    for(uint64_t k : ks) {
      cout << "TEST RRR 1024 " << k << "\n"; 
      rea_1024_rrr run(seq, k);
      cout << "run rrr 15 1024 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      /*
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t q = 0;
      start = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q += run.rank(vrank[q]);
      }
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();
      cout << "run rrr 1024 blocksize 512 top-" << k << " rank: " << q << " " << total_time << "\n";

      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q_s += run.select(vselect[q]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run rrr 1024 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s << "\n";
      */
    }
  }
  {
    for(uint64_t k : ks) {
      cout << "TEST RRR 256 " << k << "\n"; 
      rea_256_rrr_31 run(seq, k);
      cout << "run rrr 31 256 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      /*
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
      cout << "run rrr 256 blocksize 512 top-" << k << " rank: " << q << " " << total_time << "\n";

      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t i = 0; i < n_queries; i++){
        q_s += run.select(vselect[i]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run rrr 256 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s << "\n";
      */
    }
  }
  {
    for(uint64_t k : ks) {
      cout << "TEST RRR 512 " << k << "\n"; 
      rea_512_rrr_31 run(seq, k);
      cout << "run rrr 31 512 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      /*
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t q = 0;
      start = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q += run.rank(vrank[q]);
      }
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();
      cout << "run rrr 512 blocksize 512 top-" << k << " rank: " << q << " " << total_time << "\n";

      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q_s += run.select(vselect[q]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run rrr 512 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s << "\n";
      */
    }
  }
  {
    for(uint64_t k : ks) {
      cout << "TEST RRR 1024 " << k << "\n"; 
      rea_1024_rrr_31 run(seq, k);
      cout << "run rrr 31 1024 blocksize 512 top-" << k << " space: " << (double) run.bits_tunstall_seq() / seq.size() << "\n";
      /*
      chrono::high_resolution_clock::time_point start, stop;
      double total_time = 0;
      chrono::duration< double > time_span;
      uint64_t q = 0;
      start = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q += run.rank(vrank[q]);
      }
      stop = chrono::high_resolution_clock::now();
      time_span = chrono::duration_cast< chrono::microseconds >(stop - start);
      total_time = time_span.count();
      cout << "run rrr 1024 blocksize 512 top-" << k << " rank: " << q << " " << total_time << "\n";

      chrono::high_resolution_clock::time_point start_s, stop_s;
      double total_time_s = 0;
      chrono::duration< double > time_span_s;
      uint64_t q_s = 0;
      start_s = chrono::high_resolution_clock::now();
      for(uint64_t q = 0; q < n_queries; q++){
        q_s += run.select(vselect[q]);
      }
      stop_s = chrono::high_resolution_clock::now();
      time_span_s = chrono::duration_cast< chrono::microseconds >(stop_s - start_s);
      total_time_s = time_span_s.count();
      cout << "run rrr 1024 blocksize 512 top-" << k << " select: " << q_s << " " << total_time_s << "\n";
      */
    }
  }
  return 0;
}

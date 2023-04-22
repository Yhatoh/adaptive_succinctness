#ifndef UTILS
#define UTILS

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <ctime>

uint64_t alphabet_size(std::vector<uint32_t> &V);
float H0(std::vector<uint32_t> &V, bool print);

namespace utils {
    void read_input_file(const char *input_file, std::vector<uint64_t> &vec_ref);
    void write_results(bool header, const std::string s,
                   const std::vector<std::string> &l,
                   const std::vector<uint64_t> &t_seq_bytes,
                   const std::vector<uint64_t> &h_seq_bytes,
                   const std::vector<uint64_t> &bit_vec_bytes,
                   const std::vector<float> &t_seq_p,
                   const std::vector<float> &h_seq_p,
                   const std::vector<float> &bit_vec_p,
                   const std::vector<uint64_t> &total_bytes,
                   const std::vector<float> &symbol_p,
                   const std::vector<std::vector<float>> &bits_per_post);

    void separate_runs_write_results(bool header, const std::string s, int number_of_rows,
                   const std::vector<std::string> &l,
                   const std::vector<uint64_t> &total_bytes,
                   const std::vector<std::vector<float>> &symbol_p,
                   const std::vector< std::vector< std::vector<float> > > &pct,
                   const std::vector< std::vector< std::vector<float> > > &bits_per_post,
                   const std::vector<std::string> &B_labels);
                   
    std::string current_time2str();
    void write_gapfreqs2file(const std::string outputfile, const std::map<uint32_t, uint64_t> &freqs);
}

#endif
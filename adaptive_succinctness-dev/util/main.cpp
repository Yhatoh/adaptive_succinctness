
#include <cmath>
#include <ctime>
#include <iostream>
#include <ratio>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/vectors.hpp>
#include <vector>
#include <string>

//#include "huffman_coder_noblockcompression.cpp"
//#include <stdio.h>

//#include <chrono>
//#include <ratio>

//#include "hybrid_coder.hpp"
//#include "randomer.hpp"
#include "utils.hpp"
#include "h0run.hpp"

using namespace std::chrono;

using namespace std;
//using namespace sdsl;


// creating new vector v from w -> v[i] = w[i] - w[i-1]
void create_diff_vector(const vector<uint64_t> &seq_ref2, vector<uint64_t> &seq_diff_ref, uint64_t &sumOrig_r) {
    cout << "Creating difference vector..." << endl;
    for (uint64_t i = 1; i < seq_ref2.size(); i++) {
        if (seq_ref2[i] > seq_ref2[i - 1]) {
            seq_diff_ref[i] = seq_ref2[i] - seq_ref2[i - 1];
            sumOrig_r += (seq_ref2[i] - seq_ref2[i - 1]);
        } else {
            seq_diff_ref[i] = seq_ref2[i];
            sumOrig_r += seq_ref2[i];
        }
    }
}

sdsl::bit_vector create_bit_vector(const vector<uint64_t> &vec_ref, const uint64_t &n) {
    cout << "Creating bit vector of " << n << " bits" << endl;
    sdsl::bit_vector bv = sdsl::bit_vector(n + 1);

    for (uint64_t i = 0; i < n; i++) bv[i] = 0;

    uint64_t sum = vec_ref[0];
    bv[sum] = 1;
    for (uint64_t i = 1; i < vec_ref.size(); i++) {
        sum += vec_ref[i];
        bv[sum] = 1;
    }

    return bv;
}

int main(int argc, char **argv) {
    vector<uint64_t> seq_diff;

    if (argc < 3) {
        cout << "Error occured, not enough program parameters \n";
        cout << ">> Program usage: " << argv[1] << " <input_file> <mode>\n";
        cout << "\t|| mode = 1 -> input file is original (sequence of 64-bit integers)\n";
        cout << "\t|| mode = 2 -> input file is preprosessed (gaps of orig file)\n";
        exit(1);
    }

    string curr_time = utils::current_time2str();

    cout << "Running on mode: " << argv[2] << "\n";
    uint64_t sumOrig, sum;

    string mode = (string) argv[2];

    if (mode == "1") { // reading input file and preprosessing it to gap form
        vector<uint64_t> seq;
        utils::read_input_file(argv[1], seq);
        seq_diff.resize(seq.size());
        seq_diff[0] = seq[0];
        sumOrig = seq[0];
        create_diff_vector(seq, seq_diff, sumOrig);
        seq.clear();
        seq.shrink_to_fit();
    } else if (mode == "2") { // input data already gap form
        utils::read_input_file(argv[1], seq_diff);
        sumOrig = seq_diff[0];
        for (uint64_t i = 1; i < seq_diff.size(); i++) sumOrig += seq_diff[i];
    }

    uint64_t n = seq_diff.size(), u;

    // creating bit_vector
    sdsl::bit_vector bv = create_bit_vector(seq_diff, sumOrig);
    seq_diff.clear();
    seq_diff.shrink_to_fit();

    cout << "Done..." << endl;

    h0run::separate_runs(bv, n, curr_time, false);
    //h0run::runs_together(bv, n, curr_time);

    // cout << "Huffman compression..." << endl;

    // high_resolution_clock::time_point start, stop;
    // double total_time = 0.0;
    // duration<double> time_span;

    // uint8_t b;
    // uint64_t block_size;

    // Randomer randomer_G{0, bv_G.size() - 1};
    // Randomer randomer_L{0, bv_L.size() - 1};

    // vector<uint64_t> random_position_G, random_position_L;
    // uint32_t random_vec_size = 5000000;
    // random_position_G.resize(random_vec_size);
    // random_position_L.resize(random_vec_size);

    // for (uint64_t i = 0; i < 5000000; i++) {
    //     random_position_G[i] = randomer_G();
    //     random_position_L[i] = randomer_L();
    // }

    /*    for (b = 3; b <= 16; ++b) {
            total_time = 0.0;

            block_size = std::pow(2,  b);

            cout << "--------------------------------------------" << endl;
            cout << "Block Size = " << block_size << endl;


            {
                huffman_coder setG, setL;

                setG.encode(bv_G, block_size);
                setL.encode(bv_L, block_size);

                cout << "Huffman(G)+Huffman(L) = " << ((float)setG.size_bits() + setL.size_bits())/n << endl;

                for (i = 0; i < random_position_G.size(); i++) {
                    start = high_resolution_clock::now();

                        j = setG.decode(random_position_G[i]);

                    stop = high_resolution_clock::now();
                    time_span = duration_cast<microseconds>(stop - start);
                    total_time += time_span.count();

                    if (bv_G[random_position_G[i]] != j) { // dummy check, to avoid the optimizer
                        continue;
                    }
                }

                for (i = 0; i < random_position_L.size(); i++) {
                    start = high_resolution_clock::now();

                        j = setL.decode(random_position_L[i]);

                    stop = high_resolution_clock::now();
                    time_span = duration_cast<microseconds>(stop - start);
                    total_time += time_span.count();

                    if (bv_L[random_position_L[i]] != j) { // dummy check to avoid the optimizer
                        continue;
                    }
                }

                cout << "Time per access = " << total_time/(random_position_G.size() + random_position_L.size()) << endl;
            }
        }
    */

    return 0;
}

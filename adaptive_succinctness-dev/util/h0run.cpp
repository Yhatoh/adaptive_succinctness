#include "h0run.hpp"

using namespace std;
// using namespace sdsl;

void create_R0_and_R1(sdsl::bit_vector &X, vector<uint32_t> &R0, vector<uint32_t> &R1) {
    cout << "Creating sets R0 and R1" << endl;
    uint64_t u = X.size();
    uint64_t g = 0;
    uint64_t nz = 0, no = 0;
    uint8_t current = X[0];

    current = X[0];
    for (uint64_t i = 0; i < u; ++i) {
        if (current == 0 && X[i] == 1) {
            current = 1;
            R0.push_back(nz);
            no = 1;
            g++;
        } else if (current == 0 && X[i] == 0) {
            // current = 0;
            nz++;
        } else if (current == 1 && X[i] == 1) {
            // current = 1;
            no++;
        } else if (current == 1 && X[i] == 0) {
            current = 0;
            R1.push_back(no);
            nz = 1;
        }
    }
    if (current == 1 && no != 1) R1.push_back(no);
    else if (current == 0 && nz != 1) R0.push_back(nz);
    cout << "Done..." << endl;
}

namespace h0run {
void runs_together(sdsl::bit_vector &X, uint64_t posting_list_size, const string curr_time) {
    cout << "Creating set R" << endl;
    uint64_t u = X.size();
    uint64_t g = 0;
    uint64_t counter = 0;
    uint8_t current = X[0];

    vector<uint32_t> bv_R;

    current = X[0];
    for (uint64_t i = 0; i < u; ++i) {
        if (current == 0 && X[i] == 1) {
            current = 0;
            bv_R.push_back(counter);
            counter = 1;
        } else if (current == 0 && X[i] == 0) {
            // current = 0;
            counter++;
        } else if (current == 1 && X[i] == 1) {
            // current = 1;
            counter++;
        } else if (current == 1 && X[i] == 0) {
            current = 0;
            bv_R.push_back(counter);
            counter = 1;
        }
    }

    // last run
    if (counter != 1) bv_R.push_back(counter);

    cout << "Done..." << endl;

    cout << "Set R \n      + nElems = " << bv_R.size() << " alphabet size = " << alphabet_size(bv_R) << " H0 = "\
        << H0(bv_R, false) / posting_list_size << endl;

    cout << "Tunstall compression..." << endl;

    vector<uint64_t> top_ks = {4, 8, 16, 32, 64};

    vector<uint64_t> tunstall_seq_bytes;
    vector<uint64_t> huffman_seq_bytes;
    vector<uint64_t> bit_vec_b_bytes;
    vector<float> tunstall_seq_p;
    vector<float> huffman_seq_p;
    vector<float> bit_vec_b_p;
    vector<float> tunstall_bpp;
    vector<float> huffman_bpp;
    vector<float> bit_vec_b_bpp;
    vector<float> total_bpp;
    vector<uint64_t> hybrid_bytes;
    vector<string> labels;
    vector<float> symbols_in_tuns;
    float curr_symbol_pct;

    for (uint64_t i = 0; i < top_ks.size(); i++) {
        cout << "\n##------##" << endl;
        cout << "Top-" << top_ks[i] << endl;
        cout << "T_R:" << endl;
        hybrid_coder T_R(bv_R, 512, top_ks[i], curr_symbol_pct, true);
        T_R.print_info();

        cout << "Total size Tunstal top-" << top_ks[i] << ": " << ((float)T_R.size_bits()) / posting_list_size << endl;

        labels.push_back("R-top-" + to_string(top_ks[i]));
        T_R.fetch_results(tunstall_seq_bytes, tunstall_seq_p, huffman_seq_bytes, huffman_seq_p, bit_vec_b_bytes,
                          bit_vec_b_p, hybrid_bytes);
        tunstall_bpp.push_back(T_R.bits_tunstall_seq());
        huffman_bpp.push_back(T_R.bits_huffman_seq());
        bit_vec_b_bpp.push_back(T_R.bits_bit_vec_b());
        total_bpp.push_back(((float)T_R.size_bits()) / posting_list_size);
        symbols_in_tuns.push_back(curr_symbol_pct);
    }

    vector< vector<float>> bpp = {tunstall_bpp, huffman_bpp,
        bit_vec_b_bpp, total_bpp};

    utils::write_results(false, "./data/result_output" + curr_time + ".csv",
                         labels,
                         tunstall_seq_bytes,
                         huffman_seq_bytes,
                         bit_vec_b_bytes,
                         tunstall_seq_p,
                         huffman_seq_p,
                         bit_vec_b_p,
                         hybrid_bytes, symbols_in_tuns, bpp);

    bv_R.clear();
}

void separate_runs(sdsl::bit_vector &X, uint64_t posting_list_size, const string curr_time, bool print_res) {
    vector<uint32_t> bv_R0;
    vector<uint32_t> bv_R1;
    create_R0_and_R1(X, bv_R0, bv_R1);
    

    cout << "Set R0 \n\t+ nElems = " << bv_R0.size() << " alphabet size = " << alphabet_size(bv_R0) << " H0 = "\
        << H0(bv_R0, false) / posting_list_size << endl;
    cout << "Set R1 \n\t+ nElems = " << bv_R1.size() << " alphabet size = " << alphabet_size(bv_R1) << " H0 = "\
        << H0(bv_R1, false) / posting_list_size << endl;

    cout << "Tunstall compression..." << endl;

    vector<uint64_t> top_ks = {4, 8, 16, 32, 64};
    //vector<uint64_t> top_ks = {64};
    //vector<uint64_t> top_ks;
    //for (uint64_t i = 4; i < 64; i++) top_ks.push_back(i);


    vector< vector<float> > tunstall_seq_pct = {{}, {}};
    vector< vector<float> > huffman_seq_pct = {{}, {}};
    vector< vector<float> > bit_vec_b_pct = {{}, {}};
    vector< vector<float> > tunstall_bpp = {{}, {}};
    vector< vector<float> > huffman_bpp = {{}, {}};
    vector< vector<float> > bit_vec_b_bpp = {{}, {}};
    vector <vector<float> > total_bpp = {{}, {}, {}};
    vector< vector<float>> symbol_pct = {{}, {}};
    vector<uint64_t> hybrid_bytes;
    vector<string> labels;
    float curr_symbols_in_tuns;
    vector<float> all_symbols_pct;
    vector<string> bit_vec_B;

    total_bpp.push_back({});
    int rows = 0;
    // hybrid coder with elias-fano bit vector B
    for (uint64_t i = 0; i < top_ks.size(); i++) {
        hybrid_coder T_R0(bv_R0, 512, top_ks[i], curr_symbols_in_tuns, true);
        float r0_tuns_sybmols_pct = curr_symbols_in_tuns;
        for (uint64_t j = 0; j < top_ks.size(); j++) {
            rows++;
            cout << "T_R0: Top-" << top_ks[i] << " T_R1: Top-" << top_ks[j] << " ELIAS_FANO" << endl;
            if (print_res) {
                cout << "T_R0 info:" << endl;
                T_R0.print_info();
                cout << "T_R1 info:" << endl;
            }
            hybrid_coder T_R1(bv_R1, 512, top_ks[i], curr_symbols_in_tuns, true);
            if (print_res) {
                T_R1.print_info();
                cout << endl;
                cout << "Tunstall of top-" << top_ks[i] << " most-frequent runs of 0: " << (float)T_R0.size_bits() / posting_list_size << endl;
                cout << "Tunstall of top-" << top_ks[i] << " most-frequent runs of 1: " << (float)T_R1.size_bits() / posting_list_size << endl;
                cout << "T_R0 size from total (%): " << ((float)T_R0.size_bits()) / ((float)T_R0.size_bits() + T_R1.size_bits()) << endl;
                cout << "T_R1 size from total (%): " << ((float)T_R1.size_bits()) / ((float)T_R0.size_bits() + T_R1.size_bits()) << endl;
                cout << "Total size Tunstal: " << ((float)T_R0.size_bits() + T_R1.size_bits()) / posting_list_size << endl;
            }
            
            tunstall_bpp[0].push_back((float)T_R0.bits_tunstall_seq() / posting_list_size);
            huffman_bpp[0].push_back((float)T_R0.bits_huffman_seq() / posting_list_size);
            bit_vec_b_bpp[0].push_back((float)T_R0.bits_bit_vec_b() / posting_list_size);

            tunstall_bpp[1].push_back((float)T_R1.bits_tunstall_seq() / posting_list_size);
            huffman_bpp[1].push_back((float)T_R1.bits_huffman_seq() / posting_list_size);
            bit_vec_b_bpp[1].push_back((float)T_R1.bits_bit_vec_b() / posting_list_size);
            
            total_bpp[0].push_back( ( (float)T_R0.size_bits() ) / posting_list_size );
            total_bpp[1].push_back( ( (float)T_R1.size_bits() ) / posting_list_size );
            total_bpp[2].push_back( ( (float)T_R0.size_bits() + T_R1.size_bits() ) / posting_list_size );

            tunstall_seq_pct[0].push_back((float)T_R0.bits_tunstall_seq() / T_R0.size_bits());
            huffman_seq_pct[0].push_back((float)T_R0.bits_huffman_seq() / T_R0.size_bits());
            bit_vec_b_pct[0].push_back((float)T_R0.bits_bit_vec_b() / T_R0.size_bits());

            tunstall_seq_pct[1].push_back((float)T_R1.bits_tunstall_seq() / T_R1.size_bits());
            huffman_seq_pct[1].push_back((float)T_R1.bits_huffman_seq() / T_R1.size_bits());
            bit_vec_b_pct[1].push_back((float)T_R1.bits_bit_vec_b() / T_R1.size_bits());
            hybrid_bytes.push_back(T_R0.size_bytes() + T_R1.size_bytes());

            labels.push_back("( " + to_string(top_ks[i]) + " ; " + to_string(top_ks[j]) + " )");
            bit_vec_B.push_back("elias-fano");
            symbol_pct[0].push_back(r0_tuns_sybmols_pct);
            symbol_pct[1].push_back(curr_symbols_in_tuns);
        }
    }

    // hybrid coder with hybrid bit vector B
    for (uint64_t i = 0; i < top_ks.size(); i++) {
        hybrid_coder T_R0(bv_R0, 512, top_ks[i], curr_symbols_in_tuns, false);
        float r0_tuns_sybmols_pct = curr_symbols_in_tuns;
        for (uint64_t j = 0; j < top_ks.size(); j++) {
            rows++;
            cout << "T_R0: Top-" << top_ks[i] << " T_R1: Top-" << top_ks[j] << " HYBRID" << endl;
            
            if (print_res) {
                cout << "T_R0 info:" << endl;
                T_R0.print_info();
                cout << "T_R1 info:" << endl;
            }
            hybrid_coder T_R1(bv_R1, 512, top_ks[j], curr_symbols_in_tuns, false);
            if (print_res) {
                T_R1.print_info();
                cout << endl;
                cout << "Tunstall of top-" << top_ks[i] << " most-frequent runs of 0: " << (float)T_R0.size_bits() / posting_list_size << endl;
                cout << "Tunstall of top-" << top_ks[i] << " most-frequent runs of 1: " << (float)T_R1.size_bits() / posting_list_size << endl;
                cout << "T_R0 size from total (%): " << ((float)T_R0.size_bits()) / ((float)T_R0.size_bits() + T_R1.size_bits()) << endl;
                cout << "T_R1 size from total (%): " << ((float)T_R1.size_bits()) / ((float)T_R0.size_bits() + T_R1.size_bits()) << endl;
                cout << "Total size Tunstal: " << ((float)T_R0.size_bits() + T_R1.size_bits()) / posting_list_size << endl;
            }

            // saving results
            tunstall_bpp[0].push_back((float)T_R0.bits_tunstall_seq() / posting_list_size);
            huffman_bpp[0].push_back((float)T_R0.bits_huffman_seq() / posting_list_size);
            bit_vec_b_bpp[0].push_back((float)T_R0.bits_bit_vec_b() / posting_list_size);

            tunstall_bpp[1].push_back((float)T_R1.bits_tunstall_seq() / posting_list_size);
            huffman_bpp[1].push_back((float)T_R1.bits_huffman_seq() / posting_list_size);
            bit_vec_b_bpp[1].push_back((float)T_R1.bits_bit_vec_b() / posting_list_size);
            
            total_bpp[0].push_back( ( (float)T_R0.size_bits() ) / posting_list_size );
            total_bpp[1].push_back( ( (float)T_R1.size_bits() ) / posting_list_size );
            total_bpp[2].push_back( ( (float)T_R0.size_bits() + T_R1.size_bits() ) / posting_list_size );

            tunstall_seq_pct[0].push_back((float)T_R0.bits_tunstall_seq() / T_R0.size_bits());
            huffman_seq_pct[0].push_back((float)T_R0.bits_huffman_seq() / T_R0.size_bits());
            bit_vec_b_pct[0].push_back((float)T_R0.bits_bit_vec_b() / T_R0.size_bits());

            tunstall_seq_pct[1].push_back((float)T_R1.bits_tunstall_seq() / T_R1.size_bits());
            huffman_seq_pct[1].push_back((float)T_R1.bits_huffman_seq() / T_R1.size_bits());
            bit_vec_b_pct[1].push_back((float)T_R1.bits_bit_vec_b() / T_R1.size_bits());
            hybrid_bytes.push_back(T_R0.size_bytes() + T_R1.size_bytes());

            labels.push_back("( " + to_string(top_ks[i]) + " ; " + to_string(top_ks[j]) + " )");
            bit_vec_B.push_back("hybrid");
            symbol_pct[0].push_back(r0_tuns_sybmols_pct);
            symbol_pct[1].push_back(curr_symbols_in_tuns);
        }
    }

    

    vector< vector< vector<float> > > bpp = {tunstall_bpp, huffman_bpp, 
        bit_vec_b_bpp, total_bpp};

    vector< vector< vector<float> > > size_pct = {tunstall_seq_pct, huffman_seq_pct, bit_vec_b_pct};

    utils::separate_runs_write_results(true, "./data/prod_results/output" + curr_time + ".csv", rows,
                         labels,
                         hybrid_bytes,
                         symbol_pct,
                         size_pct,
                         bpp,
                         bit_vec_B);

    bv_R1.clear();
    bv_R0.clear();
}
}  // namespace h0run
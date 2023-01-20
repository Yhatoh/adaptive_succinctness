#include "hybrid_coder.hpp"

using namespace std;

hybrid_coder::hybrid_coder() { ; }

hybrid_coder::hybrid_coder(const vector<uint32_t> &seq, const uint32_t block_size, 
    const uint64_t top_k, float &symbol_p, bool eliasfano) {
    top_most_freq = top_k;
    map<uint32_t, uint64_t> frequence_map;
    uint64_t i;

    if (!eliasfano) eliasfano_B = false; 

    // calc frequencies of the symbols
    n = seq.size();
    for (i = 0; i < seq.size(); i++) {
        if (frequence_map.count(seq[i]) == 0)
            frequence_map[seq[i]] = 1;
        else
            frequence_map[seq[i]] += 1;
    }
    // string outputfile = "./data/test_freq_top-" + to_string(top_k) + ".txt";
    // utils::write_gapfreqs2file(outputfile, frequence_map);

    map<uint32_t, uint64_t>::iterator it;
    typedef pair<uint64_t, uint32_t> heap_node;
    priority_queue<heap_node> H;

    // build heap of the (symbols, freq) pairs
    for (it = frequence_map.begin(); it != frequence_map.end(); it++)
        H.push(heap_node(it->second, it->first));

    set<uint32_t> tunstall_alphabet;

    // create tunstall alphabet
    for (i = 0; i < top_most_freq; i++) {
        heap_node aux = H.top();
        H.pop();
        tunstall_alphabet.insert(aux.second);
    }

    vector<uint32_t> t_seq, h_seq;

    sdsl::bit_vector bv(seq.size(), 0);

    // separate symbols to tunstall and huffman sequences
    for (i = 0; i < seq.size(); i++)
        if (tunstall_alphabet.count(seq[i]) == 1)
            t_seq.push_back(seq[i]);
        else {
            h_seq.push_back(seq[i]);
            bv[i] = 1;
        }
    //cout << "Percentage of symbols in Tunstall sequence: " << (float)t_seq.size() / seq.size() << endl;
    symbol_p = (float)t_seq.size() / seq.size() ;

    tunstall_seq = tunstall_coder(t_seq, block_size, 65536);
    huffman_seq.encode(h_seq, block_size);
    if (eliasfano_B) {
        B = sdsl::sd_vector<>(bv);
    } else {
        hyb_B = sdsl::hyb_vector<>(bv);
    }
}

uint64_t hybrid_coder::bits_bit_vec_b() {
    return eliasfano_B ? 8 * sdsl::size_in_bytes(B) : 8 * sdsl::size_in_bytes(hyb_B);
}

uint64_t hybrid_coder::bytes_bit_vec_b() {
    return eliasfano_B ? sdsl::size_in_bytes(B) : sdsl::size_in_bytes(hyb_B);
}

uint64_t hybrid_coder::size_bits() {
    return tunstall_seq.size() * 8 + huffman_seq.size_bits() + bits_bit_vec_b();
}

uint64_t hybrid_coder::size_bytes() {
    return tunstall_seq.size() + huffman_seq.size_bits() / 8 + bytes_bit_vec_b();
}

uint64_t hybrid_coder::bytes_tunstall_seq() {
    return tunstall_seq.size();
}

uint64_t hybrid_coder::bytes_huffman_seq() {
    return huffman_seq.size_bits() / 8;
}

uint64_t hybrid_coder::bits_tunstall_seq() {
    return 8 * tunstall_seq.size();
}

uint64_t hybrid_coder::bits_huffman_seq() {
    return huffman_seq.size_bits();
}



void hybrid_coder::print_info() {
    cout << "Tunstall_seq size in bytes: " << bytes_tunstall_seq() << ", (% from total) "
         << (float)bytes_tunstall_seq() / size_bytes() << endl;

    cout << "Huffman_seq size in bytes: " << bytes_huffman_seq() << ", (% from total) "
         << (float)bytes_huffman_seq() / size_bytes() << endl;

    cout << "Bit vector B size in bytes: " << bytes_bit_vec_b() << ", (% from total) "
         << (float)bytes_bit_vec_b() / size_bytes() << endl;
}

// saving size results of hybrid coder
void hybrid_coder::fetch_results(vector<uint64_t> &t_seq_bytes,
                   vector<float> &t_seq_p,
                   vector<uint64_t> &h_seq_bytes,
                   vector<float> &h_seq_p,
                   vector<uint64_t> &bit_vec_bytes,
                   vector<float> &bit_vec_p,
                   vector<uint64_t> &total_bytes) {
    t_seq_bytes.push_back(bytes_tunstall_seq());
    h_seq_bytes.push_back(bytes_huffman_seq());
    bit_vec_bytes.push_back(bytes_bit_vec_b());

    t_seq_p.push_back((float)bytes_tunstall_seq() / size_bytes());
    h_seq_p.push_back((float)bytes_huffman_seq() / size_bytes());
tra.inf.santiago.usm.cl   bit_vec_p.push_back((float)bytes_bit_vec_b() / size_bytes());

    total_bytes.push_back(size_bytes());
}

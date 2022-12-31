#include "utils.hpp"

using namespace std;

uint64_t alphabet_size(vector<uint32_t> &V) {
    uint64_t i;
    map<uint32_t, uint32_t> alphabet_map;

    for (i = 0; i < V.size(); i++)
        alphabet_map[V[i]] = 0;

    return alphabet_map.size();
}

float H0(vector<uint32_t> &V, bool print) {
    uint64_t i;
    map<uint32_t, uint32_t> alphabet_map;

    for (i = 0; i < V.size(); i++)
        alphabet_map[V[i]] = 0;


    for (i = 0; i < V.size(); i++)
        alphabet_map[V[i]] = alphabet_map[V[i]] + 1;

    float H0 = 0.0;

    for (map<uint32_t, uint32_t>::iterator it = alphabet_map.begin(); it != alphabet_map.end(); it++) {

        if (print) cout << it->first << "\t" << it->second << endl;

        H0 = H0 + (float)it->second  * log2((double)V.size() / it->second);
    }

    return H0;
}

namespace utils {
// read inputfile integers to vector
void read_input_file(const char *input_file, vector<uint64_t> &vec_ref) {
    FILE *fp = fopen(input_file, "r");
    uint64_t u, n;

    fseek(fp, 0, SEEK_END);
    n = ftell(fp) / sizeof(uint64_t);
    fseek(fp, 0, SEEK_SET);

    cout << ">> Number of integers: " << n << endl;
    vec_ref.resize(n);

    for (uint64_t i = 0; i < n; i++) {
        fread(&u, sizeof(uint64_t), 1, fp);
        vec_ref[i] = u;
    }
    fclose(fp);  // closing inputfile
    cout << ">> Text already read\n";
}

void write_results(bool header, const string s,
                   const vector<string> &l,
                   const vector<uint64_t> &t_seq_bytes,
                   const vector<uint64_t> &h_seq_bytes,
                   const vector<uint64_t> &bit_vec_bytes,
                   const vector<float> &t_seq_p,
                   const vector<float> &h_seq_p,
                   const vector<float> &bit_vec_p,
                   const vector<uint64_t> &total_bytes,
                   const vector<float> &symbol_p,
                   const vector<vector<float>> &bits_per_post) {

    ofstream ofs;
    header ? ofs = ofstream(s, ios::trunc) : ofs = ofstream(s, ios::app);
    vector<string> s_vec = {"label", "tunstall seq per post size", "huffman seq per post size", "bit vector B per post size", "total size per post size",
                    "tunstall seq size", "huffman seq size", "bit vector B size",
                    "tunstall seq size (%)", "huffman seq size (%)", "bit vector B size (%)", "total byte size",
                    "% symbols in tunstall seq"};
    if (header) {
        for (int i = 0; i < s_vec.size()-1; i++) ofs << s_vec[i] << ","; 
        ofs << s_vec[s_vec.size()-1] << '\n';
    }

    for (uint64_t i = 0; i < t_seq_bytes.size(); i++) {
        ofs << l[i] << ",";
        ofs << bits_per_post[0][i] << "," << bits_per_post[1][i] << "," << bits_per_post[2][i] << "," << bits_per_post[3][i] << ",";
        ofs << t_seq_bytes[i] << ",";
        ofs << h_seq_bytes[i] << ",";
        ofs << bit_vec_bytes[i] << ",";
        ofs << t_seq_p[i] << ",";
        ofs << h_seq_p[i] << ",";
        ofs << bit_vec_p[i] << ",";
        ofs << total_bytes[i] << ",";
        ofs << symbol_p[i] << "\n";
    }
    ofs.close();
}

void separate_runs_write_results(bool header, const string s, int number_of_rows,
                   const vector<string> &l,
                   const vector<uint64_t> &total_bytes,
                   const vector< vector<float> > &symbol_p,
                   const vector< vector< vector<float> > > &pct,
                   const vector< vector< vector<float> > > &bits_per_post,
                   const vector<string> &B_labels) {

    ofstream ofs;
    header ? ofs = ofstream(s, ios::trunc) : ofs = ofstream(s, ios::app);
    vector<string> s_vec = {"TOP-K (R0 R1)", 
                            "R0-tunstall seq bpp",
                            "R0-huffman seq bpp",
                            "R0-bit vector B bpp",
                            "tunstall seq size (% R0)",
                            "huffman seq size (% R0)",
                            "bit vector B size (% R0)", 
                            "R1-tunstall seq bpp",
                            "R1-huffman seq bpp",
                            "R1-bit vector B bpp",
                            "tunstall seq size (% R1)",
                            "huffman seq size (% R1)",
                            "bit vector B size (% R1)",
                            "R0 total bpp",
                            "R1 total bpp",
                            "total bpp",
                            "total byte size",
                            "R0 % symbols in tunstall seq",
                            "R1 % symbols in tunstall seq",
                            "bit vector B"};
    if (header) {
        for (int i = 0; i < s_vec.size()-1; i++) ofs << s_vec[i] << ",";
        ofs << s_vec[s_vec.size()-1] << '\n';
    }
    
    for (int i = 0; i < number_of_rows; i++) {
        ofs << l[i] << ","; // label
        ofs << bits_per_post[0][0][i] << ","; // R0 tuntall bpp
        ofs << bits_per_post[1][0][i] << ","; // R0 huffman bpp
        ofs << bits_per_post[2][0][i] << ","; // R0 B bpp
        ofs << pct[0][0][i] << ","; // R0 tunstall seq size %
        ofs << pct[1][0][i] << ","; // R0 huffman seq size %
        ofs << pct[2][0][i] << ","; // R0 B size %
        ofs << bits_per_post[0][1][i] << ","; // R1 tuntall bpp
        ofs << bits_per_post[1][1][i] << ","; // R1 huffman bpp
        ofs << bits_per_post[2][1][i] << ","; // R1 B bpp
        ofs << pct[0][1][i] << ","; // R1 tunstall seq size %
        ofs << pct[1][1][i] << ","; // R1 huffman seq size %
        ofs << pct[2][1][i] << ","; // R1 B size %
        ofs << bits_per_post[3][0][i] << ","; // R0 bpp
        ofs << bits_per_post[3][1][i] << ","; // R1 bpp
        ofs << bits_per_post[3][2][i] << ","; // total bpp
        ofs << total_bytes[i] << ",";
        ofs << symbol_p[0][i] << ",";
        ofs << symbol_p[1][i] << ",";
        ofs << B_labels[i] << "\n";
    }
    ofs.close();
}

void write_gapfreqs2file(const string outputfile, const map<uint32_t, uint64_t> &freqs) {
    ofstream ofs = ofstream(outputfile, ios::app);
    for (auto const& p : freqs) {
        ofs << "(" << p.first << "," << p.second << ")" << ",";
    }
    ofs << "\n";
    ofs.close();
}

string current_time2str() {
    time_t now = time(0);
    tm *ltm = localtime(&now);
    string curr_time = to_string(1900 + ltm->tm_year) + "_" + to_string(1 + ltm->tm_mon) + "_" + to_string(ltm->tm_mday) +
        "_" + to_string(ltm->tm_hour) + to_string(ltm->tm_min) + to_string(ltm->tm_sec);
    return curr_time;
}
}  // namespace utils
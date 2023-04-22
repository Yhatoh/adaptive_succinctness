#ifndef _HUFFMAN_CODER_NOBLOCKCOMPRESSION_
#define _HUFFMAN_CODER_NOBLOCKCOMPRESSION_

#include <vector>
#include <map>
#include <inttypes.h>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>

#include "nqsort.h"
#include "mysort.h"

using namespace std;
using namespace sdsl;

#define L 31

#define LUT_BITS 8
#define LUT_SIZE (1 << LUT_BITS)
#define MAX_IT 0x00ffffff // ulong without top

#define MAX_ULONG 0xffffffff // maximum value for a ulong

const int32_t BUFF_BITS = sizeof(uint32_t) << 3; //# bits in a buffer element

#define MAX(a, b) ((a) < (b) ? (b) : (a))


struct blockElement {
    uint32_t prefix_sum;  // prefix sum
    uint32_t starting_position;
};

class huffman_coder {
    uint32_t cw_lens[L + 1];  // OJO cuidado con el valor de L
    /* Canonical coding arrays */
    uint32_t min_code[L];
    uint32_t lj_base[L]; 
    uint32_t offset[L];
    uint32_t max_cw_len;
    uint32_t min_cw_len; /* used to start the linear search if lut==NULL */

    uint32_t* freq_table;
    uint32_t max_symbol;    
    
    uint32_t n_symbols;
    
    uint32_t* syms;  // symbol map
    
    int32_t* lens;
   
    uint32_t* lut[LUT_SIZE]; /* canonical decode array */

    uint64_t length;  // length of the original sequence
    
    vector<uint32_t> compressed_seq;  // array with the compressed sequence
        
    uint32_t buff;  // current uint32_t being decompressed/compressed
   
    uint32_t buff_btg;
    
    uint32_t block_size;
    
    //enc_vector<> block_info_prefix_sum;
    //enc_vector<> block_info_starting_position;
    
    vector<blockElement> block_info;

    void build_freq_table(vector<uint32_t> &seq);    
    
    void OUTPUT_ULONG(uint32_t n, char len);
    
    uint32_t output(uint32_t i, uint32_t mapping[], uint32_t cwlens[]);

    void build_lut();
 
    void build_codes();

    void build_canonical_arrays(uint32_t max_cw_length);
    
    void calculate_minimum_redundancy(int32_t n);
    
    void generate_mapping(uint32_t max_cw_length, uint32_t n);

    uint32_t INPUT_ULONG(uint32_t &cur_int, int32_t len);

public:
    huffman_coder() {;};

    ~huffman_coder() 
    {
        delete [] syms;     
    };
   
    void encode(vector<uint32_t> &seq, uint32_t block_size);    
    
    uint32_t decode(uint64_t i);
    
    float size_bpe();

    uint64_t size_bits();

    uint64_t size();
    
    uint32_t sigma() {return n_symbols;};

};
#endif

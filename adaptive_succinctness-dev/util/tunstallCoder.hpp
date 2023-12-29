#ifndef _TUNSTALL_CODER_HPP_
#define _TUNSTALL_CODER_HPP_

#include <iostream> 
#include <queue> 
#include <vector>
#include <algorithm> 
#include <utility>       
#include <map>
#include <cmath>
#include <inttypes.h>
#include <fstream>
#include <string>
#include <tuple>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/suffix_arrays.hpp>

#include "block_element.hpp"


typedef std::pair<float, uint64_t>   heap_node; 
typedef std::tuple<uint64_t, uint64_t, std::string>   heap_node_str; 

typedef std::pair<uint64_t, int64_t> tree_node;
typedef std::tuple<uint64_t, int64_t, std::string> tree_node_str;

template< uint16_t w >
class tunstall_coder {
public:
    std::vector<std::vector<uint64_t>> D;  // Dictionary
    std::vector<uint32_t> map_table;
    //std::vector<uint16_t> compressed_seq;
    sdsl::int_vector<w> compressed_seq;
    std::vector<blockElement> block;
    //enc_vector<> block_info_prefix_sum;
    //enc_vector<> block_info_starting_position;
    uint32_t bSize;
    uint32_t sigma;
    //uint64_t D_size=65536 by default
    uint64_t D_size;
    
    void traverse(std::vector<tree_node>& tree, uint64_t curnode, uint64_t& curindex, uint64_t sigma,
                  std::vector<uint64_t>& currcode); 
    void traverse(std::vector<tree_node_str>& tree, uint64_t curnode, uint64_t& curindex, uint64_t sigma,
                  std::vector<uint64_t>& currcode); 

    
    tunstall_coder();
    tunstall_coder(std::vector<uint32_t> &seq, uint32_t block_size, uint64_t D_size_init);
    
    uint64_t decode(uint64_t i);  
    uint64_t dict_size(); // Tunstall dictionary size, in bytes  
    uint64_t compressed_seq_size();
    uint64_t block_vec_size();
    uint64_t size(); // compressed size, in bytes 
    uint64_t nCodewords();
};
#endif

#ifndef INCLUDE_INTSEQUENCE
#define INCLUDE_INTSEQUENCE
 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <vector>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/suffix_arrays.hpp>

#include <string.h>
#include <time.h>
#include <bits/stdc++.h>

#include "huffman_coder.hpp"

using namespace std;
using namespace sdsl;

template <class T_bv, uint32_t bSize>
class intSeqHuff
{
   uint64_t n;       // number of elements in the sequence
   uint64_t n_plus;
   uint64_t n_minus;
   uint64_t first_elem;

   huffman_coder sPlus;
   huffman_coder sMinus;
   T_bv     B;
   typename T_bv::rank_1_type B_rank; 
   
   enc_vector<> blockP;
   enc_vector<> blockM;
   uint32_t blockS;  // block size
         
public:
   intSeqHuff() {;};
	intSeqHuff(std::vector<int64_t> &X);
	~intSeqHuff() {;};
   uint64_t size() {return n;}; // number of elements in the sequence
   uint64_t sizePlus() {return n_plus;}; // number of positive differences
   uint64_t sizeMinus() {return n_minus;}; // number of negative differences
   uint64_t sigmaMinus() {return sPlus.sigma;}; // number of negative differences
   uint64_t sigmaPlus() {return sMinus.sigma;}; // number of negative differences
   
   double size_bpe(); // size of the data structure in bits per element
   
   uint32_t operator[](uint64_t i);   
};
#endif
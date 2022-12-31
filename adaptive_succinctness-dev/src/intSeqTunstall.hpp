#ifndef INCLUDE_INT_TUNSTALL_SEQUENCE
#define INCLUDE_INT_TUNSTALL_SEQUENCE
 
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

#include "tunstallCoder.hpp"

using namespace std;
using namespace sdsl;

template <class T_bv, uint32_t bSize>
class intSeqTunstall
{
   uint64_t n;       // number of elements in the sequence
   uint64_t n_plus;
   uint64_t n_minus;
   uint64_t first_elem;

   tunstall_coder<> sPlus;
   tunstall_coder<> sMinus;
   T_bv     B;
   typename T_bv::rank_1_type B_rank; 
   
   enc_vector<> blockP;
   enc_vector<> blockM;
   uint32_t blockS;  // block size
         
public:
   intSeqTunstall() {;};
	intSeqTunstall(std::vector<int64_t> &X)
	{
	    n = X.size();
       bit_vector bv(n);
       first_elem = X[0];
    
       blockS = bSize;    
       
       uint64_t i;
    
       vector<uint32_t> sP, sM;   
    
       sP.push_back(X[0]);
       bv[0] = 0;
        
       for (i = 1; i < n; ++i) {
           if (X[i] >= X[i-1]) {
               bv[i] = 0;
               sP.push_back(X[i] - X[i-1]);
           }      
           else {
               bv[i] = 1;
               sM.push_back(X[i-1] - X[i]);
           }   
       } 

       B = T_bv(bv);
       B_rank = typename T_bv::rank_1_type(&B);

       sPlus = tunstall_coder<>(sP, blockS);       
       sMinus = tunstall_coder<>(sM, blockS);
	}
	~intSeqTunstall() {;};
   uint64_t size() {return n;}; // number of elements in the sequence
   uint64_t sizePlus() {return n_plus;}; // number of positive differences
   uint64_t sizeMinus() {return n_minus;}; // number of negative differences
   uint64_t sigmaMinus() {return sPlus.sigma();}; // number of negative differences
   uint64_t sigmaPlus() {return sMinus.sigma();}; // number of negative differences
   
   double size_bpe()  // size of the data structure in bits per element
   {
       uint64_t size1 = size_in_bytes(B) + size_in_bytes(B_rank);     
       return ((double)8*(size1 + sPlus.size() + sMinus.size()))/n;
   }
   
   
   uint32_t operator[](uint64_t i)
   {
       if (i == 0) return first_elem;
        
       int64_t r1 = B_rank(i) + B[i];
       int64_t r0 = i - r1 + 1;
    
       --r0; 
       --r1;
       
       uint32_t sumP=0;
       uint32_t sumM=0;
        
       if (r0 >= 0)
           sumP = sPlus.decode(r0);
    
       if (r1 >= 0)
           sumM = sMinus.decode(r1);
    
       return sumP - sumM;
   }

};
#endif
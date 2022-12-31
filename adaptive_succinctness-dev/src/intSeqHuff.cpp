#include "intSeqHuff.hpp"
#include "huffman_coder.cpp"

template <class T_bv, uint32_t bSize>
intSeqHuff<T_bv,bSize>::intSeqHuff(std::vector<int64_t> &X)
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

    sPlus.encode(sP, blockS);       
    sMinus.encode(sM, blockS);    

}


template <class T_bv, uint32_t bSize>
double intSeqHuff<T_bv,bSize>::size_bpe()
{
 	 /*cout << endl << n_plus << " " << size_in_bytes(sPlus) << endl;
 	 cout << n_minus << " " << size_in_bytes(sMinus) << endl;
 	 cout << size_in_bytes(B) << endl;
 	 cout << size_in_bytes(B_rank) << endl;
 	 */
 	 uint64_t size1 = size_in_bytes(B) + size_in_bytes(B_rank);     
    return ((double)8*(size1 + sPlus.size() + sMinus.size()))/n;
}


template <class T_bv, uint32_t bSize>
uint32_t intSeqHuff<T_bv,bSize>::operator[](uint64_t i)
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
 

#include "intSeqHuff.cpp"
#include "intSeqTunstall.hpp"
#include "randomer.hpp"

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
using namespace std::chrono; 

int main(int argc, char** argv)
 {
    vector<int64_t> seq;
    vector<uint32_t> seq2;
    
    FILE *fp = fopen(argv[1], "r");
    unsigned u, n, sigma;
    uint64_t i;
   
    fseek(fp, 0, SEEK_END);
    n = ftell(fp)/sizeof(unsigned);    
    fseek(fp, 0, SEEK_SET);
    
    cout << n << endl;    
    
    for (i = 0; i < n; i++) {
       fread(&u, sizeof(unsigned), 1, fp);
       seq.push_back(u);   
    }
     
    cout << ">> Text already read" << endl;
      
    uint64_t max = 0;
    for (i=1; i < seq.size(); i++)
       if (seq[i] > seq[max])
          max = i;

    uint64_t *freq = new uint64_t[seq[max]+1];
      
    for (i=0; (int64_t)i < seq[max]+1; i++) 
       freq[i] = 0;   

    for (i=0; i < seq.size(); i++) 
       freq[seq[i]]++;   

    double H0_X = 0.0;
 
    sigma = 0;
    for (i=0; (int64_t)i < seq[max]+1; i++) {
       if (freq[i] != 0) {
          sigma++;
          H0_X += freq[i] * log2((double)seq.size()/freq[i]);   
       }   
    }
   
    H0_X = H0_X / seq.size();
   
    cout << "Sigma of X = " << sigma << " [" << (uint)ceil(log2(sigma)) << "]" << endl;
    cout << "H0 of X = " << H0_X << endl;

    Randomer randomer{0, n-1};
    
    vector<uint64_t> random_position;
    
    for (i=0; i < 50000; i++)
       random_position.push_back(randomer());    

   cout << "Done with random position generation" << endl;

   for (i=0; i < seq.size(); i++)
       seq2.push_back(seq[i]);    


/*   cout << "AP WT for X" << endl;

   {
       wt_ap<wt_huff<rrr_vector<15>>, wt_int<rrr_vector<15>>> X_wt;   
      
       construct(X_wt, argv[1], 4);
                  
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }    

    {
       wt_ap<wt_huff<rrr_vector<31>>, wt_int<rrr_vector<31>>> X_wt;   
      
       construct(X_wt, argv[1], 4);
                  
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }    

    {
       wt_ap<wt_huff<rrr_vector<63>>, wt_int<rrr_vector<63>>> X_wt;   
      
       construct(X_wt, argv[1], 4);
                  
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }    
    {
       wt_ap<wt_huff<rrr_vector<127>>, wt_int<rrr_vector<127>>> X_wt;   
      
       construct(X_wt, argv[1], 4);
                  
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }    

    {
       wt_ap<wt_huff<rrr_vector<255>>, wt_int<rrr_vector<255>>> X_wt;   
      
       construct(X_wt, argv[1], 4);
                  
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }    
        
    cout << "wt_int<rrr_vector<...>> for X" << endl;
    
    {
       wt_int<rrr_vector<15>> X_wt;   

       construct(X_wt, argv[1], 4);
           
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }

    {
       wt_int<rrr_vector<31>> X_wt;   

       construct(X_wt, argv[1], 4);
           
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }

    {
       wt_int<rrr_vector<63>> X_wt;   

       construct(X_wt, argv[1], 4);
           
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }

    {
       wt_int<rrr_vector<127>> X_wt;   

       construct(X_wt, argv[1], 4);
           
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }

    {
       wt_int<rrr_vector<255>> X_wt;   

       construct(X_wt, argv[1], 4);
           
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X_wt[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(X_wt)*8)/n  << " " << total_time/random_position.size() << endl;   
    }

    cout << "csa_wt<wt_rlmn" << endl;
    {
       csa_wt<wt_rlmn<rrr_vector<15>, rrr_vector<15>::rank_1_type, rrr_vector<15>::select_1_type, wt_int<rrr_vector<15>>>> csa;

       construct(csa, argv[1], 4);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;
       duration<double> time_span;
       uint64_t pos;

       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();
          extract(csa, pos, pos);
          stop = high_resolution_clock::now();

          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();

          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;

    }


 
  {
       csa_wt<wt_rlmn<rrr_vector<31>, rrr_vector<31>::rank_1_type, rrr_vector<31>::select_1_type, wt_int<rrr_vector<31>>>> csa;

       construct(csa, argv[1], 4);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;
       duration<double> time_span;
       uint64_t pos;

       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();
          extract(csa, pos, pos);
          stop = high_resolution_clock::now();

          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();

          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;

    }

  {
       csa_wt<wt_rlmn<rrr_vector<63>, rrr_vector<63>::rank_1_type, rrr_vector<63>::select_1_type, wt_int<rrr_vector<63>>>> csa;

       construct(csa, argv[1], 4);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;
       duration<double> time_span;
       uint64_t pos;

       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();
          extract(csa, pos, pos);
          stop = high_resolution_clock::now();

          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();

          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;

    }
  {
       csa_wt<wt_rlmn<rrr_vector<127>, rrr_vector<127>::rank_1_type, rrr_vector<127>::select_1_type, wt_int<rrr_vector<127>>>> csa;

       construct(csa, argv[1], 4);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;
       duration<double> time_span;
       uint64_t pos;

       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();
          extract(csa, pos, pos);
          stop = high_resolution_clock::now();

          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();

          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;

    }
  {
       csa_wt<wt_rlmn<rrr_vector<255>, rrr_vector<255>::rank_1_type, rrr_vector<255>::select_1_type, wt_int<rrr_vector<255>>>> csa;

       construct(csa, argv[1], 4);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;
       duration<double> time_span;
       uint64_t pos;

       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();
          extract(csa, pos, pos);
          stop = high_resolution_clock::now();

          time_span = duration_cast<microseconds>(stop - start);

          total_time += time_span.count();

          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;

    }



    cout << "csa_wt<wt_ap<wt_int<rrr_vector<63>>, wt_int<rrr_vector<63>>>,16384,...> for X" << endl;

    {
       csa_wt<wt_ap<wt_int<rrr_vector<63>>, wt_int<rrr_vector<63>>>,16384,2> csa;
         
       construct(csa, argv[1], 4);
      
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }

    {
       csa_wt<wt_ap<wt_int<rrr_vector<63>>, wt_int<rrr_vector<63>>>,16384,4> csa;
         
       construct(csa, argv[1], 4);
      
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }

    {
       csa_wt<wt_ap<wt_int<rrr_vector<63>>, wt_int<rrr_vector<63>>>,16384,8> csa;
         
       construct(csa, argv[1], 4);
      
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }


    {
       csa_wt<wt_ap<wt_int<rrr_vector<63>>, wt_int<rrr_vector<63>>>,16384,16> csa;
         
       construct(csa, argv[1], 4);
      
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }

    {
       csa_wt<wt_ap<wt_int<rrr_vector<63>>, wt_int<rrr_vector<63>>>,16384,32> csa;
         
       construct(csa, argv[1], 4);
      
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }

    {
       csa_wt<wt_ap<wt_int<rrr_vector<63>>, wt_int<rrr_vector<63>>>,16384,64> csa;
         
       construct(csa, argv[1], 4);
      
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }


    cout << "csa_wt<wt_int<rrr_vector<63>>,16384,...> for X " << endl;

    {
       csa_wt<wt_int<rrr_vector<63>>,4096,2> csa;
         
       construct(csa, argv[1], 4);
     
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }

    {
       csa_wt<wt_int<rrr_vector<63>>,4096,4> csa;
         
       construct(csa, argv[1], 4);
     
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }

    {
       csa_wt<wt_int<rrr_vector<63>>,4096,8> csa;
         
       construct(csa, argv[1], 4);
     
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }

    {
       csa_wt<wt_int<rrr_vector<63>>,4096,16> csa;
         
       construct(csa, argv[1], 4);
     
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }

    {
       csa_wt<wt_int<rrr_vector<63>>,4096,32> csa;
         
       construct(csa, argv[1], 4);
     
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }

    {
       csa_wt<wt_int<rrr_vector<63>>,4096,64> csa;
         
       construct(csa, argv[1], 4);
     
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          extract(csa, pos, pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i != seq[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << ((double)size_in_bytes(csa)*8)/n << " " << total_time/random_position.size() << endl;
    }
*/
   
    cout << "intSeqHuff<rrr_vector<63>, ...>" << endl;
    {
       intSeqHuff<rrr_vector<255>, 16> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq2[pos]) {
          	 cout << j  << " " << pos << " " << i << " " << seq2[pos] << endl;
             exit(1);          
          }     
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqHuff<rrr_vector<255>, 32> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqHuff<rrr_vector<255>, 64> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqHuff<rrr_vector<255>, 128> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqHuff<rrr_vector<255>, 256> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqHuff<rrr_vector<255>, 512> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqHuff<rrr_vector<255>, 1024> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    
 
    {
       intSeqHuff<rrr_vector<255>, 2048> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

/*
    cout << "huffman_coder" << endl;
            
    {
       huffman_coder X;

       X.encode(seq2, 16);
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X.decode(pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i == seq2[pos]) {
          //   cout << j  << " " << pos << " " << i << " " << seq2[pos] << endl;             
          //   exit(1);          
          //}      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       huffman_coder X;

       X.encode(seq2, 32);
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X.decode(pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i == seq2[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       huffman_coder X;

       X.encode(seq2, 64);
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X.decode(pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i == seq2[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       huffman_coder X;

       X.encode(seq2, 128);
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X.decode(pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i == seq2[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       huffman_coder X;

       X.encode(seq2, 256);
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X.decode(pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i == seq2[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       huffman_coder X;

       X.encode(seq2, 512);
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X.decode(pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i == seq2[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       huffman_coder X;

       X.encode(seq2, 1024);
       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X.decode(pos);          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          //if (i == seq2[random_position[j]]) {
          //   exit(1);          
          //}      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    */

     cout << "intSeqTunstall<rrr_vector<255>, ...>" << endl;
 
    {
       intSeqTunstall<rrr_vector<31>, 16> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                   
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqTunstall<rrr_vector<31>, 32> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    
    {
       intSeqTunstall<rrr_vector<31>, 64> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    
    {
       intSeqTunstall<rrr_vector<31>, 128> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

 
    {
       intSeqTunstall<rrr_vector<31>, 256> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqTunstall<rrr_vector<31>, 512> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqTunstall<rrr_vector<31>, 1024> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    {
       intSeqTunstall<rrr_vector<31>, 2048> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }     

    {
       intSeqTunstall<rrr_vector<31>, 4096> X(seq);

       high_resolution_clock::time_point start, stop;
       double total_time = 0.0;       
       duration<double> time_span;
       uint64_t pos; 
             
       for (uint64_t j= 0; j < random_position.size(); j++) {
          pos = random_position[j];
          start = high_resolution_clock::now();          
          i = X[pos];          
          stop = high_resolution_clock::now();          
          
          time_span = duration_cast<microseconds>(stop - start);
          //duration = duration_cast<microseconds>(stop - start);
                    
          total_time += time_span.count();
          
          //cout << j  << " " << random_position[j] << " " << i << " " << seq[random_position[j]] << endl;
          //fflush(stdout);                    
          if (i != seq[random_position[j]]) {
             exit(1);          
          }      
       }
       cout << X.size_bpe() << " " << total_time/random_position.size() << endl;
    }    

    return 0;
 
 }
 

#ifndef _H0RUN_HPP_
#define _H0RUN_HPP_

#include <sdsl/bit_vectors.hpp>
#include "hybrid_coder.hpp"
#include <vector>
#include <iostream>
#include "utils.hpp"

namespace h0run {
    void runs_together(sdsl::bit_vector &X, uint64_t posting_list_size, const std::string curr_time);
    void separate_runs(sdsl::bit_vector &X, uint64_t posting_list_size, const std::string curr_time, bool print_res);
}

#endif
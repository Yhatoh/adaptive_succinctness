#ifndef RANDOMER
#define RANDOMER

class Randomer {
  // random seed by default
  std::mt19937 gen_;
  std::uniform_int_distribution< size_t > dist_;

  public:
    
    Randomer(size_t min, size_t max, unsigned int seed = std::random_device{}()) 
      : gen_{seed}, dist_(min, max) {}

    // if you want predictable numbers 
    void SetSeed(unsigned int seed) {
      gen_.seed(seed);
    }

    size_t operator()() {
      return dist_(gen_);
    }
    
};

#endif

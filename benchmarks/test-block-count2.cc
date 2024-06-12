#define XOR_FUSE_EXPERIMENT_1

#include <iostream>
#include <iomanip>
#include <cmath>

#include "filterapi.h"
#include "random.h"
#include "timing.h"

// only 1 main can be, so it won't be included twice
struct test_params
{
  size_t start_size = 1000;
  double size_multiplier = 1.5;
  size_t size_steps = 18; // 1.5**18 = 1477.8918800354004

  size_t start_block = 10;
  size_t end_block   = 200;
  size_t block_step = 10;

  int one_test_tries = 10;
};

template <typename Filter = xorfusefilter_vanilla::XorFuseFilter<
              uint64_t, uint8_t>>
int check_block_count(size_t add_count, size_t block_count, int tries = 10) 
{
  int successes = 0;
  int failures = 0;
  for (int cur_try = 1; cur_try <= tries; cur_try++)
  {
    std::vector<uint64_t> source = GenerateRandom64(add_count);
    Filter filter(add_count, block_count);

    try
    {
      filter.AddAll(source, 0, add_count);

      successes++;
      failures += filter.hashIndex;
    }
    catch (...)
    {
      /** filter supposed to use a lot of attempts
       *  and probably it's infinite cycle
       *  so this case is bad
       *  and we can make one of 2 conclusions:
       *  1. It's just unluck, ignore it
       *  2. Probability is really about 0
      */
      failures += filter.hashIndex;
    }
    if (float(successes + (tries - cur_try)) / (successes + failures + (tries - cur_try)) < 0.9)
    {
      return 0; // bad count
    } 
  }  
  return 1; // good count
}


template <typename Filter = xorfusefilter_vanilla::XorFuseFilter<
              uint64_t, uint8_t>>
void stream_best_block_counts(const test_params &params)
{
  std::cout << "-";
  int add_count = params.start_size;

  for (size_t i = 1; i <= params.size_steps; i++)
  {
    int l = 1;
    int r = std::sqrt(add_count);
    
    int block_count = 1; // search for max block_count that generating probability is high
    // [l, r)
    while (l + 1 < r) // very accurate, do we need such??
    {
      int mid = (l + r) / 2;
      if (!check_block_count<Filter>(add_count, mid, params.one_test_tries))
      {
        r = mid; // bad count
      }
      else // good count
      {
        block_count = mid; 
        l = mid;
      }
    }
    std::cout << add_count << " " << block_count << "\n";
    add_count *= params.size_multiplier;
  }
  std::cout << "\n";
}

int main(int argc, char * argv[])
{
  test_params params;
  if (argc > 1)
  {
    params.start_size = std::atoi(argv[1]);
  }
  if (argc > 2)
  {
    params.size_steps = std::atoi(argv[2]);
  }
  if (argc > 3)
  {
    params.size_multiplier = std::atof(argv[3]);
  }
  if (argc > 4)
  {
    params.one_test_tries = std::atoi(argv[4]);
  }
  
  std::cout << "xorfusefilter_vanilla\n";
  stream_best_block_counts<xorfusefilter_vanilla::XorFuseFilter<uint64_t, uint8_t>>(params);

  std::cout << "xorfusefilter_vanilla4wise\n";
  stream_best_block_counts<xorfusefilter_vanilla4wise::XorFuseFilter<uint64_t, uint8_t>>(params);

  return EXIT_SUCCESS;
}
#define XOR_FUSE_EXPERIMENT_1

#include <iostream>
#include <iomanip>

#include "filterapi.h"
#include "random.h"
#include "timing.h"

// here are default params
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
float test_block_count(size_t add_count, size_t block_count, int tries = 10) 
{
  int successes = 0;
  int failures = 0;
  for (int cur_try = 0; cur_try < tries; cur_try++)
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
  }  
  return float(successes) / (successes + failures);
}


template <typename Filter = xorfusefilter_vanilla::XorFuseFilter<
              uint64_t, uint8_t>>
void stream_block_count(const test_params &params)
{
  std::cout << "-";
  for (size_t cur = params.start_block; cur <= params.end_block; cur += params.block_step)
  {
    std::cout << "\t" << cur;
  }
  cout << "\n";

  size_t size = params.start_size;
  for (size_t cur_iter = 1; cur_iter <= params.size_steps; cur_iter++)
  {
    std::cout << size;
    for (size_t block_count = params.start_block; block_count <= params.end_block; block_count += params.block_step)
    {
      float test_result = test_block_count<Filter>(size, block_count, params.one_test_tries);
      std::cout << "\t" << std::setprecision(3) << std::fixed << test_result;
    }
    std::cout << "\n";
    size *= params.size_multiplier;
  }
  std::cout << "\n";
}

int main()
{
  test_params params;
  // params.size_steps = 4;
  params.one_test_tries = 10;
  std::cout << "xorfusefilter_vanilla\n";
  stream_block_count<xorfusefilter_vanilla::XorFuseFilter<uint64_t, uint8_t>>(params);

  std::cout << "xorfusefilter_vanilla4wise\n";
  stream_block_count<xorfusefilter_vanilla4wise::XorFuseFilter<uint64_t, uint8_t>>(params);

  return EXIT_SUCCESS;
}
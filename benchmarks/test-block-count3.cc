#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>
#include <mutex>
#include <thread>

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


// mutex here is for output
template <typename Filter = xorfusefilter_vanilla::XorFuseFilter<
              uint64_t, uint8_t>>
void stream_best_block_counts(const test_params &params, std::mutex &mutex, const std::string name)
{
  size_t max_size = params.start_size * std::pow(params.size_multiplier, params.size_steps - 1);
  
  std::map<int, int> results; // results[size] = block_count
  // will be generated once, so its not a problem using that way of generating random
  std::vector<uint64_t> source = GenerateRandom64(max_size);

  auto check_block_count = [&](size_t add_count, int block_count, int tries)->bool
  {
    // with tries equals to 1 it looks weird
    int allowed_failures = std::max(int(tries * 0.1), 1);
    allowed_failures -= (tries == 1); // fix of problem above
  
    int failures = 0;
    for (int i = 0; i < tries; i++)
    {
      Filter filter(add_count, block_count);
      filter.AddAll(source, size_t(0), add_count);
      failures += filter.hashIndex;
      if (failures > allowed_failures)
      {
        return false;
      }
    }
    return true;
  };


  // search for first block size, want it to be more accurate
  size_t add_count = params.start_size;
  const int max_iterations = 5;
  int iterations_left = 2 * max_iterations; // only for first
  int l = 1;
  int r = std::sqrt(add_count);
  int block_count = l;
  // find first block count
  while (l + 1 < r && iterations_left--)
  {
    int mid = (l + r) / 2;
    if (!check_block_count(add_count, mid, params.one_test_tries))
    {
      r = mid; // bad count
    }
    else // good count
    {
      block_count = mid; 
      l = mid;
    }
  }
  results[add_count] = block_count;
  
  int cur_block_count = block_count;
  int prev_block_count = cur_block_count;
  for (size_t i = 2; i <= params.size_steps; i++)
  {
    iterations_left = max_iterations;
    // optimization to look less
    // suppose that next block_count <= cur_block_count * size_multiplier
    // and block_count >= cur_block_count + (cur_block_count - prev_block_count)
    l = std::max(cur_block_count + (cur_block_count - prev_block_count) - 10, 1);
    r = cur_block_count * params.size_multiplier + 10;
    
    block_count = l; // search for max block_count that generating probability is high
    // [l, r)
    while (l + 1 < r && iterations_left--)
    {
      int mid = (l + r) / 2;
      if (!check_block_count(add_count, mid, params.one_test_tries))
      {
        r = mid; // bad count
      }
      else // good count
      {
        block_count = mid; 
        l = mid;
      }
    }
    results[add_count] = block_count;

    prev_block_count = cur_block_count;
    cur_block_count = block_count;
    add_count *= params.size_multiplier;
  }

  // locking for output
  mutex.lock();
  std::cout << name << "\n-\n";
  for (const auto &[ size, block_count ] : results)
  {
    std::cout << size << " " << block_count << "\n";
  }
  std::cout << "\n";
  mutex.unlock();
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
  
  std::mutex output_mutex;
  std::thread t1(stream_best_block_counts<xorfusefilter_vanilla::XorFuseFilter<uint64_t, uint8_t>>, 
                std::ref(params), std::ref(output_mutex), "xorfusefilter_vanilla");
  std::thread t2(stream_best_block_counts<xorfusefilter_vanilla4wise::XorFuseFilter<uint64_t, uint8_t>>, 
                std::ref(params), std::ref(output_mutex), "xorfusefilter_vanilla4wise");

  t1.join();
  t2.join();

  return EXIT_SUCCESS;
}
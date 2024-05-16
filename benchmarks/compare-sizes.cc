#include "filterapi.h"
#include "random.h"
#include "timing.h"


void stream_size_benchmark(size_t start_size, size_t max_size,
                      size_t iteration = 100) {
  std::vector<uint64_t> source = GenerateRandom64Fast(max_size, rand());
  std::cout << "n, fuse, fuse4wise, binfuse,  binfuse4wise" << std::endl;
  std::cout << std::setprecision(5) << std::fixed;

  using Filter1 = xorfusefilter_vanilla::XorFuseFilter<uint64_t, uint8_t>;
  using Filter2 = xorfusefilter_vanilla4wise::XorFuseFilter<uint64_t, uint8_t>;
  using Filter3 = xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<uint64_t, uint8_t>;
  using Filter4 = xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<uint64_t, uint8_t>;

  double gap = exp(log(max_size / start_size) / iteration);
  for (double n = start_size; n <= max_size; n *= gap) {
    std::cout << size_t(n) << "\t";

    Filter1 filter1 = FilterAPI<Filter1>::ConstructFromAddCount(n);
    std::cout << double(filter1.SizeInBytes()) * 8.0 / n;
    std::cout.flush();
    std::cout << "\t";

    Filter2 filter2 = FilterAPI<Filter2>::ConstructFromAddCount(n);
    std::cout << double(filter2.SizeInBytes()) * 8.0 / n;
    std::cout.flush();
    std::cout << "\t";

    Filter3 filter3 = FilterAPI<Filter3>::ConstructFromAddCount(n);
    std::cout << double(filter3.SizeInBytes()) * 8.0 / n;
    std::cout.flush();
    std::cout << "\t";

    Filter4 filter4 = FilterAPI<Filter4>::ConstructFromAddCount(n);
    std::cout << double(filter4.SizeInBytes()) * 8.0 / n;
    std::cout.flush();

    std::cout << std::endl;
  }
}
int main() {
  stream_size_benchmark(100'000, 100'000'000, 200);
  return EXIT_SUCCESS;
}
#ifndef THREEWISE_XOR_FUSE_FILTER_XOR_FILTER_LOWMEM_H_
#define THREEWISE_XOR_FUSE_FILTER_XOR_FILTER_LOWMEM_H_

#include "xor_fuse_filter.h"

/**
 * As of July 2021, the lowmem versions of the binary fuse filters are
 * the recommended defaults.
 */

/*
 * Modified in May 2024
 */
namespace xorfusefilter_lowmem {
// status returned by a xor filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

__attribute__((always_inline)) inline uint32_t reduce(uint32_t hash,
                                                      uint32_t n) {
  // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
  return (uint32_t)(((uint64_t)hash * n) >> 32);
}

__attribute__((always_inline)) inline uint8_t mod3(uint8_t x) {
    if (x > 2) {
        x -= 3;
    }
    return x;
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily = SimpleMixSplit>
class XorFuseFilter {
public:
  size_t size;
  size_t arrayLength;
  size_t segmentCount;
  size_t segmentCountLength;
  size_t segmentLength;
  size_t segmentLengthMask;
  static constexpr size_t arity = 3;
  /* Taken from
   * Martin Dietzfelbinger and Stefan Walzer. 2019.
   * Dense Peelable Random Uniform Hypergraphs.
   * In 27th Annual European Symposium on Algorithms (ESA 2019)
   */
  static constexpr double dencity = 0.915; // threshold id 0.917935
  FingerprintType *fingerprints;

  HashFamily *hasher;
  size_t hashIndex{0};

  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType)hash ^ (hash >> 32);
  }

  explicit XorFuseFilter(const size_t size) {
    hasher = new HashFamily();
    this->size = size;
    // max supported segment size is 2**18
    // but it's not gonna be reached if you use size <= 3 * 10^8
    // work stability at sizes > 10^8 wasn't checked
    // you can try to choose different formula
    this->segmentCount = 0.216 * pow(size, 0.44);
    this->segmentCount = std::max(size_t(1), segmentCount);

    double sizeFactor = (this->segmentCount + arity - 1) 
                            / (this->segmentCount * this->dencity);
    size_t capacity = size * sizeFactor;
    this->segmentLength = (capacity + this->segmentCount - 1) // it's not bug
                          / this->segmentCount;
                            
    if (sizeFactor > 1.23) // xorfilter
    {
      this->segmentCount = 1;
      this->segmentLength = (size * 1.23 + 32) / 3; 
    }

    this->arrayLength = (this->segmentCount + arity - 1) * this->segmentLength;
    this->segmentCountLength = this->segmentCount * this->segmentLength;
    fingerprints = new FingerprintType[arrayLength]();
    std::fill_n(fingerprints, arrayLength, 0);
  }

  ~XorFuseFilter() {
    delete[] fingerprints;
    delete hasher;
  }

  Status AddAll(const vector<ItemType> &data, const size_t start,
                const size_t end) {
    return AddAll(data.data(), start, end);
  }

  Status AddAll(const ItemType *data, const size_t start, const size_t end);

  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  /* methods for providing stats  */
  // summary infomation
  std::string Info() const;

  // number of current inserted items;
  size_t Size() const { return size; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return arrayLength * sizeof(FingerprintType); }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status XorFuseFilter<ItemType, FingerprintType, HashFamily>::AddAll(
    const ItemType *keys, const size_t start, const size_t end) {

  uint64_t *reverseOrder = new uint64_t[size+1];
  uint8_t *reverseH = new uint8_t[size];
  size_t reverseOrderPos;

  // the lowest 2 bits are the h index (0, 1, or 2)
  // so we only have 6 bits for counting;
  // but that's sufficient
  uint8_t *t2count = new uint8_t[arrayLength];
  
  uint64_t *t2hash = new uint64_t[arrayLength];

  size_t *alone = new size_t[arrayLength];
  hashIndex = 0;
  
  // the array h0, h1, h2, h0, h1, h2
  size_t h012[5];
  while (true) {
    memset(t2count, 0, sizeof(uint8_t) * arrayLength);
    memset(t2hash, 0, sizeof(uint64_t) * arrayLength);

    // counting sort

    memset(reverseOrder, 0, sizeof(uint64_t) * size);
    reverseOrder[size] = 1;

    int blockBits = 1;
    while((size_t(1)<<blockBits) < segmentCount) { blockBits++; }
    size_t block = size_t(1) << blockBits;

    size_t *startPos = new size_t[block];
    for(uint32_t i = 0; i < block; i++) { startPos[i] = i * size / block; }
    for (size_t i = start; i < end; i++) {
      uint64_t k = keys[i];
      uint64_t hash = (*hasher)(k);
      size_t segment_index = hash >> (64 - blockBits);
      while(reverseOrder[startPos[segment_index]] != 0) {
        segment_index++;
        segment_index &= (size_t(1) << blockBits) - 1;
      }
      reverseOrder[startPos[segment_index]] = hash;
      startPos[segment_index]++;
    }
    uint8_t countMask = 0;
    for (size_t i = 0; i < size; i++) {
      uint64_t hash = reverseOrder[i];
      __uint128_t x = (__uint128_t)hash * (__uint128_t)segmentCountLength;
      h012[0] = (uint64_t)(x >> 64);
      assert(h012[0] < segmentCountLength);
      size_t index = h012[0] / segmentLength;
      
      h012[1] = (hash >> 18) & ((1 << 18) - 1);
      h012[1] = ((uint64_t)h012[1] * segmentLength) >> 18; // apply reduce
      // assert(h012[1] < segmentLength);
      h012[1] += (index + 1) * segmentLength;
      
      h012[2] = hash & ((1 << 18) - 1);
      h012[2] = ((uint64_t)h012[2] * segmentLength) >> 18; // apply reduce
      // assert(h012[2] < segmentLength);
      h012[2] += (index + 2) * segmentLength;
      
      for (int hi = 0; hi < 3; hi++) {
        index = h012[hi];
        t2count[index] += 4;
        t2count[index] ^= hi;
        t2hash[index] ^= hash;
        countMask |= t2count[index];
      }
    }
    delete[] startPos;
    if (countMask >= 0x80) {
      // we have a possible counter overflow
      // this branch is never taken except if there is a problem in the hash code
      // in which case construction fails
      memset(fingerprints, ~0, arrayLength * sizeof(FingerprintType));
      return Ok;
    }

    reverseOrderPos = 0;
    size_t alonePos = 0;
    for (size_t i = 0; i < arrayLength; i++) {
      alone[alonePos] = i;
      int inc = (t2count[i] >> 2) == 1 ? 1 : 0;
      alonePos += inc;
    }

    while (alonePos > 0) {
      alonePos--;
      size_t index = alone[alonePos];
      if ((t2count[index] >> 2) == 1) {
        // It is still there!
        uint64_t hash = t2hash[index];
        int found = t2count[index] & 3;
        
        reverseH[reverseOrderPos] = found;
        reverseOrder[reverseOrderPos] = hash;
        
        __uint128_t x = (__uint128_t)hash * (__uint128_t)segmentCountLength;
        h012[0] = (uint64_t)(x >> 64);
        index = h012[0] / segmentLength;
        
        h012[1] = (hash >> 18) & ((1 << 18) - 1);
        h012[1] = ((uint64_t)h012[1] * segmentLength) >> 18; // apply reduce
        // assert(h012[1] < segmentLength);
        h012[1] += (index + 1) * segmentLength;

        h012[2] = hash & ((1 << 18) - 1);
        h012[2] = ((uint64_t)h012[2] * segmentLength) >> 18; // apply reduce
        // assert(h012[2] < segmentLength);
        h012[2] += (index + 2) * segmentLength;

        size_t index3 = h012[mod3(found + 1)];
        alone[alonePos] = index3;
        alonePos += ((t2count[index3] >> 2) == 2 ? 1 : 0);
        t2count[index3] -= 4;
        t2count[index3] ^= mod3(found + 1);
        t2hash[index3] ^= hash;

        index3 = h012[mod3(found + 2)];
        alone[alonePos] = index3;
        alonePos += ((t2count[index3] >> 2) == 2 ? 1 : 0);
        t2count[index3] -= 4;
        t2count[index3] ^= mod3(found + 2);
        t2hash[index3] ^= hash;

        reverseOrderPos++;
      }
    }

    if (reverseOrderPos == size) {
      break;
    }
    hashIndex++;

    // use a new random numbers
    delete hasher;
    hasher = new HashFamily();
  }
  delete[] alone;
  delete[] t2count;
  delete[] t2hash;
  
  for (int i = reverseOrderPos - 1; i >= 0; i--) {
    uint64_t hash = reverseOrder[i];
    int found = reverseH[i];
    FingerprintType xor2 = fingerprint(hash);

    __uint128_t x = (__uint128_t)hash * (__uint128_t)segmentCountLength;
    h012[0] = (uint64_t)(x >> 64);
    size_t index = h012[0] / segmentLength;
    
    h012[1] = (hash >> 18) & ((1 << 18) - 1);
    h012[1] = ((uint64_t)h012[1] * segmentLength) >> 18; // apply reduce
    h012[1] += (index + 1) * segmentLength;
    
    h012[2] = hash & ((1 << 18) - 1);
    h012[2] = ((uint64_t)h012[2] * segmentLength) >> 18; // apply reduce
    h012[2] += (index + 2) * segmentLength;

    h012[3] = h012[0];
    h012[4] = h012[1];
    fingerprints[h012[found]] = xor2 ^ fingerprints[h012[found + 1]] ^ fingerprints[h012[found + 2]];
  }
  delete[] reverseOrder;
  delete[] reverseH;

  return Ok;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status XorFuseFilter<ItemType, FingerprintType, HashFamily>::Contain(
    const ItemType &key) const {
  uint64_t hash = (*hasher)(key);
  FingerprintType f = fingerprint(hash);
  __uint128_t x = (__uint128_t)hash * (__uint128_t)segmentCountLength;
  
  uint64_t hh = hash;
  uint64_t h0 = (uint64_t)(x >> 64);
  int index = h0 / segmentLength;
  uint64_t h1 = (hh >> 18) & ((1 << 18) - 1);
  uint64_t h2 = hh & ((1 << 18) - 1);

  h1 = (h1 * segmentLength) >> 18; // apply reduce
  h1 += (index + 1) * segmentLength;
  h2 = (h2 * segmentLength) >> 18; // apply reduce
  h2 += (index + 2) * segmentLength;

  f ^= fingerprints[h0] ^ fingerprints[h1] ^ fingerprints[h2];
  return f == 0 ? Ok : NotFound;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
std::string
XorFuseFilter<ItemType, FingerprintType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "XorFuseFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  return ss.str();
}
} // namespace xorfusefilter_lowmem
#endif
#ifndef FOURWISE_XOR_FUSE_FILTER_XOR_FILTER_LOWMEM_H_
#define FOURWISE_XOR_FUSE_FILTER_XOR_FILTER_LOWMEM_H_
#include "xor_fuse_filter.h"
/**
 * As of July 2021, the lowmem versions of the binary fuse filters are
 * the recommended defaults.
 */

/*
 * Modified in May 2024 
 */
namespace xorfusefilter_lowmem4wise {
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

__attribute__((always_inline)) inline uint64_t reduce(uint64_t hash,
                                                      uint32_t n)
{
  // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
  return (uint64_t)(((__uint128_t)hash * n) >> 64);
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
  static constexpr size_t arity = 4;
  /* Taken from
   * Martin Dietzfelbinger and Stefan Walzer. 2019.
   * Dense Peelable Random Uniform Hypergraphs.
   * In 27th Annual European Symposium on Algorithms (ESA 2019)
   */
  static constexpr double dencity = 0.96; // threshold is 0.97677
  FingerprintType *fingerprints;

  HashFamily *hasher;
  size_t hashIndex{0};

  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType)hash ^ (hash >> 32);
  }

  static inline __attribute__((always_inline)) uint64_t rotateLeft(uint64_t n,
                                                               unsigned int c) {
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    c &= mask;
    return (n << c) | (n >> ((-c) & mask));
  }

  static inline __attribute__((always_inline)) uint64_t rotateRight(uint64_t n,
                                                               unsigned int c) {
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    c &= mask;
    return (n >> c) | (n << ((-c) & mask));
  }

  inline __attribute__((always_inline)) size_t
  getHashFromHash(uint64_t hash, int index) const {
    __uint128_t x = (__uint128_t)hash * (__uint128_t)segmentCountLength;
    uint64_t h = (uint64_t)(x >> 64);
    if (index > 0)
    {
      uint64_t hh = h / segmentLength;
      // h = (hash >> ((index - 1) * 16)) & 65535;
      // h = (h * segmentLength) >> 16; // apply reduce
      h = rotateLeft(hash, index * 16);
      h = reduce(h, segmentLength);
      h += (hh + index) * segmentLength;
    }
    return h;
  }

  explicit XorFuseFilter(const size_t size) {
    hasher = new HashFamily();
    this->size = size;

    // constuction stability wasnt checked at size > 10^8
    this->segmentCount = 0.003 * pow(size, 0.77);
    this->segmentCount = std::max(size_t(1), segmentCount);

    double sizeFactor = (this->segmentCount + arity - 1) 
                  / (this->segmentCount * this->dencity);
    size_t capacity = size * sizeFactor;
    this->segmentLength = (capacity + this->segmentCount - 1) 
                            / this->segmentCount;

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
  
  // the array h0, h1, h2, h3, h0, h1, h2, h3
  size_t h0123[7];
  
  size_t hi0123[7];
  hi0123[0] = 0;
  hi0123[1] = 1;
  hi0123[2] = 2;
  hi0123[3] = 3;
  hi0123[4] = 0;
  hi0123[5] = 1;
  hi0123[6] = 2;

  while (true) {
    // cout << "attempt\n";
    memset(t2count, 0, sizeof(uint8_t) * arrayLength);
    memset(t2hash, 0, sizeof(uint64_t) * arrayLength);

    // counting sort

    memset(reverseOrder, 0, sizeof(uint64_t) * size);
    reverseOrder[size] = 1;

    int blockBits = 1;
    while((size_t(1)<<blockBits) < segmentCount) { blockBits++; }
    size_t block = size_t(1) << blockBits;
    size_t *startPos = new size_t[block];
    for(uint32_t i = 0; i < uint32_t(1) << blockBits; i++) { startPos[i] = i * size / block; }
    for (size_t i = start; i < end; i++) {
      uint64_t k = keys[i];
      uint64_t hash = (*hasher)(k);
      size_t segment_index = hash >> (64 - blockBits);
      // We only overwrite when the hash was zero. Zero hash values
      // may be misplaced (unlikely).
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
      for (int hi = 0; hi < 4; hi++) {
        int index = getHashFromHash(hash, hi);
        t2count[index] += 4;
        t2count[index] ^= hi;
        t2hash[index] ^= hash;
        // this branch is never taken except if there is a problem in the hash code
        // in which case construction would fail
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
        
        h0123[1] = getHashFromHash(hash, 1);
        h0123[2] = getHashFromHash(hash, 2);
        h0123[3] = getHashFromHash(hash, 3);
        h0123[4] = getHashFromHash(hash, 0);
        h0123[5] = h0123[1];
        h0123[6] = h0123[2];

        size_t index3 = h0123[found + 1];
        alone[alonePos] = index3;
        alonePos += ((t2count[index3] >> 2) == 2 ? 1 : 0);
        t2count[index3] -= 4;
        t2count[index3] ^= hi0123[found + 1];
        t2hash[index3] ^= hash;

        index3 = h0123[found + 2];
        alone[alonePos] = index3;
        alonePos += ((t2count[index3] >> 2) == 2 ? 1 : 0);
        t2count[index3] -= 4;
        t2count[index3] ^= hi0123[found + 2];
        t2hash[index3] ^= hash;

        index3 = h0123[found + 3];
        alone[alonePos] = index3;
        alonePos += ((t2count[index3] >> 2) == 2 ? 1 : 0);
        t2count[index3] -= 4;
        t2count[index3] ^= hi0123[found + 3];
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
    h0123[0] = getHashFromHash(hash, 0);
    h0123[1] = getHashFromHash(hash, 1);
    h0123[2] = getHashFromHash(hash, 2);
    h0123[3] = getHashFromHash(hash, 3);
    h0123[4] = h0123[0];
    h0123[5] = h0123[1];
    h0123[6] = h0123[2];
    fingerprints[h0123[found]] = xor2 ^ fingerprints[h0123[found + 1]] ^ fingerprints[h0123[found + 2]] ^ fingerprints[h0123[found + 3]];
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
  // Is this thing.
  // uint64_t h1 = hash;
  // uint64_t h2 = (hash >> 16) & 65535;
  // uint64_t h3 = (hash >> 32) & 65535;
  uint64_t h1 = rotateLeft(hash, 16);
  uint64_t h2 = rotateLeft(hash, 32);
  uint64_t h3 = rotateLeft(hash, 48);

  __uint128_t x = (__uint128_t)hash * (__uint128_t)segmentCountLength;
  uint64_t h0 = (uint64_t)(x >> 64);
    
  uint64_t hh = h0 / segmentLength;
  // h1 = (h1 * segmentLength) >> 16;
  // h2 = (h2 * segmentLength) >> 16;
  // h3 = (h3 * segmentLength) >> 16;
  h1 = reduce(h1, segmentLength);
  h2 = reduce(h2, segmentLength);
  h3 = reduce(h3, segmentLength);

  h1 += (hh + 1) * segmentLength; 
  h2 += (hh + 2) * segmentLength; 
  h3 += (hh + 3) * segmentLength; 

  f ^= fingerprints[h0] ^ fingerprints[h1] ^ fingerprints[h2] ^ fingerprints[h3];
  // FingerprintType f = fingerprint(hash);
  // for (int hi = 0; hi < 4; hi++) {
  //   size_t h = getHashFromHash(hash, hi);
  //   f ^= fingerprints[h];
  // }
  return f == 0 ? Ok : NotFound;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
std::string
XorFuseFilter<ItemType, FingerprintType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "4-wise XorFuseFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  return ss.str();
}
} // namespace xorfusefilter_lowmem4wise
#endif
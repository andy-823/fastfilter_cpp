#ifndef THREEWISE_XOR_FUSE_FILTER_XOR_FILTER_VANILLA_H_
#define THREEWISE_XOR_FUSE_FILTER_XOR_FILTER_VANILLA_H_
#include "xor_binary_fuse_filter.h"

/**
 * This is vanilla verion of xor fuse filter
 * It is supposed to be used in theoretical analisys
 * where speed is not the issue
 *  
 * It's much slower and takes more space during construction
 * than 3wise_xor_binary_fuse_filter_lowmem 
 */

namespace xorfusefilter_vanilla
{
// status returned by a xor filter operation
enum Status
{
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

inline uint64_t rotl64(uint64_t n, unsigned int c) {
  // assumes width is a power of 2
  const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
  // assert ( (c<=mask) &&"rotate by type width or more");
  c &= mask;
  return (n << c) | (n >> ((-c) & mask));
}

// for 32-bit integers
__attribute__((always_inline)) inline uint32_t reduce(uint32_t hash,
                                                      uint32_t n)
{
  // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
  return (uint32_t)(((uint64_t)hash * n) >> 32);
}

// for 64-bit integers
__attribute__((always_inline)) inline uint64_t reduce64(uint64_t hash,
                                                      uint32_t n)
{
  // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
  return (uint64_t)(((__uint128_t)hash * n) >> 64);
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily = SimpleMixSplit>
class XorFuseFilter
{
 public:
  size_t size;
  size_t arrayLength;
  size_t segmentCount;
  size_t segmentCountLength;
  size_t segmentLength;
  static constexpr size_t arity = 3;
  /* Taken from
   * Martin Dietzfelbinger and Stefan Walzer. 2019.
   * Dense Peelable Random Uniform Hypergraphs.
   * In 27th Annual European Symposium on Algorithms (ESA 2019)
   */
  static constexpr double dencity = 0.915; // threshold id 0.917935

  FingerprintType *fingerprints;

  // hash functions are supposed to be different every time they are constructed
  HashFamily *h;

  size_t hashIndex{0};

  inline FingerprintType fingerprint(const uint64_t hash) const
  {
    return (FingerprintType)hash;
  }

  explicit XorFuseFilter(const size_t size, int segmentCount_ = -1)
  // explicit XorFuseFilter(const size_t size, const int segmentCount_ = 100)
  {
    h = new HashFamily[3]();
    this->size = size;

    if (segmentCount_ > 0) // segmentCount has been set
    {
      this->segmentCount = segmentCount_;
    }
    else // default
    {
      this->segmentCount = 0.216 * pow(size, 0.44);
      this->segmentCount = std::max(size_t(1), segmentCount);
    }

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

  ~XorFuseFilter()
  {
    delete[] fingerprints;
    delete[] h;
  }

  Status AddAll(const vector<ItemType> &data, const size_t start,
                const size_t end)
  {
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
__attribute__((noinline))
Status XorFuseFilter<ItemType, FingerprintType, HashFamily>::AddAll(
    const ItemType *keys, const size_t start, const size_t end)
{
  // that means only sizes < 2 * 10**9 avaliable
  // for greater sizes use uint64_t
  struct edge
  {
    uint32_t i[3];
    FingerprintType fingerprint;
  };
  
  edge *reverseOrder = new edge[size];
  uint8_t *reverseH = new uint8_t[size];  // hash index (0, 1, or 2)
  size_t reverseOrderPos;

  // the lowest 2 bits are the h index (0, 1, or 2)
  // so we only have 14 bits for counting;
  // but that's sufficient
  uint16_t *edgeCount = new uint16_t[arrayLength];
  edge *edgeXor = new edge[arrayLength];

  size_t *alone = new size_t[arrayLength];
  hashIndex = 0;

  // size_t max_iter = 1;
  // while (max_iter--)  
  while (true)
  {
    memset(edgeCount, 0, sizeof(uint16_t) * arrayLength);
    memset(edgeXor, 0, sizeof(edge) * arrayLength);
    memset(reverseOrder, 0, sizeof(uint64_t) * size);

    edge curEdge;
    uint16_t countMask = 0;
    uint32_t index;
    uint64_t hash[3];

    #ifndef XOR_FUSE_FILTER_IGNORE_SORT
    int blockBits = 1;
    while((size_t(1)<<blockBits) < segmentCount)
    { 
      blockBits++;
    }
    size_t blockCnt = size_t(1) << blockBits;
    size_t startBlock = 0;
    size_t endBlock;
    size_t *startPos = new size_t[blockCnt];  // want try _malloca(block * sizeof(size_t));
    size_t *threshold = new size_t[blockCnt];  
    for(uint32_t i = 0; i < blockCnt; i++)
    { 
      endBlock = (i + 1) * size / blockCnt;
      startPos[i] = startBlock;
      threshold[i] = endBlock;
      startBlock = endBlock;
    }

    blockCnt--; // now its a mask
    for (size_t i = start; i < end; i++)
    {
      hash[0] = h[0](keys[i]);
      hash[1] = h[1](keys[i]);
      hash[2] = h[2](keys[i]);
      
      size_t blockIndex = hash[0] >> (64 - blockBits);
      while (startPos[blockIndex] == threshold[blockIndex])
      {
        blockIndex++;
        blockIndex &= blockCnt; // run in loop
      }

      curEdge.fingerprint = fingerprint(hash[0]);
      hash[0] = reduce64(hash[0], segmentCountLength);
      hash[1] = reduce64(hash[1], segmentLength);
      hash[2] = reduce64(hash[2], segmentLength);

      index = hash[0] / segmentLength;
      curEdge.i[0] = hash[0];
      curEdge.i[1] = (index + 1) * this->segmentLength + hash[1];
      curEdge.i[2] = (index + 2) * this->segmentLength + hash[2];

      reverseOrder[startPos[blockIndex]] = curEdge;
      startPos[blockIndex]++; 
    }
    delete[] startPos;
    delete[] threshold;
    for (size_t i = 0; i < size; ++i)
    {
      curEdge = reverseOrder[i];
      for (int hi = 0; hi < 3; ++hi)
      {
        index = curEdge.i[hi];
        edgeCount[index] += 4;
        edgeCount[index] ^= hi; // hash index

        edgeXor[index].i[0] ^= curEdge.i[0];
        edgeXor[index].i[1] ^= curEdge.i[1];
        edgeXor[index].i[2] ^= curEdge.i[2];
        edgeXor[index].fingerprint ^= curEdge.fingerprint;

        countMask |= edgeCount[index];
      }
    }
    #else // XOR_FUSE_FILTER_IGNORE_SORT is defined
    for (size_t i = start; i < end; ++i)
    {
      hash[0] = h[0](keys[i]);
      hash[1] = h[1](keys[i]);
      hash[2] = h[2](keys[i]);

      curEdge.fingerprint = fingerprint(hash[0]);
    
      hash[0] = reduce64(hash[0], segmentCountLength);
      hash[1] = reduce64(hash[1], segmentLength);
      hash[2] = reduce64(hash[2], segmentLength);

      index = hash[0] / segmentLength;
      curEdge.i[0] = hash[0];
      curEdge.i[1] = (index + 1) * this->segmentLength + hash[1];
      curEdge.i[2] = (index + 2) * this->segmentLength + hash[2];

      for (int hi = 0; hi < 3; hi++)
      {
        index = curEdge.i[hi];
        edgeCount[index] += 4;
        edgeCount[index] ^= hi; // hash index

        edgeXor[index].i[0] ^= curEdge.i[0];
        edgeXor[index].i[1] ^= curEdge.i[1];
        edgeXor[index].i[2] ^= curEdge.i[2];
        edgeXor[index].fingerprint ^= curEdge.fingerprint;

        countMask |= edgeCount[index];
      }
    }
    #endif

    if (countMask >= 65535)
    {
      // we have a possible counter overflow
      // this branch is never taken except if there is a problem in the hash code
      // in which case construction fails
      memset(fingerprints, ~0, arrayLength * sizeof(FingerprintType));
      return Ok;
    }

    reverseOrderPos = 0;
    size_t alonePos = 0;
    for (size_t i = 0; i < arrayLength; i++)
    {
      alone[alonePos] = i;
      alonePos += (edgeCount[i] >> 2) == 1;
    }

    // constexpr uint16_t hi012[5] = { 0, 1, 2, 0, 1 };
    uint16_t hi012[5] = { 0, 1, 2, 0, 1 };
    while (alonePos > 0)
    {
      alonePos--;
      size_t index = alone[alonePos];
      if ((edgeCount[index] >> 2) == 1)
      {
        // It is still there!
        curEdge = edgeXor[index];
        int found = edgeCount[index] & 3;

        reverseH[reverseOrderPos] = found;
        reverseOrder[reverseOrderPos] = curEdge;


        size_t index3 = curEdge.i[hi012[found + 1]];
        alone[alonePos] = index3;
        alonePos += (edgeCount[index3] >> 2) == 2;
        edgeCount[index3] -= 4;
        edgeCount[index3] ^= hi012[found + 1]; // hash index

        edgeXor[index3].i[0] ^= curEdge.i[0];
        edgeXor[index3].i[1] ^= curEdge.i[1];
        edgeXor[index3].i[2] ^= curEdge.i[2];
        edgeXor[index3].fingerprint ^= curEdge.fingerprint;


        index3 = curEdge.i[hi012[found + 2]];
        alone[alonePos] = index3;
        alonePos += (edgeCount[index3] >> 2) == 2;
        edgeCount[index3] -= 4;
        edgeCount[index3] ^= hi012[found + 2]; // hash index

        edgeXor[index3].i[0] ^= curEdge.i[0];
        edgeXor[index3].i[1] ^= curEdge.i[1];
        edgeXor[index3].i[2] ^= curEdge.i[2];
        edgeXor[index3].fingerprint ^= curEdge.fingerprint;


        ++reverseOrderPos;
      }
    }
    if (reverseOrderPos == size)
    {
      break;
    }
    ++hashIndex;
    // use a new random numbers
    delete[] h;
    h = new HashFamily[3]();
  }
  delete[] alone;
  delete[] edgeXor;
  delete[] edgeCount;

  // the array h0, h1, h2, h0, h1, h2
  uint32_t h012[5];
  edge *curEdge = &reverseOrder[size];
  for (int i = reverseOrderPos - 1; i >= 0; --i)
  {
    curEdge--;
    int found = reverseH[i];

    h012[0] = curEdge->i[0];
    h012[1] = curEdge->i[1];
    h012[2] = curEdge->i[2];
    h012[3] = h012[0];
    h012[4] = h012[1];

    fingerprints[h012[found]] = curEdge->fingerprint ^ fingerprints[h012[found + 1]] ^ fingerprints[h012[found + 2]];
  }

  delete[] reverseOrder;
  delete[] reverseH;

  return Ok;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status XorFuseFilter<ItemType, FingerprintType, HashFamily>::Contain(
    const ItemType &key) const
{
  uint64_t h0 = h[0](key);
  uint64_t h1 = h[1](key);
  uint64_t h2 = h[2](key);
  
  FingerprintType f = fingerprint(h0);

  h0 = reduce64(h0, segmentCountLength);
  h1 = reduce64(h1, segmentLength);
  h2 = reduce64(h2, segmentLength);

  uint32_t ib = h0 / segmentLength;
  h1 += (ib + 1) * segmentLength;
  h2 += (ib + 2) * segmentLength;

  f ^= fingerprints[h0] ^ fingerprints[h1] ^ fingerprints[h2];
  return f == 0 ? Ok : NotFound;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
std::string
XorFuseFilter<ItemType, FingerprintType, HashFamily>::Info() const
{
  std::stringstream ss;
  ss << "XorFuseFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  return ss.str();
}

} // namespace xorfusefilter_vanilla
#endif // THREEWISE_XOR_FUSE_FILTER_XOR_FILTER_VANILLA_H_
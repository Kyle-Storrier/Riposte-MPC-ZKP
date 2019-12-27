#include <bitset>
#include <vector>
#include <string>
#include<iostream>
#include <x86intrin.h>

static constexpr unsigned numofboxes = 63; // Number of Sboxes
static constexpr unsigned blocksize = 256; // Block size in bits
static constexpr unsigned keysize = 128; // Key size in bits
static constexpr unsigned rounds = 14; // Number of rounds
static constexpr unsigned identitysize = blocksize - 3*numofboxes;

static const __m128i lsb128_mask[2] = {
   _mm_setzero_si128(),                                     // 0b00...0000
   _mm_set_epi64x(0,1)                                      // 0b00...0001
};
static const __m128i lsb128_mask_inv = _mm_set_epi64x(-1,-2); // 0b11...1110
static const __m128i if128_mask[2] = {
   _mm_setzero_si128(),                                     // 0b00...0000
   _mm_set1_epi8(-1)                                        // 0b11...1111
 };

static const __m256i lsb256_mask[2] = {
   _mm256_setzero_si256(),                                     // 0b00...0000
   _mm256_set_epi64x(0,0,0,1)                                      // 0b00...0001
};
static const __m256i lsb256_mask_inv = _mm256_set_epi64x(-1,-1,-1,-2); // 0b11...1110
static const __m256i if256_mask[2] = {
   _mm256_setzero_si256(),                                     // 0b00...0000
   _mm256_set1_epi8(-1)                                        // 0b11...1111
};

template <typename __mX>
union blocks
{
 public:
  blocks(const unsigned long long input = 0ULL) : bits(input) { }
  blocks(const __mX & val) : mX(val) { }
  blocks(std::string bit_string) : bits(bit_string) { }
  inline operator __mX() const { return mX; }
  inline blocks<__mX> & operator=(const __mX & val) { mX = val; return *this; }
  //inline std::bitset<blocksize>::reference operator[](const size_t pos) { return bits[pos]; }
  inline const bool operator[] (const size_t pos) const { return bits[pos]; }
  constexpr inline size_t size() const { return 8 * sizeof(__mX); }
  inline const unsigned parity() const { return bits.count() % 2; } 
  inline const blocks<__mX> shiftr_bits(const size_t pos) const { return bits >> pos; } 

 //private:
  __mX mX;
  std::bitset<8 * sizeof(__mX)> bits;

  inline typename std::bitset<sizeof(__mX) * 8>::reference operator[](const size_t pos) { return bits[pos]; }

    private:
  blocks(const std::bitset<8 * sizeof(__mX)> & inp) : bits(inp)  {}
};

template<typename __mX>
inline blocks<__mX> xor_if(const blocks<__mX> & block1, blocks<__mX> & block2, bool flag);
template<>
inline blocks<__m128i> xor_if(const blocks<__m128i> & block1, blocks<__m128i> & block2, bool flag)
{
  return _mm_xor_si128(block1, _mm_and_si128(block2, if128_mask[flag ? 1 : 0]));
}
template<>
inline blocks<__m256i> xor_if(const blocks<__m256i> & block1, blocks<__m256i> & block2, bool flag)
{
  return _mm256_xor_si256(block1, _mm256_and_si256(block2, if256_mask[flag ? 1 : 0]));
}

inline bool xor_if(const bool & block1, bool & block2, bool flag)
{
  if(flag) 
    {
    return (block1 ^ block2);
    }
  else
  {
    return block1;
  }
  //return _mm256_xor_si256(block1, _mm256_and_si256(block2, if256_mask[flag ? 1 : 0]));
}

template<typename __mX>
inline uint8_t get_lsb(const blocks<__mX> & block);
template<>
inline uint8_t get_lsb(const blocks<__m128i> & block)
{
   __m128i vcmp = _mm_xor_si128(_mm_and_si128(block, lsb128_mask[1]), lsb128_mask[1]);
   return static_cast<uint8_t>(_mm_testz_si128(vcmp, vcmp));
}
template<>
inline uint8_t get_lsb(const blocks<__m256i> & block)
{
   __m256i vcmp = _mm256_xor_si256(_mm256_and_si256(block, lsb256_mask[1]), lsb256_mask[1]);
   return static_cast<uint8_t>(_mm256_testz_si256(vcmp, vcmp));
}

template<typename __mX>
blocks<__mX> clear_lsb(const blocks<__mX> & block);
template<>
inline blocks<__m128i> clear_lsb(const blocks<__m128i> & block)
{
   return _mm_and_si128(block, lsb128_mask_inv);
}
template<>
inline blocks<__m256i> clear_lsb(const blocks<__m256i> & block)
{
   return _mm256_and_si256(block, lsb256_mask_inv);
}

template<typename __mX>
blocks<__mX> set_lsb(const blocks<__mX> & block, const bool bit = true);
template<>
inline blocks<__m128i> set_lsb(const blocks<__m128i> & block, const bool bit)
{
  return _mm_or_si128(clear_lsb(block), lsb128_mask[bit ? 1 : 0]);
}
template<>
inline blocks<__m256i> set_lsb(const blocks<__m256i> & block, const bool bit)
{
  return _mm256_or_si256(clear_lsb(block), lsb256_mask[bit ? 1 : 0]);;
}

template<typename __mX>
inline blocks<__mX> operator|(const blocks<__mX> & block1, const blocks<__mX> & block2);
template<>
inline blocks<__m256i> operator|(const blocks<__m256i> & block1, const blocks<__m256i> & block2)
{
  return _mm256_or_si256(block1, block2);
}
template<>
inline blocks<__m128i> operator|(const blocks<__m128i> & block1, const blocks<__m128i> & block2)
{
  return _mm_or_si128(block1, block2);
}



template<typename __mX>
inline blocks<__mX> operator&(const blocks<__mX> & block1, const blocks<__mX> & block2);
template<>
inline blocks<__m256i> operator&(const blocks<__m256i> & block1, const blocks<__m256i> & block2)
{
  return _mm256_and_si256(block1, block2);
}
template<>
inline blocks<__m128i> operator&(const blocks<__m128i> & block1, const blocks<__m128i> & block2)
{
  return _mm_and_si128(block1, block2);
}



template<typename __mX>
inline blocks<__mX> operator^(const blocks<__mX> & block1, const blocks<__mX> & block2);
template<>
inline blocks<__m256i> operator^(const blocks<__m256i> & block1, const blocks<__m256i> & block2)
{
  return _mm256_xor_si256(block1, block2);
}

template<>
inline blocks<__m128i> operator^(const blocks<__m128i> & block1, const blocks<__m128i> & block2)
{
  return _mm_xor_si128(block1, block2);
}



template<typename __mX>
inline blocks<__mX> & operator^=(blocks<__mX> & block1, const blocks<__mX> & block2);
template<>
inline blocks<__m256i> & operator^=(blocks<__m256i> & block1, const blocks<__m256i> & block2)
{
  block1 = _mm256_xor_si256(block1, block2);
  return block1;
}
template<>
inline blocks<__m128i> & operator^=(blocks<__m128i> & block1, const blocks<__m128i> & block2)
{
  block1 = _mm_xor_si128(block1, block2);
  return block1;
}



template<typename __mX>
inline blocks<__mX> operator<<(const blocks<__mX> & block, const long & shift);
template<>
inline blocks<__m256i> operator<<(const blocks<__m256i> & block , const long & shift)
{
  return _mm256_or_si256(_mm256_slli_epi64(block, shift), _mm256_blend_epi32(_mm256_setzero_si256(), _mm256_permute4x64_epi64(_mm256_srli_epi64(block, 64 - shift), _MM_SHUFFLE(2,1,0,0)), _MM_SHUFFLE(3,3,3,0)));
}
template<>
inline blocks<__m128i> operator<<(const blocks<__m128i> & block, const long & shift)
{
  return _mm_or_si128(_mm_slli_epi64(block, shift), _mm_srli_epi64(_mm_slli_si128(block, 8), 64 - shift));
}


template<typename __mX>
inline blocks<__mX> operator>>(const blocks<__mX> & block, const long & shift);
template<>
inline blocks<__m256i> operator>>(const blocks<__m256i> & block, const long & shift)
{
  return _mm256_or_si256(_mm256_srli_epi64(block, shift), _mm256_blend_epi32(_mm256_setzero_si256(), _mm256_permute4x64_epi64(_mm256_slli_epi64(block, 64 - shift), _MM_SHUFFLE(0,3,2,1)), _MM_SHUFFLE(0,3,3,3)));
}
template<>
inline blocks<__m128i> operator>>(const blocks<__m128i> & block, const long & shift)
{
  return _mm_or_si128(_mm_srli_epi64(block, shift), _mm_slli_epi64(_mm_srli_si128(block, 8), 64 - shift));
}


template<typename __mX>
inline blocks<__mX> & operator>>=(blocks<__mX> & block, const long & shift)
{
  block = block >> shift;
  return block;
}


template<typename __mX>
inline blocks<__mX> & operator<<=(blocks<__mX> & block, const long & shift)
{
  block = block << shift;
  return block;
}

template<typename __mX>
inline std::ostream & operator<<(std::ostream & out, const blocks<__mX> & block)
{
    return out << block.bits;
}
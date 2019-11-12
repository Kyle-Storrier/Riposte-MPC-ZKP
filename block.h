#include "key.h"


template <typename __mX>
union blocks
{
    __mX a;
    std::bitset<8 * sizeof(__mX)> b;
    blocks(unsigned long long input = 0ULL) : b(input) { }
    //inline operator __m256i() const { return a; }
    inline blocks<__mX> & operator=(__mX val) { a = val; return *this; }
 
 public:
    inline operator __mX() const { return a; }
    std::bitset<blocksize>::reference operator[](size_t pos) { return b[pos]; }
    const bool operator[] (const size_t pos) const { return b[pos]; }
    long to_ulong() const { return b.to_ulong(); }
    size_t count() const { return b.count(); }
    size_t size() const { return b.size(); }
    void clear_lsb() { b.reset(0); }
    void set_lsb() { b.set(0); }
    uint8_t get_lsb(){return b[0]; }
    blocks set_lsb(uint8_t bit){b[0] = bit; return *this; }
    void set() {b.set();}
};

// template<typename __mX = __m256i>
// union block
// {
//     __mX a;
//     std::bitset<blocksize> b;
//     // block() { }
//     // block(long input) 
//     // {
//     //  b = input;
//     // }

//    block(unsigned long long input = 0ULL) : b(input) { } 

//  public:

//     std::bitset<blocksize>::reference operator[](size_t pos) { return b[pos]; }
//     const std::bitset<blocksize> operator[] (const size_t idx) const
//     {
//     return b[idx];
//     }

//     long to_ulong()
//     {
//         return b.to_ulong();
//     }

//     size_t count()
//     {
//         return b.count();
//     }

//     size_t size()
//     {
//         return b.size();
//     }

// };

template<typename __mX>
inline blocks<__mX> operator & (const blocks<__mX>& a1, const blocks<__mX> & a2)
{
  blocks<__mX> result;
  result.b = a1.b & a2.b;

  return result;
}

template<typename __mX>
inline blocks<__mX> operator^(const blocks<__mX>& a1, const blocks<__mX> & a2)
{
  blocks<__mX> result;
  result.b = a1.b ^ a2.b;
  return result;
}

template<typename __mX>
inline blocks<__mX> & operator^=(blocks<__mX> & result, const blocks<__mX>& a2)
{
  result.b ^= a2.b;
  return result;
}
template<typename __mX>
inline blocks<__mX> & operator^=(blocks<__mX> & result, const long & a2)
{
  result.b ^= a2;
  return result;
}
template<typename __mX>
inline blocks<__mX> operator<<(const blocks<__mX>& a1, const long & a2)
{
  blocks<__mX> result;
  result.b = a1.b << a2;
  return result;
}

template<typename __mX>
inline blocks<__mX> operator>>(const blocks<__mX>& a1, const long & a2)
{
  blocks<__mX> result;
  result.b = a1.b >> a2;
  return result;
}


template<typename __mX>
inline blocks<__mX> & operator>>=(blocks<__mX> & result, const long & a2)
{
  result.b >>= a2;
  return result;
}


template<typename __mX>
inline blocks<__mX> & operator<<=(blocks<__mX> & result, const long & a2)
{
  result.b <<= a2;
  return result;
}

template<typename __mX>
inline std::ostream &  operator<<(std::ostream &out, blocks<__mX> & result)
{
    out << result.b;

    return out;
}
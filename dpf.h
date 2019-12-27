#ifndef DPF_H__
#define DPF_H__

#include <cstdint>
#include <cmath>
//#include <x86intrin.h>
#include <bsd/stdlib.h>
#include<memory.h>

#include "LowMC.h"

#include "aes.h"

#include "prg.h"

#define L 0
#define R 1

struct cwbits { private: uint8_t t[2]; public: uint8_t & operator[](bool b) { return t[b ? 1 : 0]; } };

template<typename __mX, size_t nitems>
struct dpf_key
{
  static_assert(nitems % 64 == 0, "nitems must be a multiple of 64");
  static constexpr size_t depth = static_cast<size_t>(std::ceil(std::log2(nitems)));

  blocks<__mX> root;
  blocks<__mX> cw[depth+1];
  cwbits t[depth+1];
  blocks<__mX> leaf;
};


template<typename KEY_TYPE, typename __mX>
inline void expand(KEY_TYPE & key, const blocks<__mX> & seed, blocks<__mX> s[2], uint8_t t[2])
{
  //std::cout << "expand: seed: " << seed << std::endl;
  blocks<__mX> seedL = clear_lsb(seed); // _mm_clearlsb_si128(seed);
//  std::cout << "expand: seedL: " << (seedL ^ seed) << std::endl;
  blocks<__mX> seedR = set_lsb(seed); //  _mm_setlsb_si128(seed);
 // std::cout << "expand: seedR: " << (seedR ^ seed) << std::endl;
  s[L] = seedL;
  s[R] = seedR;
 
  #ifdef AES
  __mX ss[2] = { s[0], s[1] };
  AES_ecb_encrypt_blks(ss, 2, &key);
  s[0] = ss[0];
  s[1] = ss[1];
  #endif

  #ifdef LOWMC
    const block ss0 = s[0]; const block ss1 = s[1];
    s[0] = key.encrypt(ss0);
    s[1] = key.encrypt(ss1);
  #endif

  s[L] ^= seedL;
  t[L] = get_lsb(s[L]);
  // std::cout << "s_L = " << s[L] << std::endl;
  // std::cout << "t[L] = " << (unsigned) t[L] << std::endl;
  s[L] = clear_lsb(s[L]);

  s[R] ^= seedR;
  t[R] = get_lsb(s[R]);
  s[R] = clear_lsb(s[R]);
  
}

template<typename KEY_TYPE, typename __mX>
inline void expand_nonmpc(KEY_TYPE & key, const blocks<__mX> & seed, blocks<__mX> s[2], bool t[2])
{
  
  blocks<__mX> seedL = clear_lsb(seed); // _mm_clearlsb_si128(seed);
  blocks<__mX> seedR = set_lsb(seed); //  _mm_setlsb_si128(seed);

  s[L] = seedL;
  s[R] = seedR;
 
  // #ifdef AES
  // __mX ss[2] = { s[0], s[1] };
  // AES_ecb_encrypt_blks(ss, 2, &key);
  // s[0] = ss[0];
  // s[1] = ss[1];
  // #endif

   #ifdef LOWMC
     const block ss0 = s[0]; const block ss1 = s[1];
     s[0] = key.encrypt(ss0);
     s[1] = key.encrypt(ss1);
   #endif

   s[L] ^= seedL;
   t[L] = get_lsb(s[L]);
   s[L] = clear_lsb(s[L]);

   s[R] ^= seedR;
   t[R] = get_lsb(s[R]);
   s[R] = clear_lsb(s[R]);
  
}




 // seed0 and seed1 are the shares of the seed being expanded. s0 and s1 are the a


template<typename KEY_TYPE, typename __mX>
inline void nextlayer(KEY_TYPE & prgkey, blocks<__mX> s[2], uint8_t t[2], uint8_t bit,
  blocks<__mX> & cw, cwbits & cwt)
{
   
  blocks<__mX> s0[2], s1[2];
  uint8_t t0[2], t1[2];
  
  expand(prgkey, s[0], s0, t0);
  expand(prgkey, s[1], s1, t1);

  const uint8_t keep = (bit == 0) ? L : R, lose = 1 - keep;
  cw = clear_lsb(s0[lose] ^ s1[lose]);
  cwt[L] = t0[L] ^ t1[L] ^ bit ^ 1;
  cwt[R] = t0[R] ^ t1[R] ^ bit;

  s[L] = xor_if(s0[keep], cw, !t[L]);
  t[L] = t0[keep] ^ (t[L] & cwt[keep]);
  
  s[R] = xor_if(s1[keep], cw, !t[R]);
  t[R] = t1[keep] ^ (t[R] & cwt[keep]);
}

template<typename KEY_TYPE, typename __mX>
inline void leafayer(KEY_TYPE & prgkey, blocks<__mX> s[2], uint8_t t[2], uint8_t bit,
  blocks<__mX> & cw, cwbits & cwt, blocks<__mX> out)
{
   

  cw = s[0] ^ s[1] ^ out;

}

template <typename KEY_TYPE, size_t nitems, typename __mX>
void gen(KEY_TYPE & prgkey, size_t point, dpf_key<__mX, nitems> dpfkey[2])
{

  blocks<__mX> s[2], s0[2], s1[2];
  uint8_t t[2], t0[2], t1[2];

   arc4random_buf(s, 2 * sizeof(blocks<__mX>));
   // srand(1);
   // s[0] = rand(); s[1] = rand();
   t[0] = get_lsb(s[0]);

   dpfkey[0].root = s[0];
   t[1] = !t[0];

   dpfkey[1].root = set_lsb(s[1], t[1]);

   for (size_t i = 0; i <= dpf_key<__mX, nitems>::depth; ++i)
   {
     const uint8_t bit = (point >> (dpf_key<__mX, nitems>::depth - i - 1)) & 1U;

     if(i < dpf_key<__mX, nitems>::depth) 
     {
      nextlayer(prgkey, s, t, bit, dpfkey[0].cw[i], dpfkey[0].t[i]);
     }   
     else
     {
      leafayer(prgkey, s, t, bit, dpfkey[0].cw[i], dpfkey[0].t[i], blocks<__mX> (_mm256_set1_epi64x(-1)) );
     }
   }
    
   __mX ss[2] = { s[0], s[1] };

   memcpy(&dpfkey[1].cw, &dpfkey[0].cw, sizeof(dpf_key<__mX, nitems>::cw));
   memcpy(&dpfkey[1].t, &dpfkey[0].t, sizeof(dpf_key<__mX, nitems>::t));
}




template <typename KEY_TYPE, size_t nitems, typename __mX>
inline void evalfull3(KEY_TYPE & prgkey, dpf_key<__mX, nitems> & dpfkey, blocks<__mX> * s, uint8_t * t)
{
   constexpr size_t depth = dpf_key<__mX, nitems>::depth;

   blocks<__mX> * s_[2] = { s, s + nitems/2 };
   uint8_t * t_[2] = { t, t + nitems/2 };

   int curlayer = depth % 2;

    s_[curlayer][0] = dpfkey.root;
    t_[curlayer][0] = get_lsb(dpfkey.root);

   blocks<__mX> child[2];
   uint8_t ts[2];
   for (size_t i = 0; i < depth; ++i)
   {
     curlayer = 1 - curlayer;
     const size_t itemnumber = std::max((nitems / (1ULL << (depth - i))), 1ULL);
     
     for (size_t j = 0; j < itemnumber; ++j)
      {
        expand(prgkey, s_[1-curlayer][j], child, ts);
         // std::cout << (unsigned) t_[1-curlayer][j] << " "; 
        s_[curlayer][2*j] = xor_if(child[L], dpfkey.cw[i], !t_[1-curlayer][j]);
        t_[curlayer][2*j] = ts[L] ^ dpfkey.t[i][L] & t_[1-curlayer][j];
        //std::cout << "[L: " << (unsigned) t_[curlayer][2*j] << " ]" ;
      
        if (2*j+1 < 2*itemnumber)
        {
          //std::cout << (unsigned) t_[1-curlayer][j] << " "; 
          s_[curlayer][2*j+1] = xor_if(child[R], dpfkey.cw[i], !t_[1-curlayer][j]);
          t_[curlayer][2*j+1] = ts[R] ^ dpfkey.t[i][R] & t_[1-curlayer][j];
      //    std::cout << "[R: " << (unsigned) t_[curlayer][2*j + 1] << " ]" ;
          
        }
     }
   }

   for(size_t i = 0; i < nitems; ++i)
   {
    if(t_[curlayer][i]) s[i] = s[i] ^ dpfkey.cw[depth];
   }
}



#endif

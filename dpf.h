#ifndef DPF_H__
#define DPF_H__

#include <omp.h>
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


// static const __m128i lsb_mask[3] = {
//   _mm_setzero_si128(),                                     // 0b00...0000
//   _mm_set_epi64x(0,1)                                      // 0b00...0001
// };
// static const __m128i lsb_mask_inv = _mm_set_epi64x(-1,-2); // 0b11...1110
// static const __m128i if_mask[2] = {
//   _mm_setzero_si128(),                                     // 0b00...0000
//   _mm_set1_epi8(-1)                                        // 0b11...1111
// };

// inline uint8_t _mm_getlsb_si128(const __m128i & x)
// {
//   __m128i vcmp = _mm_xor_si128(_mm_and_si128(x, lsb_mask[1]), lsb_mask[1]);
//   return static_cast<uint8_t>(_mm_testz_si128(vcmp, vcmp));
// }
// inline __m128i _mm_clearlsb_si128(const __m128i & x)
// {
//   return _mm_and_si128(x, lsb_mask_inv);
// }
// inline __m128i _mm_setlsb_si128(const __m128i & x, bool b = true)
// {
//   return _mm_or_si128(_mm_clearlsb_si128(x), lsb_mask[b ? 1 : 0]);
// }
// inline __m128i _mm_xorif_si128(const __m128i & x, const __m128i & y, bool b)
// {
//   return _mm_xor_si128(x, _mm_and_si128(y, if_mask[b ? 1 : 0]));
// }

struct cwbits { private: uint8_t t[2]; public: uint8_t & operator[](bool b) { return t[b ? 1 : 0]; } };

template<typename __mX, size_t nitems> struct dpf_key
{
  static_assert(nitems % 64 == 0, "nitems must be a multiple of 64");
  static constexpr size_t depth = static_cast<size_t>(std::ceil(std::log2(nitems)));

  blocks<__mX> root;
  blocks<__mX> cw[depth];
  cwbits t[depth];
  //profile<3*precision> leaf;
};


template<typename KEY_TYPE, typename __mX>
inline void expand(KEY_TYPE & key, const blocks<__mX> & seed, blocks<__mX> s[2], uint8_t t[2])
{
  
  blocks<__mX> seedL = seed; // _mm_clearlsb_si128(seed);
  seedL.clear_lsb();

  blocks<__mX> seedR = seed; //  _mm_setlsb_si128(seed);
  seedR.set_lsb();

  s[L] = seedL;
  s[R] = seedR;
 
  #ifdef AES
  __mX ss[2] = {s[0], s[1]};
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
  t[L] = s[L].get_lsb();
  s[L].clear_lsb();

  s[R] ^= seedR;
  t[R]  = s[R].get_lsb();
  s[R].clear_lsb();
  
}

template<typename KEY_TYPE, typename __mX>
inline void nextlayer(KEY_TYPE & aeskey, blocks<__mX> s[2], uint8_t t[2], uint8_t bit,
  blocks<__mX> & cw, cwbits & cwt)
{
   
  blocks<__mX> s0[2], s1[2];
  uint8_t t0[2], t1[2];
  
  expand(aeskey, s[0], s0, t0);
  expand(aeskey, s[1], s1, t1);

  const uint8_t keep = (bit == 0) ? L : R, lose = 1 - keep;
  blocks<__mX> s0xors1 = s0[lose] ^ s1[lose];
  s0xors1.clear_lsb();
  cw = s0xors1;
  cwt[L] = t0[L] ^ t1[L] ^ bit ^ 1;
  cwt[R] = t0[R] ^ t1[R] ^ bit;

  if(t[L])
  {
    blocks<__mX> y = 0;
    s[L] = s0[keep] ^ (cw & y);
  }
  else
  {  
    blocks<__mX> y;
    y.set();
    s[L] = s0[keep] ^ (cw & y);
  }

  t[L] = t0[keep] ^ (t[L] & cwt[keep]);
  
  if(t[R])
  {
    blocks<__mX> y = 0;
    s[R] = s1[keep] ^ (cw & y);
  }
  else
  {  
    blocks<__mX> y;
    y.set();
    s[R] = s1[keep] ^ (cw & y);
  }

  t[R] = t1[keep] ^ (t[R] & cwt[keep]);
}

// template<typename KEY_TYPE, typename __mX>
// inline void leaflayer(KEY_TYPE & aeskey, __mX s[2], uint8_t t[2],
//   profile<3*precision> & lastL, profile<3*precision> & lastR)
// {
//   constexpr size_t len = sizeof(profile<precision>) / sizeof(__m128i);

//   profile<3*precision> s0, s1;
//   PRG(aeskey, s[L], (__m128i*)&s0, dim128);
//   PRG(aeskey, s[R], (__m128i*)&s1, dim128);
   
//   arc4random_buf(&lastL, sizeof(profile<precision>));  
  
//   lastR = (t[L] ? (s1 - s0) - lastL :  (s0 - s1) - lastL);
  
//   lastL = -lastL; lastR = -lastR;

// }

template <typename KEY_TYPE, size_t nitems, typename __mX>
void gen(KEY_TYPE & aeskey, size_t point, dpf_key<__mX, nitems> dpfkey[2])
{

  blocks<__mX> s[2], s0[2], s1[2];
  uint8_t t[2], t0[2], t1[2];
 

   arc4random_buf(s, 2 * sizeof(blocks<__mX>));
   
   t[0] = s[0].get_lsb();

   dpfkey[0].root = s[0];
   t[1] = !t[0];

   dpfkey[1].root = s[1].set_lsb(t[1]);

   for (size_t i = 0; i < dpf_key<__mX, nitems>::depth; ++i)
   {
     const uint8_t bit = (point >> (dpf_key<__mX, nitems>::depth - i - 1)) & 1U;
     nextlayer(aeskey, s, t, bit, dpfkey[0].cw[i], dpfkey[0].t[i]);
   }
    
   __mX ss[2] = {s[0], s[1]};
  
  // leaflayer(aeskey, ss, t, dpfkey[0].leaf, dpfkey[1].leaf);

   memcpy(&dpfkey[1].cw, &dpfkey[0].cw, sizeof(dpf_key<__mX, nitems>::cw));
   memcpy(&dpfkey[1].t, &dpfkey[0].t, sizeof(dpf_key<__mX, nitems>::t));
}



 template <typename KEY_TYPE, size_t nitems, typename __mX>
 inline void evalfull3(KEY_TYPE & aeskey, dpf_key<__mX, nitems> & dpfkey, blocks<__mX> * s, uint8_t * t)
 {
   constexpr size_t depth = dpf_key<__mX, nitems>::depth;

   blocks<__mX> * s_[2] = { s, s + nitems/2 };
   uint8_t * t_[2] = { t, t + nitems/2 };

   int curlayer = depth % 2;

    s_[curlayer][0] = dpfkey.root;
    t_[curlayer][0] = dpfkey.root.get_lsb();

   blocks<__mX> child[2];
   uint8_t ts[2];
   for (size_t i = 0; i < depth; ++i)
   {
     curlayer = 1 - curlayer;
     const size_t itemnumber = std::max((nitems / (1ULL << (depth - i))), 1ULL);
     
     size_t sub_tree_size = 1ULL << depth; 
     for (size_t j = 0; j < itemnumber; ++j)
      {
       expand(aeskey, s_[1-curlayer][j], child, ts);

        if(t_[1-curlayer][j])
        {
          blocks<__mX> y = 0;
          s_[curlayer][2*j] = child[L] ^ (dpfkey.cw[i] & y);
        }
        else
        {  
          blocks<__mX> y;
          y.set();
          s_[curlayer][2*j] = child[L] ^ (dpfkey.cw[i] & y);
        }

       t_[curlayer][2*j] = ts[L] ^ dpfkey.t[i][L] & t_[1-curlayer][j];
       if (2*j+1 < 2*itemnumber)
       {
        if(t_[1-curlayer][j])
        {
          blocks<__mX> y = 0;
          s_[curlayer][2*j+1] = child[R] ^ (dpfkey.cw[i] & y);
        }
        else
        {  
          blocks<__mX> y;
          y.set();
          s_[curlayer][2*j+1] = child[R] ^ (dpfkey.cw[i] & y);
        }

         t_[curlayer][2*j+1] = ts[R] ^ dpfkey.t[i][R] & t_[1-curlayer][j];
       }
     }
   }
 }


// inline void expand(AES_KEY & aeskey, const __m128i & seed, __m128i s[2],
//   uint8_t t[2])
// {
//   const __m128i seedL = _mm_clearlsb_si128(seed);
//   const __m128i seedR = _mm_setlsb_si128(seed);

//   s[L] = seedL;
//   s[R] = seedR;

//   AES_ecb_encrypt_blks(s, 2, &aeskey);

//   s[L] = _mm_xor_si128(s[L], seedL);
//   t[L] = _mm_getlsb_si128(s[L]);
//   s[L] = _mm_clearlsb_si128(s[L]);

//   s[R] = _mm_xor_si128(s[R], seedR);
//   t[R] = _mm_getlsb_si128(s[R]);
//   s[R] = _mm_clearlsb_si128(s[R]);
// }




// template <size_t nitems>
// inline void dolayer(AES_KEY & aeskey,__m128i * s, uint8_t * t,
//   const size_t stepsize, dpf_key<nitems> dpfkey, size_t i, __m128i child[2],
//   uint8_t ts[2])
// {
//   for (size_t j = 0; j < nitems; j += 2*stepsize)
//   {
//     expand(aeskey, s[j], child, ts);
//     if (j+stepsize < nitems)
//     {
//       s[j+stepsize] = _mm_xorif_si128(child[R], dpfkey.cw[i], t[j]);
//       t[j+stepsize] = ts[R] ^ dpfkey.t[i][R] & t[j];
//     }
//     s[j] = _mm_xorif_si128(child[L], dpfkey.cw[i], t[j]);
//     t[j] = ts[L] ^ dpfkey.t[i][L] & t[j];
//   }
// }












// template <size_t nitems>
// inline void evalfull(AES_KEY & aeskey, dpf_key<nitems> & dpfkey, __m128i * s,
//   uint8_t * t)
// {
//   s[0] = dpfkey.root;
//   t[0] = _mm_getlsb_si128(dpfkey.root);

//   __m128i child[2];
//   uint8_t ts[2];
//   size_t stepsize = 1ULL << dpf_key<nitems>::depth - 1;
//   for (size_t i = 0; i < dpf_key<nitems>::depth; ++i, stepsize /= 2)
//   {
//     dolayer(aeskey, s, t, stepsize, dpfkey, i, child, ts);
//   }
// }

// template <size_t nitems>
// inline void evalfull2(AES_KEY & aeskey, dpf_key<nitems> & dpfkey, __m128i ** s,
//   uint8_t ** t)
// {
//   s[0][0] = dpfkey.root;
//   t[0][0] = _mm_getlsb_si128(dpfkey.root);

//   int curlayer = 1;

//   for (size_t i = 0; i < dpf_key<nitems>::depth; ++i)
//   {
//     const size_t itemnumber = std::max(nitems / (1ULL << (dpf_key<nitems>::depth - i)), 1ULL);
//     for (size_t j = 0; j < itemnumber; ++j)
//     {
//       __m128i * child = &s[curlayer][2*j];
//       uint8_t * ts = &t[curlayer][2*j];

//       blocks<__m128i> s__ = s[1-curlayer][j];
//       blocks<__m128i> *child_ = &s[curlayer][2*j];

//       expand<AES_KEY, __m128i>(aeskey, s__, child_ , ts);
      
//       expand(aeskey, s[1-curlayer][j], child, ts);


//       child[L] = _mm_xorif_si128(child[L], dpfkey.cw[i], t[1-curlayer][j]);
//       ts[L] ^= dpfkey.t[i][L] & t[1-curlayer][j];
//       child[R] = _mm_xorif_si128(child[R], dpfkey.cw[i], t[1-curlayer][j]);
//       ts[R] ^= dpfkey.t[i][R] & t[1-curlayer][j];
//     }
//     curlayer = 1 - curlayer;
//   }
// }





// template <size_t nitems>
// inline void evalfull_interval ( 
//                                 AES_KEY & aeskey, dpf_key<nitems> & dpfkey,
//                                __m128i * s, uint8_t * t, size_t from, size_t to
//                               )
// {
//   size_t interval = 0;

//   size_t interval_len = to - from;
//   constexpr size_t depth = dpf_key<nitems>::depth;

//   __m128i * s_[2] = { s, s + interval_len};
//   uint8_t * t_[2] = { t, t + interval_len};

//   int curlayer = depth % 2;

//   s_[curlayer][0] = dpfkey.root;
//   t_[curlayer][0] = _mm_getlsb_si128(dpfkey.root);
 
  
//   __m128i child[2];
  
//   uint8_t ts[2];

//   size_t expanded = 0;
//   size_t skipped_prev_layer = 0;
//   size_t book_keep[2][2* nitems] = {0};
//   for (size_t i = 0; i < depth  ; ++i)
//   {
   
//    curlayer = 1 - curlayer;
//    const size_t itemnumber = std::max((nitems / (1ULL << (depth - i))), 1ULL);

//    size_t sub_tree_size = 1ULL << (depth - i);

//    interval = 0;
//    size_t skipped = 0;
//    expanded = 0;
//    for (size_t j = 0; j < itemnumber; ++j)
//     {
     
//      size_t sub_tree_from = j * sub_tree_size;
//      size_t sub_tree_to   = (j + 1) * sub_tree_size -  1;
//      if(i != 0 && (to >= sub_tree_to && from >= sub_tree_to || to <= sub_tree_from && from <= sub_tree_from))
//       {
//         ++skipped;
//         std::cout << "skipped the sub tree at depth " << i << " and breadth " << j << std::endl;
//         std::cout << "subtree skipped: " << sub_tree_from << " to " << sub_tree_to << std::endl << std::endl;
//       }     
//      else
//       {  
        
//         expand(aeskey, s_[1-curlayer][book_keep[1-curlayer][j]], child, ts);
        
//         book_keep[curlayer][2 * j] = 2 * interval;

//         s_[curlayer][2*interval] = _mm_xorif_si128(child[L], dpfkey.cw[i], t_[1-curlayer][book_keep[1-curlayer][j]]);         
//         t_[curlayer][2*interval] = ts[L] ^ dpfkey.t[i][L] & t_[1-curlayer][book_keep[1-curlayer][j]];

//         std::cout << "t_[" << curlayer << "][" << 2 * interval << "] = " << (double)t_[curlayer][2*interval] << std::endl;
//         std::cout << "s_[" << curlayer << "][" << 2 * interval << "] = " << (double)s_[curlayer][2*interval][0] << " " 
//                   << (double)s_[curlayer][2*interval][1] << std::endl;

//         if (2 * j + 1 < 2 * itemnumber)
//         {
          
//           book_keep[curlayer][2 * j + 1] = 2 * interval + 1;  

//           s_[curlayer][2*interval+1] = _mm_xorif_si128(child[R], dpfkey.cw[i], t_[1-curlayer][book_keep[1-curlayer][j]]);              
//           t_[curlayer][(2*interval) +1 ] = ts[R] ^ dpfkey.t[i][R] & t_[1-curlayer][book_keep[1-curlayer][j]];

//           std::cout << "t_[" << curlayer << "][" << 2 * interval + 1 << "] = " << (double)t_[curlayer][2*interval + 1] << std::endl;
//           std::cout << "s_[" << curlayer << "][" << 2 * interval + 1 << "] = " << (double)s_[curlayer][2*interval + 1][0] << " " 
//                   << (double)s_[curlayer][2*interval + 1][1] << std::endl;
//         }

//         ++interval;
//         ++expanded;
//       }

//     }
//      skipped_prev_layer = skipped;
//   }
// }


























// template <size_t nitems, bool negate = false>
// inline void dpf_mux(AES_KEY & aeskey,
//   dpf_key<nitems> & dpfkey, profile<precision> * iprofiles,
//   profile<precision> & itm, __m128i * s, uint8_t * t, fixed_t<2*precision>& cnt2)
// {
//   //evalfull(aeskey, dpfkey, s[0], t[0]); uint8_t * tt = t[0];
//   // std::cout << "mux called" << std::endl;
//   evalfull3(aeskey, dpfkey, s, t); 
//   uint8_t * tt = t;// t[dpf_key<nitems>::depth % 2];
//   memset(&itm, 0, sizeof(profile<precision>));
//   uint64_t cnt = 0;
  
//   for (size_t i = 0; i < nitems; ++i)
//   {
//     if (tt[i])
//     {
//       if constexpr(negate)
//       {
//         itm -= iprofiles[i];
//       }
//       else
//       {
//         itm += iprofiles[i];
//       }
//       cnt++;
//     }
//   }
//   cnt2 = negate ? -static_cast<double>(cnt) : static_cast<double>(cnt);
//     //return cnt;
// }


// template <size_t nitems, bool negate = false>
// inline void dpf_mux2(AES_KEY & aeskey,
//   dpf_key<nitems> & dpfkey, profile<precision> * iprofiles,
//   profile<precision> & itm, __m128i * s, uint8_t * t, fixed_t<2*precision>& cnt2, size_t j)
// {
  
//   //evalfull3(aeskey, dpfkey, s, t); 
  
//   uint8_t * tt = t;// t[dpf_key<nitems>::depth % 2];
  
//   memset(&itm, 0, sizeof(profile<precision>));
  
//   uint64_t cnt = 0;
  
//   for (size_t i = 0; i < nitems; ++i)
//   {  
//    if (tt[i])
//     {
//       if constexpr(negate)
//       {
//         itm -= iprofiles[i];
//   //if(j == 2)      std::cout << "added: " << i << "--->" << iprofiles[i] << std::endl;

//       }

//       else
//       {
//         itm += iprofiles[i];
//    // if(j ==2)        std::cout << "added: " << i << "--->" << iprofiles[i] << std::endl;

//       }

//    //  if(j == 2) std::cout << "itm: " << i << " - >" << itm << std::endl;
//       cnt++;
//     }
//   }
//    cnt2 = negate ? -static_cast<double>(cnt) : static_cast<double>(cnt);
// }



// template <size_t nitems, bool negate = false>
// inline void dpf_demux(AES_KEY & aeskey, dpf_key<nitems> & dpfkey,
//   const profile<3*precision> & finalblock, profile<3*precision> * iprofiles,
//   const __m128i * s, const uint8_t * t)
// {

//      profile<3*precision> tmp;
//      //std::cout << "dpf_demux: " << std::endl;
//      for (size_t i = 0; i < nitems; ++i)
//      {
//        PRG(aeskey, s[i], (__m128i *)&tmp, dim128);
       
//        if (t[i])
//        {
//         tmp -= finalblock;
//        }
//        else{

//         //std::cout << i << ": " << tmp << std::endl;
//        }

//        if constexpr(negate) { iprofiles[i] = iprofiles[i] + tmp; /*std::cout << "N -> " << iprofiles[i] << std::endl;*/ }
//        else { iprofiles[i] = iprofiles[i] - tmp; /*std::cout << "P-> " << iprofiles[i] << std::endl;*/  }      
//      }
// }

#endif

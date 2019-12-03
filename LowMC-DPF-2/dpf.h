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
  blocks<__mX> cw[depth];
  cwbits t[depth];
  //profile<3*precision> leaf;
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


template<typename KEY_TYPE, typename __mX>
inline void expand_mpc_verify(KEY_TYPE & key, const blocks<__mX> seed, const std::vector<block> c2, block blind[rounds], block gamma[rounds], bool party)
{

 // key.encrypt_MPC_verify(seed, c2, blind, gamma, party);
}


 // seed0 and seed1 are the shares of the seed being expanded. s0 and s1 are the a
 template<typename KEY_TYPE, typename __mX>
 inline std::pair<std::vector<block>,std::vector<block>> expand_mpc(KEY_TYPE & key, const blocks<__mX> & seed0, const blocks<__mX> & seed1,  blocks<__mX> s0[2], bool t0[2], 
                        blocks<__mX> s1[2], bool t1[2], const block blind0[rounds], const block blind1[rounds], block gamma[2][rounds])
 {
  
    std::vector<block> c0(2 + rounds + rounds);
    std::vector<block> c1(2 + rounds + rounds);
  

   const blocks<__mX> seed0L = clear_lsb(seed0); // _mm_clearlsb_si128(seed);
   const blocks<__mX> seed0R = set_lsb(seed0); //  _mm_setlsb_si128(seed);

   const blocks<__mX> seed1L = clear_lsb(seed1); // _mm_clearlsb_si128(seed);
   const blocks<__mX> seed1R = clear_lsb(seed1); //  _mm_setlsb_si128(seed);

   s0[L] = seed0L;
   s0[R] = seed0R;
 
   s1[L] = seed1L;
   s1[R] = seed1R;

   
    std::pair<std::vector<block>,std::vector<block>> encrypt_outL = key.encrypt_MPC_proof(seed0L, seed1L, blind0, blind1, gamma);

    for(size_t j = 0; j < rounds; ++j)
    {
      c0[2 + 2*j]     = encrypt_outL.first[j+1];
      c0[2 + 2*j + 1] = encrypt_outL.second[j+1];
    }


    std::pair<std::vector<block>,std::vector<block>> encrypt_outR = key.encrypt_MPC_proof(seed0R, seed1R, blind0, blind1, gamma);

    for(size_t j = 0; j < rounds; ++j)
    {
      c1[2 + 2*j]     = encrypt_outR.first[j+1];
      c1[2 + 2*j + 1] = encrypt_outR.second[j+1];
    }

    s0[L] = encrypt_outL.first[rounds];
    s0[R] = encrypt_outR.first[rounds];
 
    s1[L] = encrypt_outL.second[rounds];
    s1[R] = encrypt_outR.second[rounds];

    s0[L] ^= seed0L; 
    s1[L] ^= seed1L;
    
    t0[L] = get_lsb(s0[L]); 
    t1[L] = get_lsb(s1[L]);
     
    s0[L] = clear_lsb(s0[L]); 
    s1[L] = clear_lsb(s1[L]);

    s0[R] ^= seed0R; s1[R] ^= seed1R;
    t0[R] = get_lsb(s0[R]); t1[R] = get_lsb(s1[R]);
    s0[R] = clear_lsb(s0[R]); s1[R] = clear_lsb(s1[R]);   

    c0[L] = s0[L];
    c1[L] = s1[L];

    c0[R] = s0[R];
    c1[R] = s1[R];

    return std::make_pair(c0, c1);
  
 }


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

   for (size_t i = 0; i < dpf_key<__mX, nitems>::depth; ++i)
   {
     const uint8_t bit = (point >> (dpf_key<__mX, nitems>::depth - i - 1)) & 1U;
     nextlayer(prgkey, s, t, bit, dpfkey[0].cw[i], dpfkey[0].t[i]);
   }
    
   __mX ss[2] = { s[0], s[1] };
  
  // leaflayer(prgkey, ss, t, dpfkey[0].leaf, dpfkey[1].leaf);

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
    // std::cout << std::endl << "------------------" << std::endl;
   }
}


// inline void expand(AES_KEY & prgkey, const __m128i & seed, __m128i s[2],
//   uint8_t t[2])
// {
//   const __m128i seedL = _mm_clearlsb_si128(seed);
//   const __m128i seedR = _mm_setlsb_si128(seed);

//   s[L] = seedL;
//   s[R] = seedR;

//   AES_ecb_encrypt_blks(s, 2, &prgkey);

//   s[L] = _mm_xor_si128(s[L], seedL);
//   t[L] = _mm_getlsb_si128(s[L]);
//   s[L] = _mm_clearlsb_si128(s[L]);

//   s[R] = _mm_xor_si128(s[R], seedR);
//   t[R] = _mm_getlsb_si128(s[R]);
//   s[R] = _mm_clearlsb_si128(s[R]);
// }




// template <size_t nitems>
// inline void dolayer(AES_KEY & prgkey,__m128i * s, uint8_t * t,
//   const size_t stepsize, dpf_key<nitems> dpfkey, size_t i, __m128i child[2],
//   uint8_t ts[2])
// {
//   for (size_t j = 0; j < nitems; j += 2*stepsize)
//   {
//     expand(prgkey, s[j], child, ts);
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
// inline void evalfull(AES_KEY & prgkey, dpf_key<nitems> & dpfkey, __m128i * s,
//   uint8_t * t)
// {
//   s[0] = dpfkey.root;
//   t[0] = _mm_getlsb_si128(dpfkey.root);

//   __m128i child[2];
//   uint8_t ts[2];
//   size_t stepsize = 1ULL << dpf_key<nitems>::depth - 1;
//   for (size_t i = 0; i < dpf_key<nitems>::depth; ++i, stepsize /= 2)
//   {
//     dolayer(prgkey, s, t, stepsize, dpfkey, i, child, ts);
//   }
// }

// template <size_t nitems>
// inline void evalfull2(AES_KEY & prgkey, dpf_key<nitems> & dpfkey, __m128i ** s,
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

//       expand<AES_KEY, __m128i>(prgkey, s__, child_ , ts);
      
//       expand(prgkey, s[1-curlayer][j], child, ts);


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
//                                 AES_KEY & prgkey, dpf_key<nitems> & dpfkey,
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
        
//         expand(prgkey, s_[1-curlayer][book_keep[1-curlayer][j]], child, ts);
        
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
// inline void dpf_mux(AES_KEY & prgkey,
//   dpf_key<nitems> & dpfkey, profile<precision> * iprofiles,
//   profile<precision> & itm, __m128i * s, uint8_t * t, fixed_t<2*precision>& cnt2)
// {
//   //evalfull(prgkey, dpfkey, s[0], t[0]); uint8_t * tt = t[0];
//   // std::cout << "mux called" << std::endl;
//   evalfull3(prgkey, dpfkey, s, t); 
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
// inline void dpf_mux2(AES_KEY & prgkey,
//   dpf_key<nitems> & dpfkey, profile<precision> * iprofiles,
//   profile<precision> & itm, __m128i * s, uint8_t * t, fixed_t<2*precision>& cnt2, size_t j)
// {
  
//   //evalfull3(prgkey, dpfkey, s, t); 
  
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
// inline void dpf_demux(AES_KEY & prgkey, dpf_key<nitems> & dpfkey,
//   const profile<3*precision> & finalblock, profile<3*precision> * iprofiles,
//   const __m128i * s, const uint8_t * t)
// {

//      profile<3*precision> tmp;
//      //std::cout << "dpf_demux: " << std::endl;
//      for (size_t i = 0; i < nitems; ++i)
//      {
//        PRG(prgkey, s[i], (__m128i *)&tmp, dim128);
       
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

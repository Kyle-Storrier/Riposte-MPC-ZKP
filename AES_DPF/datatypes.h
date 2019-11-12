#ifndef DATATYPES_H__
#define DATATYPES_H__

#include <cstring>
#include <x86intrin.h>
#include <cmath>

#include "fixed.h"
#include "common.h"


#define PR 13

static constexpr size_t dim256 = DIM * sizeof(uint64_t) / sizeof(__m256i);
static constexpr size_t dim128 = DIM * sizeof(uint64_t) / sizeof(__m128i);
static constexpr size_t precision = 64 / 3 - std::log2(DIM) - 0.5;

template <size_t precision>
union profile
{
  __m256i m256[dim256];
  __m128i m128[dim128];
  uint64_t u64[DIM];
  fixed_t<precision> f64[DIM];
  profile() { }
  profile(double d[DIM]) { for (size_t i = 0; i < DIM; ++i) f64[i] = d[i]; }
  inline profile<precision> & operator=(double d[DIM]) { for (size_t i = 0; i < DIM; ++i) f64[i] = d[i]; return *this; }
  inline profile<precision> operator-()
  {
    profile<precision> p;
    for (size_t i = 0; i < DIM; ++i) p.u64[i] = -u64[i];
    return p;
  }
  template <size_t precision2>
  inline operator profile<precision2>() 
  {
    profile<precision2> p;

    double scale1 = 1.0L * static_cast<double>(1ULL << precision);
    double scale2 = 1.0L * static_cast<double>(1ULL << precision2);
    for (size_t i = 0; i < DIM; ++i)
    {
      if(uint64_t(scale2/scale1) == 0){
      p.f64[i] = fixed_t<precision2>(f64[i]);
      }
      else{
      p.f64[i].val = (f64[i].val) * (uint64_t(scale2/scale1));  
      }
    } 
    return p;
  }
};

// union u64_array
// {
//   uint64_t u64[DIM];
//   u64_array(){ }
//   inline u64_array & operator=(uint64_t d[DIM]) { for (size_t i = 0; i < DIM; ++i) u64[i] = d[i]; return *this; }
// };

struct p3_beavers
{
  profile<precision>   A_m; // Used for MUX
  fixed_t<3*precision> RRtCW; // Used for MUX
  fixed_t<3*precision> sanity_blinds;
 
  static constexpr size_t size = 1 * sizeof(profile<precision>) + 1 * sizeof(fixed_t<3*precision>)
  								 + sizeof(fixed_t<3*precision>);// + sizeof(bool) + 2 * sizeof(uint64_t); 
};


struct beaver3
{

  profile<precision> u0_blind, u1_blind;  
  
  fixed_t<(precision - PR)> aa0, aa1, bb0, bb1;
  
  fixed_t<precision/2> r20_blindred, r21_blindred;

  profile<precision/2> u0_blindred, u1_blindred; 

  fixed_t<precision> r20_blind, r21_blind;

  fixed_t<2*precision> alpha; // used to blind \langle \vec{u}_0, \vec{u}_1 \rangle
  profile<precision>   beta;
  fixed_t<2*(precision-PR)> gamma;
  profile<precision> sanity_blind;
  
  uint64_t b0_blind_0, b1_blind_0, b0_blind_1, b1_blind_1;

  uint64_t b0_blind_0_u[DIM], b1_blind_0_u[DIM], b0_blind_1_u[DIM], b1_blind_1_u[DIM];  //used for adjusting uprofiles

  uint64_t b0_blind_0_ur[DIM], b1_blind_0_ur[DIM], b0_blind_1_ur[DIM], b1_blind_1_ur[DIM]; //used for adjusting ur

  static constexpr size_t size =  2 *  sizeof(profile<precision>) 

                                + 4 * sizeof(fixed_t<precision-PR>)   

                                + 2 *  sizeof(fixed_t<precision/2>)

                                + 2 *  sizeof(profile<precision/2>)

                                + 2 *  sizeof(fixed_t<precision>)

                                + 1 * sizeof(fixed_t<2*precision>) // used to blind \langle \vec{u}_0, \vec{u}_1 \rangle
                                + 1 * sizeof(profile<precision>)
                                + 1 * sizeof(fixed_t<2*(precision-PR)>)
                                + 1 * sizeof(profile<precision>)

                                + 4 * sizeof(uint64_t) 

                                + 4 * DIM * sizeof(uint64_t)

                                + 4 * DIM * sizeof(uint64_t);

};

struct norm_stuff
{
  fixed_t<precision-PR> rr; 
  fixed_t<precision> r;

  profile<precision> u_blind; 
  fixed_t<2*precision> s;

  fixed_t<precision-PR>    rr_blind, dp_blind, rr_blinded; 
  fixed_t<2*(precision-PR)> t;
  
  profile<precision/2> u_blind_reduced;
  fixed_t<precision/2> r_blind;
  profile<precision> tt;

  profile<precision> u_blindA;
  fixed_t<precision> r_blindA, r_blindedA;
  profile<2 * precision> ttA;


  uint64_t   b0_blind_u0[DIM], b1_blind_u0[DIM], bb_u[DIM];

  uint64_t   b0_blind_u1[DIM], b1_blind_u1[DIM], bb_u1[DIM];

  uint64_t   b0_blind_ur0[DIM], b1_blind_ur0[DIM], bb_ur[DIM];

  uint64_t   b0_blind_ur1[DIM], b1_blind_ur1[DIM];


  static constexpr size_t size =  1 * sizeof(fixed_t<precision-PR>) + 1 * sizeof(fixed_t<precision>) 

                                + 1 * sizeof(profile<precision>) + 1 * sizeof(fixed_t<2*precision>) 

                                + 3 * sizeof(fixed_t<precision-PR>) + 1 * sizeof(fixed_t<2*(precision-PR)>)

                                + 1 * sizeof(profile<precision/2>) + 1 * sizeof(fixed_t<precision/2>) + 1 * sizeof(profile<precision>)

                                + 1 * sizeof(profile<precision>) + 2 * sizeof(fixed_t<precision>) + 1 * sizeof(profile<2*precision>)

                               
                                + 3 * DIM * sizeof(uint64_t)
                                + 3 * DIM * sizeof(uint64_t)

                                + 3 * DIM * sizeof(uint64_t)
                                + 2 * DIM * sizeof(uint64_t);

                               
                              
};


struct beaver1
{
  //norm_stuff itm_norm;
  //norm_stuff usr_norm;
  profile<precision> A_m;
  profile<precision> usr;
  profile<precision> itm;
  fixed_t<2*precision> dp;
  fixed_t<3*precision> RRtCW;

  fixed_t<3*precision> sanity_blind;

  static constexpr size_t size =// 2 * sizeof(norm_stuff) + 
                                3 * sizeof(profile<precision>) + 1 * sizeof(fixed_t<2*precision>) 
                                + 1 * sizeof(fixed_t<3*precision>) + 1 * sizeof(fixed_t<3*precision>)  ;
};

struct beaver2
{
  // norm_stuff itm_norm;
  // norm_stuff usr_norm;
  profile<precision> A_m;
  profile<3*precision> usr;
  profile<3*precision> itm;
  fixed_t<2*precision> dp;
  fixed_t<3*precision> RRtCW;

  fixed_t<3*precision> sanity_blind;

  static constexpr size_t size =  //2 * sizeof(norm_stuff) +
                                 1 * sizeof(profile<precision>) 
                                + 2 * sizeof(profile<3*precision>) 
                                + 1 * sizeof(fixed_t<2*precision>)
                                + 1 * sizeof(fixed_t<3*precision>)
                                
                                + 1 * sizeof(fixed_t<3*precision>);
};






struct blind_vecs
{
  profile<precision> usr;
  profile<precision> itm;

  static constexpr size_t size = 2 * sizeof(profile<precision>);
};

struct blind_vecs_constant
{
  profile<precision> vec;
  fixed_t<precision> sca;

  static constexpr size_t size =  sizeof(fixed_t<precision>) + sizeof(profile<precision>);
};

template <size_t precision>
inline profile<2*precision> operator*(fixed_t<precision> c, const profile<precision> & p)
{
  profile<2*precision> cp;
  for (size_t i = 0; i < DIM; ++i) cp.f64[i] = c * p.f64[i];
  return cp;
}

template <size_t precision>
inline profile<precision> operator * (const profile<precision> & p, double c)
{
  profile<precision> cp;
  for (size_t i = 0; i < DIM; ++i) cp.f64[i] = (double)p.f64[i] * c;
  return cp;
}

template <size_t precision>
inline profile<3*precision> operator*(fixed_t<2*precision> c, const profile<precision> & p)
{
  profile<3*precision> cp;
  for (size_t i = 0; i < DIM; ++i) cp.f64[i] = p.f64[i] * c;
  return cp;
}


template <size_t precision>
inline profile<3*precision> operator*(fixed_t<precision> c, const profile<2*precision> & p)
{
  profile<3*precision> cp;
  for (size_t i = 0; i < DIM; ++i) cp.f64[i] = c * p.f64[i];
  return cp;
}


 // inline u64_array operator+(const u64_array & a1, const u64_array& a2)
 // {
 //  u64_array a3;

 //  for(size_t i = 0; i < DIM; ++i) a3.u64[i] = a1.u64[i] + a2.u64[i];

 //  return a3;
 // }

 // inline u64_array operator-(const u64_array & a1, const u64_array& a2)
 // {
 //  u64_array a3;

 //  for(size_t i = 0; i < DIM; ++i) a3.u64[i] = a1.u64[i] - a2.u64[i];

 //  return a3;
 // }

 //  inline u64_array operator*(const u64_array & a1, const u64_array& a2)
 // {
 //  u64_array a3;

 //  for(size_t i = 0; i < DIM; ++i) a3.u64[i] = a1.u64[i] * a2.u64[i];

 //  return a3;
 // }

template <size_t precision>
inline profile<precision> operator+(const profile<precision> & p1, const profile<precision> & p2)
{
  profile<precision> p3;
  for (size_t i = 0; i < dim256; ++i) p3.m256[i] = _mm256_add_epi64(p1.m256[i], p2.m256[i]);
  //for (size_t i = 0; i < DIM; ++i) p3.u64[i] = p1.u64[i] + p2.u64[i];
  return p3;
}

template <size_t precision>
inline profile<precision> operator-(const profile<precision> & p1, const profile<precision> & p2)
{
  profile<precision> p3;
  for (size_t i = 0; i < dim256; ++i) p3.m256[i] = _mm256_sub_epi64(p1.m256[i], p2.m256[i]);
//  for (size_t i = 0; i < DIM; ++i) p3.u64[i] = p1.u64[i] - p2.u64[i];
  return p3;
}

template <size_t precision>
inline profile<precision> & operator+=(profile<precision> & p1, const profile<precision> & p2)
{
  for (size_t i = 0; i < dim256; ++i) p1.m256[i] = _mm256_add_epi64(p1.m256[i], p2.m256[i]);
//  for (size_t i = 0; i < DIM; ++i) p1.u64[i] += p2.u64[i];
  return p1;
}

template <size_t precision>
inline bool operator==(profile<precision> & p1, const profile<precision> & p2)
{
  return memcmp(&p1, &p2, sizeof(profile<precision>)) == 0;
}

template <size_t precision>
inline bool operator!=(profile<precision> & p1, const profile<precision> & p2)
{
  return memcmp(&p1, &p2, sizeof(profile<precision>)) != 0;
}

template <size_t precision>
inline profile<precision> & operator-=(profile<precision> & p1, const profile<precision> & p2)
{
  for (size_t i = 0; i < dim256; ++i) p1.m256[i] = _mm256_sub_epi64(p1.m256[i], p2.m256[i]);
//  for (size_t i = 0; i < DIM; ++i) p1.u64[i] -= p2.u64[i];
  return p1;
}

template <size_t precision>
inline fixed_t<2*precision> dot(const profile<precision> & p1, const profile<precision> & p2)
{
  fixed_t<2*precision> dp = 0.0;
  for (size_t i = 0; i < DIM; ++i) dp += p1.f64[i] * p2.f64[i];
  return dp;
}

template <size_t precision>
inline fixed_t<precision> profile_dot3(profile<precision> * p1, const profile<precision> * p2, size_t n)
{
  fixed_t<precision> dp = 0.0;
  for (size_t i = 0; i < n; ++i)
  {

    for(size_t j = 0; j < DIM; ++j)
    {
      dp +=  p1[i].f64[j];
    }

  }

  return dp;
}

template <size_t precision>
std::ostream & operator<<(std::ostream & os, const profile<precision> & u)
{
  os << "[ ";
  for (size_t i = 0; i < DIM; ++i) os << u.f64[i] << " ";
  return os << "]";
}

template <size_t precision>
std::ostream & operator<<(std::ostream & os, const fixed_t<precision> & f)
{
  return os << (double)f;
}

#endif

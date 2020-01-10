

#include <type_traits>
#include <set>
#include <vector>
#include "dpf.h"
#include <iostream>
#include <assert.h>
#include "picosha2.h"
using namespace dpf;

typedef uint8_t leaf_t;
typedef __m256i node_t;
typedef LowMC<__m256i> prgkey_t;
typedef block<node_t> block_t;


  const prgkey_t prgkey;

  const block_t mask   = prgkey.mask;
  const block_t maska  = prgkey.maska;
  const block_t maskb  = prgkey.maskb;
  const block_t maskc  = prgkey.maskc;
  const block_t maskbc = prgkey.maskbc;


 const size_t round_ =  prgkey.rounds;
 const size_t rounds =  14;

inline bool xor_if(const bool & block1, const bool & block2, bool flag)
{
  if(flag) 
    {
    return (block1 ^ block2);
    }
  else
  {
    return block1;
  }  
}

class MPCrandomness
{
  
  public:
  
    MPCrandomness(AES_KEY& aeskey, __m128i seed, size_t len) :
    buf((unsigned char *)malloc(len * sizeof(seed))),
    cur(buf)
    {
      PRG(aeskey, seed, (__m128i *)buf, len);
    } 
    ~MPCrandomness()
    {
    free(buf);
    }
  



  inline bool next(int & val)
  {
    memcpy(&val, cur, sizeof(int));
    val &= 1;
    cur += sizeof(int);
    bool val_b = val;

    //std::cout << "val_b = " << val_b << std::endl;

    return val_b;
  }

  inline bool next_bool()
  {
    int  val;
    memcpy(&val, cur, sizeof(int));
    val &= 1;
    cur += sizeof(int);
    bool val_b = val;

    //std::cout << "val_b = " << val_b << std::endl;

    return val_b;
  }

  template<typename T>
  inline T& next(T & val)
  {
    memcpy(&val, cur, sizeof(T));
    cur += sizeof(T);
    return val;
  }
  
  inline block_t next_block()
  {
    block_t val;
    memcpy(&val, cur, sizeof(block_t));
    cur += sizeof(block_t);
    return val;
  }
  

  template<typename T>
  inline T && next()
  {
    T val;
    memcpy(&val, cur, sizeof(T));
    cur += sizeof(T);
    return val;
  }


  private:
  unsigned char * buf;
  unsigned char * cur;
};



#include "transcripts.h"
#include "simulator.h"
#include "verifier.h"

int main(int argc, char * argv[])
{ 


  std::cout << "rounds = " << rounds << std::endl;

  bool party0 = false; bool party1 = true;
  


  const size_t nitems = 1ULL << 44;
  const size_t target = atoi(argv[1]);
 
  const leaf_t val = 0x12;//_mm_set1_epi8(0x12);

  auto [dpfkey0, dpfkey1] = dpf_key<leaf_t, node_t, prgkey_t>::gen(prgkey, nitems, target, val);
  
  const size_t depth = dpfkey1.depth(nitems);  
 
  std::bitset<depth> directions = target;
  
  for(std::size_t i = 0; i < depth/2; ++i) 
  {
    bool t = directions[i];
    directions[i] = directions[depth-i-1];
    directions[depth-i-1] = t;
  }
  
  std::cout << "directions = " << directions << std::endl;
  
  AES_KEY aeskey;
  
  from_P2 from_P2_to_P0(depth, rounds), from_P2_to_P1(depth, rounds);
  from_PB from_P0_original(depth, rounds), from_P1_original(depth, rounds);
  from_PB from_P0_decompressed(depth, rounds);
  from_PB from_P1_decompressed(depth, rounds);

  from_P2 from_P2_to_PB(depth, rounds);
  from_PB from_PB_other(depth, rounds);
  
  size_t len = 100000;
  
  __m128i seed0, seed1, seed2;
  
  arc4random_buf(&seed0, sizeof(__m128i));
  arc4random_buf(&seed1, sizeof(__m128i));
  arc4random_buf(&seed2, sizeof(__m128i));
  

  Simulator sim(aeskey, seed0, seed1, seed2, len, depth);  
  Verifier  ver0(aeskey, seed0, len, depth);
  Verifier  ver1(aeskey, seed1, len, depth);

  sim.root_layer(from_P0_original, from_P1_original, from_P2_to_P0, from_P2_to_P1, prgkey, dpfkey0, dpfkey1);

  std::cout << "dept = " << depth << std::endl;

  for(size_t index = 1; index < depth-1; ++index)
  {
   sim.middle_layers(from_P0_original, from_P1_original, from_P2_to_P0, from_P2_to_P1, prgkey, dpfkey0, dpfkey1, index);
  }

  from_PB_compressed from_P0_compressed, from_P1_compressed;

  from_P0_compressed = compressTranscript(from_P0_original);
  from_P1_compressed = compressTranscript(from_P1_original);
  
  from_P0_decompressed = decompressTranscript(from_P0_compressed);
  from_P1_decompressed = decompressTranscript(from_P1_compressed);
  
  ver0.Pdirection = sim.P0direction;
  
  ver0.root_layer(from_P2_to_P0, from_P1_decompressed, from_PB_other, prgkey, dpfkey0, party0);
  
  assert(from_PB_other.L_shares_recv == from_P0_decompressed.L_shares_recv);
  assert(from_PB_other.R_shares_recv == from_P0_decompressed.R_shares_recv);
  assert(from_PB_other.bit_L_shares_recv == from_P0_decompressed.bit_L_shares_recv);
  assert(from_PB_other.bit_R_shares_recv == from_P0_decompressed.bit_R_shares_recv);
  
  for(size_t j = 0; j < 4; ++j)
  {
   assert(from_PB_other.next_bit_L_recv[0][j] == from_P0_decompressed.next_bit_L_recv[0][j]);
   assert(from_PB_other.next_bit_R_recv[0][j] == from_P0_decompressed.next_bit_R_recv[0][j]);
   assert(from_PB_other.blinds_recv[0][j] == from_P0_decompressed.blinds_recv[0][j]);   
   assert(from_PB_other.next_bit_L_recv2[0][j] == from_P0_decompressed.next_bit_L_recv2[0][j]);
   assert(from_PB_other.next_bit_R_recv2[0][j] == from_P0_decompressed.next_bit_R_recv2[0][j]);
  }

   for(size_t index = 1; index < depth-1; ++index)
   { 
       ver0.middle_layers(from_P2_to_P0, from_P1_decompressed, from_PB_other, prgkey, ver0.seed0[index], ver0.seed1[index],   dpfkey0, index, party0);
   }

  for(size_t i = 1; i < depth-1; ++i)
  {
    for(size_t j = 0; j < 4; ++j)
    {
     assert(from_PB_other.next_bit_L_recv[i][j] == from_P0_decompressed.next_bit_L_recv[i][j]);
     assert(from_PB_other.next_bit_R_recv[i][j] == from_P0_decompressed.next_bit_R_recv[i][j]);
     assert(from_PB_other.blinds_recv[i][j] == from_P0_decompressed.blinds_recv[i][j]);   
     assert(from_PB_other.next_bit_L_recv2[i][j] == from_P0_decompressed.next_bit_L_recv2[i][j]);
     assert(from_PB_other.next_bit_R_recv2[i][j] == from_P0_decompressed.next_bit_R_recv2[i][j]);

    }

    for(size_t j = 0; j < rounds; ++j)
    {
      assert(from_PB_other.seed0R_encrypt[i][j] == from_P0_decompressed.seed0R_encrypt[i][j]);
    } 
   }

  std::string src = from_PB_other.to_string();
  std::string src2 = from_P0_decompressed.to_string();
  std::cout << "hash: (verifier 0) " << picosha2::hash256_hex_string(src) << "\n" << std::endl;
  std::cout << "hash: (verifier 0) " << picosha2::hash256_hex_string(src2) << "\n" << std::endl;

  ver1.Pdirection = sim.P1direction;

  std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << " ------------------------------------------------------------ " << std::endl << std::endl;

  ver1.root_layer(from_P2_to_P1, from_P0_decompressed, from_PB_other,   prgkey, dpfkey1, party1);

  assert(from_PB_other.L_shares_recv == from_P1_decompressed.L_shares_recv);
  assert(from_PB_other.R_shares_recv == from_P1_decompressed.R_shares_recv);
  assert(from_PB_other.bit_L_shares_recv == from_P1_decompressed.bit_L_shares_recv);
  assert(from_PB_other.bit_R_shares_recv == from_P1_decompressed.bit_R_shares_recv);  
 
  for(size_t j = 0; j < 4; ++j)
  {
   assert(from_PB_other.next_bit_L_recv[0][j] == from_P1_decompressed.next_bit_L_recv[0][j]);
   assert(from_PB_other.next_bit_R_recv[0][j] == from_P1_decompressed.next_bit_R_recv[0][j]);
   assert(from_PB_other.blinds_recv[0][j] == from_P1_decompressed.blinds_recv[0][j]);   
   assert(from_PB_other.next_bit_L_recv2[0][j] == from_P1_decompressed.next_bit_L_recv2[0][j]);
   assert(from_PB_other.next_bit_R_recv2[0][j] == from_P1_decompressed.next_bit_R_recv2[0][j]);
  }

  for(size_t index = 1; index < depth-1; ++index)
  {  
    ver1.middle_layers(from_P2_to_P1, from_P0_decompressed, from_PB_other, prgkey, ver1.seed0[index], ver1.seed1[index], dpfkey1, index, party1);
  }

 
  for(size_t i = 1; i < depth-1; ++i)
  {
    for(size_t j = 0; j < 4; ++j)
    {
     assert(from_PB_other.next_bit_L_recv[i][j] == from_P1_decompressed.next_bit_L_recv[i][j]);
     assert(from_PB_other.next_bit_R_recv[i][j] == from_P1_decompressed.next_bit_R_recv[i][j]);
     assert(from_PB_other.blinds_recv[i][j] == from_P1_decompressed.blinds_recv[i][j]);   
     assert(from_PB_other.next_bit_L_recv2[i][j] == from_P1_decompressed.next_bit_L_recv2[i][j]);
    } 

    for(size_t j = 0; j < rounds; ++j)
    {
      assert(from_PB_other.seed1R_encrypt[i][j] == from_P1_decompressed.seed1R_encrypt[i][j]);
    }

  }
 src = from_PB_other.to_string();
 src2 = from_P1_decompressed.to_string();
std::cout << "hash: (verifier 1) " << picosha2::hash256_hex_string(src) << "\n" << std::endl;
std::cout << "hash: (verifier 1) " << picosha2::hash256_hex_string(src2) << "\n" << std::endl;

  return 0;
}

 
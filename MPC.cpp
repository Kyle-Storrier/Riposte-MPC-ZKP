

#include <type_traits>
#include <set>
#include <vector>
#include "dpf.h"
#include <iostream>
#include <assert.h>
#include "picosha2.h"

 
using namespace dpf;

typedef uint16_t leaf_t;
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




#include "randomness.h"
#include "transcripts.h"
#include "simulator.h"
#include "verifier.h"
#include "verifier2.h"

// size_t number_of_simulations = 10;

 
int main(int argc, char * argv[])
{  
  std::cout << "rounds = " << rounds << std::endl;

  bool party0 = false; bool party1 = true;
  


  const size_t nitems = 1ULL << 40;
  const size_t target =  atoi(argv[1]);
 
  const leaf_t val =  3333; //_mm256_set1_epi8(0x01);

 
  auto [dpfkey0, dpfkey1] = dpf_key<leaf_t, node_t, prgkey_t>::gen(prgkey, nitems, target, val);
  std::cout << "leaves per node: "  << dpfkey0.leaves_per_node << std::endl;
  const size_t depth = dpfkey1.depth(nitems);  
 
  from_PB_ from_P0_[depth], from_P1_[depth], from_PB_other_[depth], from_PB_other_2[depth];
  memset(from_P0_, 0, sizeof(from_P0_[0]) * depth);
  memset(from_P1_, 0, sizeof(from_P1_[1]) * depth);
  memset(from_PB_other_, 0, sizeof(from_P0_[0]) * depth);
  memset(from_PB_other_2, 0, sizeof(from_P1_[0]) * depth);
  from_PB_compressed_ from_P0_compressed[depth], from_P1_compressed[depth];
  memset(from_P0_compressed, 0, sizeof(from_PB_compressed_) * depth);
  memset(from_P1_compressed, 0, sizeof(from_PB_compressed_) * depth);
  from_PB_ from_P0_decompressed[depth], from_P1_decompressed[depth];
  memset(from_P0_decompressed, 0, sizeof(from_P0_[0]) * depth);
  memset(from_P1_decompressed, 0, sizeof(from_P1_[0]) * depth);

  from_PB_ from_P0_V0_[depth], from_P1_V0_[depth];

  std::cout << "depth = " << depth << std::endl;
  std::bitset<depth> directions = ceil(target/dpfkey0.leaves_per_node);
  
  for(std::size_t i = 0; i < depth/2; ++i) 
  {
    bool t = directions[i];
    directions[i] = directions[depth-i-1];
    directions[depth-i-1] = t;
  } 
  
  std::bitset<depth>  P0direction = rand();
  std::bitset<depth>  P1direction =  P0direction ^ directions;
  std::cout << "directions = " << directions << std::endl;
   
  AES_KEY aeskey;
  
  from_P2 from_P2_to_P0(depth, rounds), from_P2_to_P1(depth, rounds);
  //from_PB from_P0(depth, rounds), from_P1(depth, rounds);

  from_P2 from_P2_to_PB(depth, rounds);
 // from_PB from_PB_other(depth, rounds);
  

  from_P2 from_P2_to_P0_V0(depth, rounds), from_P2_to_P1_V0(depth, rounds);
 // from_PB from_P0_V0(depth, rounds), from_P1_V0(depth, rounds);

  size_t len = 100000;
  
  __m128i seed0, seed1, seed2;

  
//for(size_t ii = 0; ii < number_of_simulations; ++ii)
{  
//  std::cout << "ii = " << ii << std::endl;
  arc4random_buf(&seed0, sizeof(__m128i));
  arc4random_buf(&seed1, sizeof(__m128i));
  arc4random_buf(&seed2, sizeof(__m128i));

  Simulator sim(aeskey, seed0, seed1, seed2, len, depth);  
  for(size_t i = 0; i < depth; ++i)
  {
    sim.P0direction[i] = P0direction[i];
    sim.P1direction[i] = P1direction[i];
  }

  sim.root_layer(from_P0_, from_P1_, from_P2_to_P0, from_P2_to_P1, prgkey, dpfkey0, dpfkey1);

  for(size_t index = 1; index < depth; ++index)
  {
   sim.middle_layers(from_P0_, from_P1_, from_P2_to_P0, from_P2_to_P1, prgkey, dpfkey0, dpfkey1, index);
  }  

  std::cout << "--------------------------------------" << std::endl << std::endl << std::endl;
}

  block_t compressionClearBits;
  compressionClearBits.bits.reset();
  for(long i = 0; i < 256 - 66; i++) { // Using 256 - 66 is due to the shifts applied to the masks in the LowMC code.
    compressionClearBits.bits.set(i);
  }

for(size_t i = 0; i < depth; i++) {
  compressTranscript(from_P0_[i], from_P0_compressed[i]);
  compressTranscript(from_P1_[i], from_P1_compressed[i]);
}

for(size_t i = 0; i < depth; i++) {
  decompressTranscript(from_P0_compressed[i], from_P0_decompressed[i]);
  decompressTranscript(from_P1_compressed[i], from_P1_decompressed[i]);
}

  Verifier2 ver2_0(aeskey, seed0, seed1, seed2, len, depth);  

  for(size_t i = 0; i < depth; ++i)
  {
    ver2_0.P0direction[i] = P0direction[i];
    ver2_0.P1direction[i] = P1direction[i];
  }

  ver2_0.root_layer(from_P0_V0_, from_P1_V0_, from_P2_to_P0_V0, from_P2_to_P1_V0, prgkey, dpfkey0, dpfkey1);
 
  for(size_t index = 1; index < depth; ++index)
  {
   ver2_0.middle_layers(from_P0_V0_, from_P1_V0_, from_P2_to_P0_V0, from_P2_to_P1_V0, prgkey, dpfkey0, dpfkey1, index);
  }
  

  Verifier  ver0(aeskey, seed0, len, depth);
  Verifier  ver1(aeskey, seed1, len, depth);
 
 

  for(size_t i = 0; i < depth; ++i) ver0.Pdirection[i] = ver2_0.P0direction[i];
  
  ver0.root_layer(from_P2_to_P0, from_P1_decompressed, from_PB_other_, prgkey, dpfkey0, party0);
  
  assert(from_PB_other_[0].L_shares_recv == from_P0_decompressed[0].L_shares_recv);
  assert(from_PB_other_[0].R_shares_recv == from_P0_decompressed[0].R_shares_recv);
  assert(from_PB_other_[0].bit_L_shares_recv == from_P0_decompressed[0].bit_L_shares_recv);
  assert(from_PB_other_[0].bit_R_shares_recv == from_P0_decompressed[0].bit_R_shares_recv);
  
  for(size_t index = 1; index < depth; ++index)
  { 
       ver0.middle_layers(from_P2_to_P0 , from_P1_decompressed, from_PB_other_, prgkey, ver0.seed0[index], ver0.seed1[index],   dpfkey0, index, party0);
  }

   for(size_t i = 0; i < depth-1; ++i)
  {
    for(size_t j = 0; j < rounds; ++j)
    {
      from_PB_other_[i].seed0L_encrypt[j] = from_PB_other_[i].seed0L_encrypt[j] & compressionClearBits;
      from_PB_other_[i].seed0R_encrypt[j] = from_PB_other_[i].seed0R_encrypt[j] & compressionClearBits;
      from_PB_other_[i].seed1L_encrypt[j] = from_PB_other_[i].seed1L_encrypt[j] & compressionClearBits;
      from_PB_other_[i].seed1R_encrypt[j] = from_PB_other_[i].seed1R_encrypt[j] & compressionClearBits;
    }
  }
  
  for(size_t i = 0; i < depth-1; ++i)
  {
      //std::cout << "i = " << i << std::endl;
      for(size_t j = 0; j < 4; ++j)
      { 
       //std::cout << "-> " << (from_PB_other_[i].blinds_recv[j] ^ from_P0_[i].blinds_recv[j]).bits << std::endl;
       assert(from_PB_other_[i].next_bit_L_recv[j] == from_P0_decompressed[i].next_bit_L_recv[j]);
       assert(from_PB_other_[i].next_bit_R_recv[j] == from_P0_decompressed[i].next_bit_R_recv[j]);
       assert(from_PB_other_[i].blinds_recv[j] == from_P0_decompressed[i].blinds_recv[j]);   
       //assert(from_PB_other_[i].next_bit_L_recv2[j] == from_P0_[i].next_bit_L_recv2[j]);
       //assert(from_PB_other_[i].next_bit_R_recv2[j] == from_P0_[i].next_bit_R_recv2[j]);
      }

      for(size_t j = 0; j < rounds; ++j)
      {
        std::cout << "i = " << i << " j = " << j << " depth = " << depth << std::endl;
        assert(from_PB_other_[i].seed0L_encrypt[j] == from_P0_decompressed[i].seed0L_encrypt[j]);
        assert(from_PB_other_[i].seed0R_encrypt[j] == from_P0_decompressed[i].seed0R_encrypt[j]);
        assert(from_PB_other_[i].seed1L_encrypt[j] == from_P0_decompressed[i].seed1L_encrypt[j]);
        assert(from_PB_other_[i].seed1R_encrypt[j] == from_P0_decompressed[i].seed1R_encrypt[j]);
      } 
  }




  const char * buffer0 = (const  char*)& from_P0_decompressed[2];
  const char * buffer_other_0 = (const  char*)&from_PB_other_[2];
  std::string str_other_0(buffer_other_0);
  std::string str0(buffer0);
  std::cout << "str_other_0 = " << picosha2::hash256_hex_string(str_other_0) << std::endl;
  std::cout << "str0        = " << picosha2::hash256_hex_string(str0) << std::endl;
  
   ver1.Pdirection = ver2_0.P1direction;

   std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << " ------------------------------------------------------------ " << std::endl << std::endl;

   ver1.root_layer(from_P2_to_P1 , from_P0_decompressed , from_PB_other_2,   prgkey, dpfkey1, party1);

 
 
 


  for(size_t index = 1; index < depth; ++index)
  {  
    ver1.middle_layers(from_P2_to_P1 , from_P0_decompressed , from_PB_other_2 , prgkey, ver1.seed0[index], ver1.seed1[index], dpfkey1, index, party1);
  }

    for(size_t i = 0; i < depth-1; ++i)
  {
    for(size_t j = 0; j < rounds; ++j)
    {
      from_PB_other_2[i].seed0L_encrypt[j] = from_PB_other_2[i].seed0L_encrypt[j] & compressionClearBits;
      from_PB_other_2[i].seed0R_encrypt[j] = from_PB_other_2[i].seed0R_encrypt[j] & compressionClearBits;
      from_PB_other_2[i].seed1L_encrypt[j] = from_PB_other_2[i].seed1L_encrypt[j] & compressionClearBits;
      from_PB_other_2[i].seed1R_encrypt[j] = from_PB_other_2[i].seed1R_encrypt[j] & compressionClearBits;
    }
  }

  for(size_t i = 0; i < depth-1; ++i)
  {
    assert(from_PB_other_2[i].L_shares_recv == from_P1_decompressed[i].L_shares_recv);
    assert(from_PB_other_2[i].R_shares_recv == from_P1_decompressed[i].R_shares_recv);
    assert(from_PB_other_2[i].final == from_P1_decompressed[i].final);
    assert(from_PB_other_2[i].bit_L_shares_recv == from_P1_decompressed[i].bit_L_shares_recv);
    assert(from_PB_other_2[i].bit_R_shares_recv == from_P1_decompressed[i].bit_R_shares_recv);
    for(size_t j = 0; j < 4; ++j)
    { 
      assert(from_PB_other_2[i].next_bit_L_recv[j] == from_P1_decompressed[i].next_bit_L_recv[j]);
      assert(from_PB_other_2[i].next_bit_R_recv[j] == from_P1_decompressed[i].next_bit_R_recv[j]);
      assert(from_PB_other_2[i].blinds_recv[j] == from_P1_decompressed[i].blinds_recv[j]);   
      assert(from_PB_other_2[i].bit_blinds_recv[j] == from_P1_decompressed[i].bit_blinds_recv[j]);   
      assert(from_PB_other_2[i].next_bit_L_recv2[j] == from_P1_decompressed[i].next_bit_L_recv2[j]);
      assert(from_PB_other_2[i].next_bit_R_recv2[j] == from_P1_decompressed[i].next_bit_R_recv2[j]);
    }

    for(size_t j = 0; j < rounds; ++j)
    {
      assert(from_PB_other_2[i].seed0L_encrypt[j] == from_P1_decompressed[i].seed0L_encrypt[j]);
      assert(from_PB_other_2[i].seed0R_encrypt[j] == from_P1_decompressed[i].seed0R_encrypt[j]);
      assert(from_PB_other_2[i].seed1L_encrypt[j] == from_P1_decompressed[i].seed1L_encrypt[j]);
      assert(from_PB_other_2[i].seed1R_encrypt[j] == from_P1_decompressed[i].seed1R_encrypt[j]);
    } 
  }
 
 
 
  const char * buffer_other = (const  char*)&from_PB_other_2[2];

  const char * buffer1 = (const  char*)&from_P1_decompressed[2];
 
 
 
  std::string str_other(buffer_other);
  std::string str1(buffer1);    
 
  std::cout << "check0 = " << picosha2::hash256_hex_string(str_other) << std::endl;
  std::cout << "check1 = " << picosha2::hash256_hex_string(str1) << std::endl;
  
 
 
  char * transcript = reinterpret_cast<char* >(from_PB_other_2);
  unsigned char hashed[picosha2::k_digest_size];
  picosha2::hash256(transcript, transcript + depth * sizeof(from_PB_other_2[0]), hashed, hashed + picosha2::k_digest_size);
 
  char * transcript1 = reinterpret_cast<char* >(from_P1_decompressed);
  unsigned char hashed1[picosha2::k_digest_size];
  picosha2::hash256(transcript1, transcript1 + depth * sizeof(from_P1_decompressed[0]), hashed1, hashed1 + picosha2::k_digest_size);
// TODO (kyle): Problems present getting this to work properly with compression.

bool match = true;

for (int i = 0; i < 32; i++) {
  printf("%x", hashed[i]);
  match = match && (hashed[i] == hashed1[i]);
}
 std::cout << std::endl << "-- " << std::endl;
  for (int i = 0; i < 32; i++) {
  printf("%x", hashed1[i]);
}
std::cout << std::endl;

if(match) {
  std::cout << "SUCCESS: The hashes match\n";
} else {
  std::cout << "FAILURE: The hashes don't match\n";
}
 // std::cout << "hash: (verifier 1) " << picosha2::hash256_hex_string(src) << "\n" << std::endl;
 // std::cout << "hash: (verifier 1) " << picosha2::hash256_hex_string(src2) << "\n" << std::endl;

  return 0;
}

 
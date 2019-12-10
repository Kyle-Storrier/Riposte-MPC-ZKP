#include <unistd.h>
#include <fcntl.h>
#include <tuple>
#include <chrono>
#include "dpf.h"
#include <string>
#include <assert.h>

using namespace std::chrono;
__m256i all_ones(void) { return _mm256_set1_epi64x(-1); }

const size_t nitems   = 1ULL << 20;
typedef __m256i __mX;
dpf_key<__mX, nitems> dpfkey[2] = { 0 }; 

const size_t depth = dpfkey[0].depth;

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

  template<typename T>
  inline T& next(T & val)
  {
    memcpy(&val, cur, sizeof(T));
    cur += sizeof(T);
    return val;
  }
  
  inline block next_block()
  {
    block val;
    memcpy(&val, cur, sizeof(block));
    cur += sizeof(block);
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


  // template<typename T>
  // inline  && T next()
  // {
  // T val;
  // next(val);
  // return val;
  // }
  
  private:
  unsigned char * buf;
  unsigned char * cur;
};


#include "transcripts.h"
P2_Transcripts P2_message[2];
PB_Transcripts P0_message, P1_message;

#include "simulator.h"
#include "verifier.h"



  void Simulator::dpf_mpc(LowMC& lowmckey, AES_KEY& aeskey, dpf_key<__mX, nitems>& k0, dpf_key<__mX, nitems> k1)
  {
    std::cout << "Simulator::dpf_mpc" << std::endl;

    const block root0 = k0.root; 
    const block root1 = k1.root;

    bool t0[2];
    block s0[2];

    bool t1[2];
    block s1[2];
 
    block share[2];
    arc4random_buf(share, 2 * sizeof(block));
    share[0] = k0.root ^ share[1];

    block gamma[2][rounds];
    block blinds0[rounds];
    block blinds1[rounds];
    block rand[rounds];
    
    size_t buflen =  rounds * sizeof(block);    
    
     for (unsigned r = 0; r < rounds; ++r)
     {
        // P0rand.next(blinds0[r]);
        // P1rand.next(blinds1[r]); 
        // P2rand.next(rand[r]);

         arc4random_buf(&blinds0[r], sizeof(block));
         arc4random_buf(&blinds1[r], sizeof(block));
         arc4random_buf(&rand[r],  sizeof(block));

        const block tmp1 = ((blinds0[r] >> 1) & blinds1[r]) ^ ((blinds1[r] >> 1) & blinds0[r]);
        const block tmp2 = ((blinds0[r] >> 2) & blinds1[r]) ^ ((blinds1[r] >> 2) & blinds0[r]);
    
        const block bc = (tmp1 << 2) & maska;
        const block ac = (tmp2 << 1) & maskb;
        const block ab = (tmp1 >> 1) & maskc;
    
        gamma[0][r] = (bc | ac | ab) ^ rand[r];
        gamma[1][r] = rand[r];// ^ roundkeysXORconstants[r+1];
     }

     const block c0 = share[0]; const block c1 = share[1];
    std::pair<std::vector<block>,std::vector<block>> encrypt_outL = lowmckey.encrypt_MPC_proof(c0, c1, blinds0, blinds1, gamma);
    std::vector<block> verify_out = lowmckey.encrypt_MPC_verify(c0, encrypt_outL.second, blinds0, gamma[0], false);
    
    for(size_t j = 0; j < 10; ++j)
    {
      std::cout << j << " : " << (verify_out[j] ^ encrypt_outL.first[j]) << std::endl;
    }

    //std::pair<std::vector<block>,std::vector<block>> expand_out = expand_mpc(lowmckey, c0, c1,  s0, t0, s1, t1, blinds0, blinds1, gamma);
    //expanded_out (first and second are the two shares of expand)
    const block m = c0 ^ c1;
   
    block s[2]; bool t[2];

    expand_nonmpc(lowmckey, m , s , t);

//    std::cout << "zero: " << ((expand_out.first[L] ^ expand_out.second[L]) ^ s[0]) << " " << ((expand_out.first[R] ^ expand_out.second[R]) ^ s[1]) << std::endl;


    // const block ss0 = expand_out.first[R];
    // const block ss1 = expand_out.second[R];

    // block results[2];
    // multiply_mpc(ss0 , t0, ss1 , !t0, results, 0, 0);
   
    // std::cout << "non-mpc: " << (ss0 ^ ss1) << std::endl;

    //std::cout << std::endl << "expand: " << s[0] << " " << s[1] << std::endl;
  }







int main(int argc, char ** argv)
{  

  uint64_t point = std::stoi(argv[1]);
  std::cout << "point = " << point << std::endl;
  std::bitset<depth> directions = point;
      for(std::size_t i = 0; i < depth/2; ++i) {
        bool t = directions[i];
        directions[i] = directions[depth-i-1];
        directions[depth-i-1] = t;
    }

  std::cout << "directions = " << directions << std::endl;
  AES_KEY aeskey;
  LowMC lowmckey;

  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);



// __m128i seed;
// arc4random_buf(&seed, sizeof(__m128i));
// size_t len = 10;
// MPCrandomness P0(aeskey, seed, len);
// block x[20];
// P0.next(x[0]);


  gen(lowmckey, point , dpfkey);

  __m128i seed0, seed1, seed2;

  arc4random_buf(&seed0, sizeof(__m128i));
  arc4random_buf(&seed1, sizeof(__m128i));
  arc4random_buf(&seed2, sizeof(__m128i));

  size_t len = 100000;

  Simulator sim(aeskey, seed0, seed1, seed2, len);  
  Verifier  ver0(aeskey, seed0, len);  
  Verifier  ver1(aeskey, seed1, len);

  sim.P0direction = rand();
  sim.P1direction = sim.P0direction ^ directions;

  std::cout << "P0direction = " << sim.P0direction << std::endl;
  std::cout << "P1direction = " << sim.P1direction << std::endl;

  std::cout << "P0direction[depth] = " << sim.P0direction[depth-1] << std::endl;
  std::cout << "P1direction[depth] = " << sim.P1direction[depth-1] << std::endl;

  block block0, block1;
  bool b0 , b1;


  //sim.dpf_mpc(lowmckey, aeskey, dpfkey[0], dpfkey[1]);

  sim.root_layer(lowmckey, dpfkey);
  
  for(size_t index = 1; index < depth-1; ++index)
  {
   sim.middle_layers(lowmckey, dpfkey, index);
  }

  bool party0 = false; bool party1 = true;
  

  std::cout << std::endl << std::endl << "=======================================VERIFIER====================================================" << std::endl << std::endl << std::endl << std::endl;
 
  
  ver0.Pdirection = sim.P0direction;

  
  ver0.root_layer(P0_message, lowmckey, dpfkey[0], party0);

 
  
   for(size_t index = 1; index < depth-1; ++index)
   { 
  //  // std::cout << std::endl << std::endl << "--------------------------------------------" << std::endl << std::endl;
      ver0.middle_layers(lowmckey, ver0.seed0[index][0], ver0.seed1[index][0], P0_message, dpfkey[0], index, party0);
   }


  ver1.Pdirection = sim.P1direction;

  std::cout << std::endl << std::endl << std::endl << " ------------------------------------------------------------ " << std::endl << std::endl;

  ver1.root_layer(P1_message, lowmckey, dpfkey[1], party1);
  for(size_t index = 1; index < 2; ++index)
  { 
   // std::cout << std::endl << std::endl << "--------------------------------------------" << std::endl << std::endl;
    ver1.middle_layers(lowmckey, ver1.seed0[index][0], ver1.seed1[index][0], P1_message, dpfkey[1], index, party1);
  }

  return 0;
}

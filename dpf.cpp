#include <unistd.h>
#include <fcntl.h>
//#include <iostream>
#include <tuple>
#include <chrono>
//#include "common.h"
#include "dpf.h"
//#include "block.h"

using namespace std::chrono;


    

int main(int argc, char ** argv)
{

  #ifdef LOWMC
    typedef __m256i __mX;
    LowMC key;
  #endif

  #ifdef AES 
    typedef __m128i __mX;
    AES_KEY key;
    AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &key);
  #endif

  const size_t nitems   =  512;

  

  blocks<__mX> * s0 , * s1;
  uint8_t * t0, * t1;
  
  t0 = (uint8_t*)malloc(nitems * sizeof(uint8_t));
  t1 = (uint8_t*)malloc(nitems * sizeof(uint8_t));
  if (posix_memalign((void**)&s0, sizeof(blocks<__mX>), nitems * sizeof(__mX)))
  {
   throw std::runtime_error("posix_memalign failed");
  }

  if (posix_memalign((void**)&s1, sizeof(blocks<__mX>), nitems * sizeof(__mX)))
  {
   throw std::runtime_error("posix_memalign failed");
  }

  dpf_key<__mX, nitems> dpfkey[2] = { 0 }; 
  
  gen(key, 2 , dpfkey);
  
  evalfull3(key, dpfkey[0], s0, t0);
  evalfull3(key, dpfkey[1], s1, t1);

 for(size_t j = 0; j < nitems; ++j)
 {
   std::cout << (double)t0[j] << " " << (double)t1[j] << std::endl;
   std::cout << blocks<__mX>(s0[j] ^ s1[j]) << std::endl;
 }

  return 0;
}

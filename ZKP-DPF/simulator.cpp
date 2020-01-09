#include "simulator.h"

  // void Simulator::Gen_Blinds(block_t blinds0[rounds], block_t blinds1[rounds], block_t gamma[2][rounds])
  // {
  //   block_t rand[rounds];
  //   for (unsigned r = 0; r < rounds; ++r)
  //    {        
  //       blinds0[r] = P0rand.next_block();
  //       blinds1[r] = P1rand.next_block();
  //       rand[r]    = P2rand.next_block();

  //       const block_t tmp1 = ((blinds0[r] >> 1) & blinds1[r]) ^ ((blinds1[r] >> 1) & blinds0[r]);
  //       const block_t tmp2 = ((blinds0[r] >> 2) & blinds1[r]) ^ ((blinds1[r] >> 2) & blinds0[r]);
    
  //       const block_t bc = (tmp1 << 2) & maska;
  //       const block_t ac = (tmp2 << 1) & maskb;
  //       const block_t ab = (tmp1 >> 1) & maskc;
    
  //       gamma[0][r] = (bc | ac | ab) ^ rand[r];
  //       gamma[1][r] = rand[r];
  //    }
  // }
  

//   void Simulator::multiply_mpc_b(const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth)
//   {
//     bool D0 = P0rand.next_bool(); 
//     bool D1 = P1rand.next_bool(); 
     
//     bool d0 = P0rand.next_bool();     
//     bool d1 = P1rand.next_bool();
      
//     const bool alpha = P2rand.next_bool(); 

//     bool c0 = xor_if(alpha, D0, d1);
//     bool c1 = xor_if(alpha, D1, d0);
 
//     from_P2_to_P0.c_bit[cur_depth][mul]       = c0;
//     from_P2_to_P1.c_bit[cur_depth][mul]       = c1;
 
//     from_P0.next_bit_L_recv[cur_depth][mul] = x0 ^ D0;
//     from_P0.next_bit_R_recv[cur_depth][mul] = b0 ^ d0;

//     from_P1.next_bit_L_recv[cur_depth][mul] = x1 ^ D1;
//     from_P1.next_bit_R_recv[cur_depth][mul] = b1 ^ d1;

//     const bool gamma0 = xor_if(c0, from_P1.next_bit_L_recv[cur_depth][mul], d0);
//     bool xx0 = x0;   
 

//     const bool gamma1 = xor_if(c1, from_P0.next_bit_L_recv[cur_depth][mul], d1);
//     bool xx1 = x1;  
  
//     results[0] = xor_if(gamma0, xx0, (b0 ^ from_P1.next_bit_R_recv[cur_depth][mul]));
  
//     results[1] = xor_if(gamma1, xx1, (b1 ^ from_P0.next_bit_R_recv[cur_depth][mul])); 
//   }


//   void Simulator::get_next_bits(bool bit0_L[2], bool bit0_R[2], bool bit1_L[2], bool bit1_R[2], size_t cur_depth)
//   {
    
//     bool notb_bit0_L[2];    
//     multiply_mpc_b(bit0_L[0], P0direction[cur_depth], bit0_L[1], !P1direction[cur_depth], notb_bit0_L, 0, cur_depth); // L_0 * (direction - 1)


//     bool notb_bit1_L[2]; 
//     multiply_mpc_b(bit1_L[0], P0direction[cur_depth], bit1_L[1], !P1direction[cur_depth], notb_bit1_L, 1, cur_depth); // L_1 * (direction - 1)

//     bool b_bit0_R[2];    
//     multiply_mpc_b(bit0_R[0], P0direction[cur_depth], bit0_R[1],  P1direction[cur_depth], b_bit0_R, 2, cur_depth);

//     bool b_bit1_R[2];
//     multiply_mpc_b(bit1_R[0], P0direction[cur_depth], bit1_R[1],  P1direction[cur_depth], b_bit1_R, 3, cur_depth); // R_1 * direction    
    
    
//     bit0[cur_depth+1][0] = notb_bit0_L[0] ^ b_bit0_R[0];
    
//     bit0[cur_depth+1][1] = notb_bit0_L[1] ^ b_bit0_R[1];    
    
//     bit1[cur_depth+1][0] = notb_bit1_L[0] ^ b_bit1_R[0];    
    
//     bit1[cur_depth+1][1] = notb_bit1_L[1] ^ b_bit1_R[1];

//   }


//     void Simulator::multiply_mpc(const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth)
//   {

//     //std::cout << "Simulator:::multiply_mpc, cur_depth = " << cur_depth << std::endl; 

//     bool D0 = P0rand.next_bool(); 
//     bool D1 = P1rand.next_bool(); 
     
//     bool d0 = P0rand.next_bool();     
//     bool d1 = P1rand.next_bool();
    


//     const bool alpha = P2rand.next_bool(); 

//     bool c0 = xor_if(alpha, D0, d1);
//     bool c1 = xor_if(alpha, D1, d0);

//     from_P2_to_P0.c_bit2[cur_depth][mul] = c0;
//     from_P2_to_P1.c_bit2[cur_depth][mul] = c1;


 

//     from_P0.next_bit_L_recv2[cur_depth][mul] = x0 ^ D0;
//    // std::cout << "sim: x0 ^ D0 = " << x0 << " ^ " << D0 << std::endl;
//     from_P0.next_bit_R_recv2[cur_depth][mul] = b0 ^ d0;

 

//     from_P1.next_bit_L_recv2[cur_depth][mul] = x1 ^ D1;
//    //  std::cout << "sim: x1 ^ D1 = " << x1 << " ^ " << D1 << std::endl;
//     from_P1.next_bit_R_recv2[cur_depth][mul] = b1 ^ d1;


//       const bool gamma0 = xor_if(c0, from_P1.next_bit_L_recv2[cur_depth][mul], d0);
  
//     bool xx0 = x0;   
    
    
//     const bool gamma1 = xor_if(c1, from_P0.next_bit_L_recv2[cur_depth][mul], d1);
   
//     bool xx1 = x1;  
   
//     results[0] = xor_if(gamma0, xx0, (b0 ^ from_P1.next_bit_R_recv2[cur_depth][mul]));
  
//     results[1] = xor_if(gamma1, xx1, (b1 ^ from_P0.next_bit_R_recv2[cur_depth][mul])); 

 
//   }



//     template<typename T>
//   void Simulator::multiply_mpc(const T& x0, bool b0, const T& x1, bool b1, T results[2], size_t mul, size_t cur_depth)
//   {   

//     T D0 = P0rand.next_block();
    
//     T D1 = P1rand.next_block();
  
    
//     bool d0 = P0rand.next_bool(); 
    
//     bool d1 = P1rand.next_bool(); 

   
//     const T alpha = P2rand.next_block();

//     T c0 = xor_if(alpha, D0, d1);
//     T c1 = xor_if(alpha, D1, d0);
 

//     from_P2_to_P0.c[cur_depth][mul] = c0;
//     from_P2_to_P1.c[cur_depth][mul] = c1;

 
//     from_P0.blinds_recv[cur_depth][mul]    = x0 ^ D0;
//     from_P0.bit_blinds_recv[cur_depth][mul] = b0 ^ d0;

 
//     from_P1.blinds_recv[cur_depth][mul]     = x1 ^ D1;    
//     from_P1.bit_blinds_recv[cur_depth][mul] = b1 ^ d1;


//     const T gamma0 = xor_if(c0, from_P1.blinds_recv[cur_depth][mul], d0);



//     T xx0 = x0;

   
//     const T gamma1 = xor_if(c1, from_P0.blinds_recv[cur_depth][mul], d1);
    
    
//     T xx1 = x1;

//     results[0] = xor_if(gamma0, xx0, (b0 ^ from_P1.bit_blinds_recv[cur_depth][mul]));
    

//     results[1] = xor_if(gamma1, xx1, (b1 ^ from_P0.bit_blinds_recv[cur_depth][mul])); 

 
//   }

// template<typename prgkey_t>
//  inline void Simulator::expand_mpc(prgkey_t & key, const block_t & seed0, const block_t & seed1,  block_t s0[2], bool t0[2], 
//                         block_t s1[2], bool t1[2], const block_t blind0[rounds], const block_t blind1[rounds], block_t gamma[2][rounds], size_t cur_depth, bool LR)
//  {
  
//     std::vector<block_t> c0(2 + rounds + rounds);
//     std::vector<block_t> c1(2 + rounds + rounds);
  

//    const block_t seed0L = clear_lsb(seed0); // _mm_clearlsb_si128(seed);
//    const block_t seed0R = set_lsb(seed0); //  _mm_setlsb_si128(seed);

//    const block_t seed1L = clear_lsb(seed1); // _mm_clearlsb_si128(seed);
//    const block_t seed1R = clear_lsb(seed1); //  _mm_setlsb_si128(seed);
 
//    s0[L] = seed0L;
//    s0[R] = seed0R;
 
//    s1[L] = seed1L;
//    s1[R] = seed1R;

   
//    std::pair<std::vector<block_t>,std::vector<block_t>> encrypt_outL = key.encrypt_MPC_proof(seed0L, seed1L, blind0, blind1, gamma);
      

//     if(!LR)
//     {


//      from_P1.seed0L_encrypt[cur_depth] = encrypt_outL.second;
//      //from_P1.seed0L_encrypt_other[cur_depth] = encrypt_outL.first;
 
     
//      from_P0.seed0L_encrypt[cur_depth] = encrypt_outL.first;
//      //from_P0.seed0L_encrypt_other[cur_depth] = encrypt_outL.second;
 
//     }

//     if(LR)
//     {

//      from_P1.seed1L_encrypt[cur_depth] = encrypt_outL.second;
//      //from_P1.seed1L_encrypt_other[cur_depth] = encrypt_outL.first;
 
     
//      from_P0.seed1L_encrypt[cur_depth] = encrypt_outL.first;
//      //from_P0.seed1L_encrypt_other[cur_depth] = encrypt_outL.second;
 
//     }

//     for(size_t j = 0; j < rounds; ++j)
//     {
//       c0[2 + 2*j]     = encrypt_outL.first[j+1];
//       c0[2 + 2*j + 1] = encrypt_outL.second[j+1];
//     }

//     std::pair<std::vector<block_t>,std::vector<block_t>> encrypt_outR = key.encrypt_MPC_proof(seed0R, seed1R, blind0, blind1, gamma);
    
//     if(!LR)
//     {

//     from_P1.seed0R_encrypt[cur_depth] = encrypt_outR.second;
//     //from_P1.seed0R_encrypt_other[cur_depth] = encrypt_outR.first;
 

//     from_P0.seed0R_encrypt[cur_depth] = encrypt_outR.first;
//     //from_P0.seed0R_encrypt_other[cur_depth] = encrypt_outR.second;
 
//     }

//     if(LR)
//     {
    
//     from_P1.seed1R_encrypt[cur_depth] = encrypt_outR.second; 
//     //from_P1.seed1R_encrypt_other[cur_depth] = encrypt_outR.first; 
 

//     from_P0.seed1R_encrypt[cur_depth] = encrypt_outR.first;
//     //from_P0.seed1R_encrypt_other[cur_depth] = encrypt_outR.second;
  
//     }

//     for(size_t j = 0; j < rounds; ++j)
//     {
//       c1[2 + 2*j]     = encrypt_outR.first[j+1];
//       c1[2 + 2*j + 1] = encrypt_outR.second[j+1];
//     }

//     s0[L] = encrypt_outL.first[rounds];
//     s0[R] = encrypt_outR.first[rounds];

//     s1[L] = encrypt_outL.second[rounds];
//     s1[R] = encrypt_outR.second[rounds];

//     s0[L] ^= seed0L; 
//     s1[L] ^= seed1L;
    
//     t0[L] = get_lsb(s0[L]); 
//     t1[L] = get_lsb(s1[L]);
     
//     s0[L] = clear_lsb(s0[L]); 
//     s1[L] = clear_lsb(s1[L]);

//     s0[R] ^= seed0R; s1[R] ^= seed1R;
//     t0[R] = get_lsb(s0[R]); t1[R] = get_lsb(s1[R]);
//     s0[R] = clear_lsb(s0[R]); s1[R] = clear_lsb(s1[R]);   

//     c0[L] = s0[L];
//     c1[L] = s1[L];

//     c0[R] = s0[R];
//     c1[R] = s1[R];


//     std::pair<std::pair<std::vector<block_t>,std::vector<block_t>>, std::pair<std::vector<block_t>,std::vector<block_t>>> encrypt_out = std::pair(encrypt_outL, encrypt_outR);
    
//     //return encrypt_out;
    
//     //std::make_pair(c0, c1);
  
//  }


//   void Simulator::conditional_swap_and_next_seed(LowMC& key, block_t s0_L[2], block_t s0_R[2], block_t s1_L[2], block_t s1_R[2], block_t left_out[2], block_t right_out[2], size_t cur_depth)
//   {
     
//      block_t bs0_L[2]; // shares of b \times s0_L
//      block_t notbs0_L[2]; // shares of (1 - b)  \times s0_L
//      block_t bs1_L[2]; // shares of b \times s1_L
//      block_t notbs1_L[2]; // shares of (1 - b) \times s1_L
//      block_t bs0_R[2]; // shares of b \times s0_R
//      block_t bs1_R[2]; // shares of b \times s1_R
//      block_t notbs1_R[2];  // shares of (1-b) \times s1_R
//      block_t notbs0_R[2]; // shares of (1-b) \times s0_R

//      multiply_mpc(s0_L[0], P0direction[cur_depth], s0_L[1], P1direction[cur_depth], bs0_L, 0, cur_depth); // L_0 * direction
//      notbs0_L[0] = s0_L[0] ^ bs0_L[0];
//      notbs0_L[1] = s0_L[1] ^ bs0_L[1];
     
//      multiply_mpc(s1_L[0], P0direction[cur_depth], s1_L[1], P1direction[cur_depth], bs1_L, 1, cur_depth); // L_1 * direction     
//      notbs1_L[0] = s1_L[0] ^ bs1_L[0];
//      notbs1_L[1] = s1_L[1] ^ bs1_L[1];
     
//      multiply_mpc(s0_R[0], P0direction[cur_depth], s0_R[1], P1direction[cur_depth], bs0_R, 2, cur_depth); // R_0 * direction      
//      notbs0_R[0] = s0_R[0] ^ bs0_R[0];
//      notbs0_R[1] = s0_R[1] ^ bs0_R[1];

//      multiply_mpc(s1_R[0], P0direction[cur_depth], s1_R[1], P1direction[cur_depth], bs1_R, 3, cur_depth); // R_1 * direction
//      notbs1_R[0] = s1_R[0] ^ bs1_R[0];
//      notbs1_R[1] = s1_R[1] ^ bs1_R[1];

//      right_out[0] = notbs0_L[0] ^ bs0_R[0] ^ notbs0_L[1] ^ bs0_R[1];
//      right_out[1] = notbs1_L[0] ^ bs1_R[0] ^ notbs1_L[1] ^ bs1_R[1];

//      left_out[0] = bs0_L[0] ^ notbs0_R[0] ^ bs1_L[0] ^ notbs1_R[0];
//      left_out[1] = bs0_L[1] ^ notbs0_R[1] ^ bs1_L[1] ^ notbs1_R[1];
    
//      seed0[cur_depth+1][0] = notbs0_L[0] ^ bs0_R[0];
//      seed0[cur_depth+1][1] = notbs0_L[1] ^ bs0_R[1];
//      seed1[cur_depth+1][0] = notbs1_L[0] ^ bs1_R[0];
//      seed1[cur_depth+1][1] = notbs1_L[1] ^ bs1_R[1];


 
//   }









//   void Simulator::middle_layers(LowMC& key, dpf_key<__mX, nitems> dpfkey[2], const size_t cur_depth) 
//   {
    
//     std::cout << "middle layer: " << cur_depth << std::endl;
//     block_t CW = dpfkey[1].cw[cur_depth];
//     block_t CW0 = dpfkey[0].cw[cur_depth];


//     block_t P0gamma[2][rounds];
//     block_t P0blinds0[rounds];
//     block_t P0blinds1[rounds];

//     Gen_Blinds(P0blinds0, P0blinds1, P0gamma);
    
//     for(size_t j = 0; j < rounds; ++j)
//     {
//       from_P2_to_P0.gamma0[cur_depth][j] = P0gamma[0][j];
//       from_P2_to_P1.gamma0[cur_depth][j] = P0gamma[1][j];
//     } 
    
//     block_t s0[2], s1[2];
//     bool  t0[2], t1[2]; 

//     expand_mpc(key, seed0[cur_depth][0], seed0[cur_depth][1],  s0, t0, s1, t1, P0blinds0, P0blinds1, P0gamma, cur_depth, false);
    
//     block_t s0_L[2];
//     block_t s0_R[2];
  
//     s0_L[0] = xor_if(s0[L], CW, bit0[cur_depth][0]);
//     s0_L[1] = xor_if(s1[L], CW, !bit0[cur_depth][1]);
//     s0_R[0] = xor_if(s0[R], CW, bit0[cur_depth][0]);
//     s0_R[1] = xor_if(s1[R], CW, !bit0[cur_depth][1]);


//     block_t ss0[2], ss1[2];

//     bool tt0[2], tt1[2]; 

//     block_t P1gamma[2][rounds];
//     block_t P1blinds0[rounds];
//     block_t P1blinds1[rounds];

//     Gen_Blinds(P1blinds0, P1blinds1, P1gamma);
    
//     for(size_t j = 0; j < rounds; ++j)
//     {
//     from_P2_to_P0.gamma1[cur_depth][j] = P1gamma[0][j];
//     from_P2_to_P1.gamma1[cur_depth][j] = P1gamma[1][j];
 
//     } 

   
//     expand_mpc(key, seed1[cur_depth][0], seed1[cur_depth][1],  ss0, tt0, ss1, tt1, P1blinds0, P1blinds1, P1gamma, cur_depth, true);

//     block_t s1_L[2];
//     block_t s1_R[2];

//     s1_L[0] = xor_if(ss0[L], CW, bit1[cur_depth][0]);
//     s1_L[1] = xor_if(ss1[L], CW, !bit1[cur_depth][1]);
//     s1_R[0] = xor_if(ss0[R], CW, bit1[cur_depth][0]);
//     s1_R[1] = xor_if(ss1[R], CW, !bit1[cur_depth][1]);

 
//     // from_P1.L_shares_recv[cur_depth] = s1_L[0];
//     // from_P1.R_shares_recv[cur_depth] = s1_R[0];
 
//     // from_P0.L_shares_recv[cur_depth] = s0_L[0];
//     // from_P0.R_shares_recv[cur_depth] = s0_R[0];

//     block_t left_out[2]; 
 
//     block_t right_out[2];
 
//     conditional_swap_and_next_seed(key, s0_L, s0_R, s1_L, s1_R, left_out, right_out, cur_depth);
    
    
//     bool bit0_L[2]; 

//     bool results[2];
    
//     multiply_mpc(dpfkey[0].t[cur_depth][L],  false, false, bit0[cur_depth][1], bit0_L, 0, cur_depth);
 
//     bit0_L[0] = bit0_L[0] ^ (t0[0]) ^ dpfkey[0].t[cur_depth][L] & (bit0[cur_depth][0]);
    
   
//     bit0_L[1] = bit0_L[1] ^ t1[0];

 
//     bool bit0_R[2]; 
//     multiply_mpc(dpfkey[0].t[cur_depth][R],  false, false, bit0[cur_depth][1], bit0_R, 1, cur_depth);

//     bit0_R[0] = bit0_R[0] ^ t0[1] ^ (dpfkey[0].t[cur_depth][R] & bit0[cur_depth][0]);
//     bit0_R[1] = bit0_R[1] ^ t1[1];

 
//     bool bit1_L[2];

 
//     multiply_mpc(bit1[cur_depth][0], false, false, dpfkey[1].t[cur_depth][L],  bit1_L, 2, cur_depth);
 
//     bit1_L[0] = bit1_L[0] ^ tt0[0];
 
//     bit1_L[1] = bit1_L[1] ^  tt1[0] ^ (dpfkey[1].t[cur_depth][L] & bit1[cur_depth][1]);


      
//     bool bit1_R[2];
     
//     multiply_mpc(bit1[cur_depth][0],  false, false, dpfkey[1].t[cur_depth][R],  bit1_R, 3, cur_depth);

//     bit1_R[0] = bit1_R[0] ^  tt0[1];
//     bit1_R[1] = bit1_R[1] ^ (tt1[1]) ^(dpfkey[1].t[cur_depth][R] & bit1[cur_depth][1]);

//     get_next_bits(bit0_L, bit0_R, bit1_L,  bit1_R, cur_depth);
    
//     std::cout << "reconstructed L = " << (left_out[0]  ^ left_out[1]) << std::endl;
     
// }

//   void Simulator::root_layer(LowMC& key, dpf_key<__mX, nitems> dpfkey[2])
//   {
//     std::cout << "root_layer: " << std::endl;
    
//     block_t seed0 = dpfkey[0].root; 
//     block_t seed1 = dpfkey[1].root;

//     size_t cur_depth = 0;
   
//     block_t CW = dpfkey[0].cw[cur_depth];
   
//     bool b0 = get_lsb(dpfkey[0].root);
    
//     bool b1 = get_lsb(dpfkey[1].root);

//     block_t s0[2], s1[2];    bool t0[2], t1[2];

//     block_t child[2];

//     expand_nonmpc(key, seed0, child, t0);
   
//     // s0[L] and s0[R] are P0's left and right children after the correction word is applied
//     s0[L] = xor_if(child[L], CW, !b0); 
//     s0[R] = xor_if(child[R], CW, !b0); 

//     expand_nonmpc(key, seed1, child, t1);
   
//     // s1[L] and s1[R] are P1's left and right children after the correction word is applied
//     s1[L] = xor_if(child[L], CW, !b1);
//     s1[R] = xor_if(child[R], CW, !b1);

//     block_t s0_L[2]; 
//     block_t s0_R[2];
   


//     s0_L[0] = P0rand.next_block();
//     s0_L[1] = s0_L[0] ^ s0[L];

//     s0_R[0] = P0rand.next_block();    
//     s0_R[1] = s0_R[0] ^ s0[R];
 
//     block_t s1_L[2]; 
//     block_t s1_R[2];
   
    
//     s1_L[0] = P1rand.next_block();
    
//     s1_L[1] = s1_L[0] ^ s1[L];
  
//     s1_R[0] = P1rand.next_block();
 
//     s1_R[1] = s1_R[0] ^ s1[R];

    
 
//     from_P1.L_shares_recv = s1_L[0];
//     from_P1.R_shares_recv = s1_R[0];

 
//     from_P0.L_shares_recv = s0_L[1];
//     from_P0.R_shares_recv = s0_R[1];


//     block_t left_out[2], right_out[2];
    
//     conditional_swap_and_next_seed(key, s0_L, s0_R, s1_L, s1_R, left_out, right_out, cur_depth);


//     std::cout << "reconstructed L = " << (left_out[0]  ^ left_out[1]) << std::endl;
   

//     // b0_L[2] are the shares of the cwbit for L-child of P0's root
//     bool bit0_L[2]; 

//     // b0_R[2] are the shares of the cwbit for R-child of P0's root
//     bool bit0_R[2]; 

//     bit0_L[0] = P0rand.next_bool();
//     bit0_L[1] = bit0_L[0] ^ (t0[L] ^ dpfkey[0].t[cur_depth][L] & (b0));
    
    

//     bit0_R[0] = P0rand.next_bool();
//     bit0_R[1] = bit0_R[0] ^ (t0[R] ^ dpfkey[0].t[cur_depth][R] & (b0));
    
//     // b1_L[2] are the shares of the cwbit for L-child of P1's root
//     bool bit1_L[2];
   
//     // b1_R[2] are the shares of the cwbit for R-child of P1's root
//     bool bit1_R[2];

//     bit1_L[0] = P1rand.next_bool();
//     bit1_L[1] = bit1_L[0] ^ (t1[L] ^ dpfkey[1].t[cur_depth][L] & (b1));
    
    

//     bit1_R[0] = P1rand.next_bool();
//     bit1_R[1] = bit1_R[0] ^ (t1[R] ^ dpfkey[1].t[cur_depth][R] & (b1));

 
//     from_P1.bit_L_shares_recv = bit1_L[0];
//     from_P1.bit_R_shares_recv = bit1_R[0];

 
//     from_P0.bit_L_shares_recv = bit0_L[1];
//     from_P0.bit_R_shares_recv = bit0_R[1];


//     get_next_bits(bit0_L, bit0_R, bit1_L, bit1_R, cur_depth);

//   }















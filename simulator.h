class  Simulator
{

  public:

  std::bitset<depth> P0direction;
  std::bitset<depth> P1direction;
  
  block seed0[depth][2];
  block seed1[depth][2];
  block sim_P0gamma[2][rounds];
  block sim_P0blinds0[rounds];
  block sim_P0blinds1[rounds];
  bool bit0_new[depth][2];
  bool bit1_new[depth][2];
  
  Simulator(AES_KEY& aeskey, __m128i seed0, __m128i seed1, __m128i seed2,  size_t len)
      : P0rand(aeskey, seed0, len), P1rand(aeskey, seed1, len), P2rand(aeskey, seed2, len)
  {


  }


  void Gen_Blinds(block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds]);

  
  void get_next_bits(bool L0_shares[2], bool R0_shares[2], bool L1_shares[2], bool R1_shares[2], size_t cur_depth);

  // void get_t_bits_for_next_layer(size_t cur_depth, size_t prev_dir, bool t0[2], bool t1[2], bool tt0[2], bool tt1[2]);
  template<typename KEY_TYPE, typename __mX>
  void expand_mpc(KEY_TYPE & key, const blocks<__mX> & seed0, const blocks<__mX> & seed1,  blocks<__mX> s0[2], bool t0[2], 
                        blocks<__mX> s1[2], bool t1[2], const block blind0[rounds], const block blind1[rounds], block gamma[2][rounds], size_t cur_depth, bool LR);
  
  void conditional_swap(LowMC& key, block L0_shares[2], block R0_shares[2], block L1_shares[2], block R1_shares[2], block left_out[2], block right_out[2], size_t cur_depth);

  void middle_layers(LowMC& key, dpf_key<__mX, nitems> dpfkey[2], size_t cur_depth);  

  void leaf_layers();

  template<typename T>
  void multiply_mpc(const T&, bool, const T&, bool, T results[2], size_t mul, size_t cur_depth);
  void multiply_mpc_b(const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth);
  void dpf_mpc(LowMC&, AES_KEY&, dpf_key<__mX, nitems>& k0, dpf_key<__mX, nitems> k1); 
  void root_layer(LowMC& key, dpf_key<__mX, nitems> dpfkey[2]);
 

  ~ Simulator()
  {

  }

//private:
  MPCrandomness P0rand, P1rand, P2rand;
  
};


 template<typename KEY_TYPE, typename __mX>
 inline void Simulator::expand_mpc(KEY_TYPE & key, const blocks<__mX> & seed0, const blocks<__mX> & seed1,  blocks<__mX> s0[2], bool t0[2], 
                        blocks<__mX> s1[2], bool t1[2], const block blind0[rounds], const block blind1[rounds], block gamma[2][rounds], size_t cur_depth, bool LR)
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
      
   if(!LR) P0_message.PB_middle.TL[cur_depth] = encrypt_outL.second;
   if(LR) P0_message.PB_middle.TL2[cur_depth] = encrypt_outL.second;
    
    
    if(!LR)
    {
     P0_message.PB_middle.TL_other[cur_depth] = encrypt_outL.first;
     P1_message.PB_middle.TL[cur_depth] = encrypt_outL.first;
    }

    if(LR)
    {
     P0_message.PB_middle.TL_other2[cur_depth] = encrypt_outL.first;
     P1_message.PB_middle.TL2[cur_depth] = encrypt_outL.first;
    }

    for(size_t j = 0; j < rounds; ++j)
    {
      c0[2 + 2*j]     = encrypt_outL.first[j+1];
      c0[2 + 2*j + 1] = encrypt_outL.second[j+1];
    }


    std::pair<std::vector<block>,std::vector<block>> encrypt_outR = key.encrypt_MPC_proof(seed0R, seed1R, blind0, blind1, gamma);
   

    if(!LR)
    {
    P0_message.PB_middle.TR[cur_depth] = encrypt_outR.second;
    P0_message.PB_middle.TR_other[cur_depth] = encrypt_outR.first;
    P1_message.PB_middle.TR[cur_depth] = encrypt_outR.first;
    }

    if(LR)
    {
    P0_message.PB_middle.TR2[cur_depth] = encrypt_outR.second;
    P0_message.PB_middle.TR_other2[cur_depth] = encrypt_outR.first;
    P1_message.PB_middle.TR2[cur_depth] = encrypt_outR.first;
    }

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


    std::pair<std::pair<std::vector<block>,std::vector<block>>, std::pair<std::vector<block>,std::vector<block>>> encrypt_out = std::pair(encrypt_outL, encrypt_outR);
    
    //return encrypt_out;
    
    //std::make_pair(c0, c1);
  
 }


  void Simulator::conditional_swap(LowMC& key, block L0_shares[2], block R0_shares[2], block L1_shares[2], block R1_shares[2], 
                                 block left_out[2], block right_out[2], size_t cur_depth)
  {

     block bL0Shares[2];
     
     multiply_mpc(L0_shares[0], P0direction[cur_depth], L0_shares[1], P1direction[cur_depth], bL0Shares,0, cur_depth); // L_0 * direction
     // if(cur_depth == 1) 
     // {  
     //     std::cout << "L0_shares[0] = " << L0_shares[0] << std::endl;
     //     std::cout << "L0_shares[1] = " << L0_shares[1] << std::endl;    
     //     std::cout << "bL0Shares[1] = " << bL0Shares[1] << std::endl;
     // }
     block bNotL0Shares[2];
     bNotL0Shares[0] = L0_shares[0] ^ bL0Shares[0];
     bNotL0Shares[1] = L0_shares[1] ^ bL0Shares[1];
     
     block bL1Shares[2];
     multiply_mpc(L1_shares[0], P0direction[cur_depth], L1_shares[1], P1direction[cur_depth], bL1Shares,1, cur_depth); // L_1 * direction

     block bNotL1Shares[2];
     bNotL1Shares[0] = L1_shares[0] ^ bL1Shares[0];
     bNotL1Shares[1] = L1_shares[1] ^ bL1Shares[1];
     
     
     block bR0Shares[2];
     multiply_mpc(R0_shares[0], P0direction[cur_depth], R0_shares[1], P1direction[cur_depth], bR0Shares,2, cur_depth); // R_0 * direction
 
     block bNotR0Shares[2];
     bNotR0Shares[0] = R0_shares[0] ^ bR0Shares[0];
     bNotR0Shares[1] = R0_shares[1] ^ bR0Shares[1];

    
     block bR1Shares[2];    
     multiply_mpc(R1_shares[0], P0direction[cur_depth], R1_shares[1], P1direction[cur_depth], bR1Shares, 3, cur_depth); // R_1 * direction
   
     block bNotR1Shares[2]; 
     bNotR1Shares[0] = R1_shares[0] ^ bR1Shares[0];
     bNotR1Shares[1] = R1_shares[1] ^ bR1Shares[1];

     right_out[0] = bNotL0Shares[0] ^ bR0Shares[0] ^ bNotL0Shares[1] ^ bR0Shares[1];
     right_out[1] = bNotL1Shares[0] ^ bR1Shares[0] ^ bNotL1Shares[1] ^ bR1Shares[1];

     left_out[0] = bL0Shares[0] ^ bNotR0Shares[0] ^ bL1Shares[0] ^ bNotR1Shares[0];
     left_out[1] = bL0Shares[1] ^ bNotR0Shares[1] ^ bL1Shares[1] ^ bNotR1Shares[1];
    
     seed0[cur_depth+1][0] = bNotL0Shares[0] ^ bR0Shares[0];
     seed0[cur_depth+1][1] = bNotL0Shares[1] ^ bR0Shares[1];
     seed1[cur_depth+1][0] = bNotL1Shares[0] ^ bR1Shares[0];
     seed1[cur_depth+1][1] = bNotL1Shares[1] ^ bR1Shares[1];


     P0_message.PB_middle.seedL[cur_depth+1] = seed0[cur_depth+1][0];
     P0_message.PB_middle.seedR[cur_depth+1] = seed1[cur_depth+1][0];
     P1_message.PB_middle.seedL[cur_depth+1] = seed0[cur_depth+1][1];
     P1_message.PB_middle.seedR[cur_depth+1] = seed1[cur_depth+1][1];
  }



  void Simulator::get_next_bits(bool L0_shares[2], bool R0_shares[2], bool L1_shares[2], bool R1_shares[2], size_t cur_depth)
  {
    
    bool bNotL0Shares[2];    
    multiply_mpc_b(L0_shares[0], P0direction[cur_depth], L0_shares[1], !P1direction[cur_depth], bNotL0Shares,  0, cur_depth); // L_0 * (direction - 1)

    bool bNotL1Shares[2];
    
    multiply_mpc_b(L1_shares[0], P0direction[cur_depth], L1_shares[1], !P1direction[cur_depth], bNotL1Shares,  1, cur_depth); // L_1 * (direction - 1)

    bool bR0Shares[2];
    
    multiply_mpc_b(R0_shares[0], P0direction[cur_depth], R0_shares[1], P1direction[cur_depth],  bR0Shares, 2  , cur_depth);

    bool bR1Shares[2];

    multiply_mpc_b(R1_shares[0], P0direction[cur_depth], R1_shares[1], P1direction[cur_depth], bR1Shares, 3  , cur_depth); // R_1 * direction    
    
    
    bit0_new[cur_depth+1][0] = bNotL0Shares[0] ^ bR0Shares[0];
    
    bit0_new[cur_depth+1][1] = bNotL0Shares[1] ^ bR0Shares[1];    
    
    bit1_new[cur_depth+1][0] = bNotL1Shares[0] ^ bR1Shares[0];    
    
    bit1_new[cur_depth+1][1] = bNotL1Shares[1] ^ bR1Shares[1];

    std::cout << "bit0_new[cur_depth+1][1] = " << bit0_new[cur_depth+1][1] << std::endl;
    std::cout << "bit1_new[cur_depth+1][1] = " << bit1_new[cur_depth+1][1] << std::endl;

  }

  void Simulator::Gen_Blinds(block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds])
  {
    block rand[rounds];
    for (unsigned r = 0; r < rounds; ++r)
     {
         block x;
         
          //P0rand.next_block(blinds0[r]);
         // P1rand.next(blinds1[r]); 
         // P2rand.next(rand[r]);

        blinds0[r] = P0rand.next_block();
        blinds1[r] = P1rand.next_block();
        rand[r]    = P2rand.next_block();
        const block tmp1 = ((blinds0[r] >> 1) & blinds1[r]) ^ ((blinds1[r] >> 1) & blinds0[r]);
        const block tmp2 = ((blinds0[r] >> 2) & blinds1[r]) ^ ((blinds1[r] >> 2) & blinds0[r]);
    
        const block bc = (tmp1 << 2) & maska;
        const block ac = (tmp2 << 1) & maskb;
        const block ab = (tmp1 >> 1) & maskc;
    
        gamma[0][r] = (bc | ac | ab) ^ rand[r];
        gamma[1][r] = rand[r];
     }
  }
  
  void Simulator::middle_layers(LowMC& key, dpf_key<__mX, nitems> dpfkey[2], const size_t cur_depth) {
    
    std::cout << "middle layer: " << cur_depth << std::endl;
    block CW = dpfkey[1].cw[cur_depth];
    block CW0 = dpfkey[0].cw[cur_depth];


    block P0gamma[2][rounds];
    block P0blinds0[rounds];
    block P0blinds1[rounds];

    Gen_Blinds(P0blinds0, P0blinds1, P0gamma);
    
    for(size_t j = 0; j < rounds; ++j)
    {
    P0_message.PB_middle.gamma0[cur_depth][j] = P0gamma[0][j];
    P1_message.PB_middle.gamma0[cur_depth][j] = P0gamma[1][j];
    } 
    
    block s0[2], s1[2];
    bool  t0[2], t1[2]; 



    expand_mpc(key, seed0[cur_depth][0], seed0[cur_depth][1],  s0, t0, s1, t1, P0blinds0, P0blinds1, P0gamma, cur_depth, false);

    std::cout << "seed0[cur_depth][1] = " << seed0[cur_depth][1] << std::endl; 

    block L0_shares[2]; 
    block R0_shares[2];
  
    L0_shares[0] = xor_if(s0[L], CW, bit0_new[cur_depth][0]);
    L0_shares[1] = xor_if(s1[L], CW, !bit0_new[cur_depth][1]);
    R0_shares[0] = xor_if(s0[R], CW, bit0_new[cur_depth][0]);
    R0_shares[1] = xor_if(s1[R], CW, !bit0_new[cur_depth][1]);
     
    block ss0[2], ss1[2];

    bool tt0[2], tt1[2]; 

    block P1gamma[2][rounds];
    block P1blinds0[rounds];
    block P1blinds1[rounds];

    Gen_Blinds(P1blinds0, P1blinds1, P1gamma);
    
    for(size_t j = 0; j < rounds; ++j)
    {
    P0_message.PB_middle.gamma1[cur_depth][j] = P1gamma[0][j];
    P1_message.PB_middle.gamma1[cur_depth][j] = P1gamma[1][j];
    } 

   
    expand_mpc(key, seed1[cur_depth][0], seed1[cur_depth][1],  ss0, tt0, ss1, tt1, P1blinds0, P1blinds1, P1gamma, cur_depth, true);


    block L1_shares[2];
    block R1_shares[2];

    L1_shares[0] = xor_if(ss0[L], CW, bit1_new[cur_depth][0]);
    L1_shares[1] = xor_if(ss1[L], CW, !bit1_new[cur_depth][1]);
    R1_shares[0] = xor_if(ss0[R], CW, bit1_new[cur_depth][0]);
    R1_shares[1] = xor_if(ss1[R], CW, !bit1_new[cur_depth][1]);


    // std::cout << "ss0[L] = " << ss0[L] << std::endl;
    // std::cout << "ss0[R] = " << ss0[R] << std::endl;
    std::cout << "LHS: " << (L0_shares[0] ^ L1_shares[0] ^ L0_shares[1] ^ L1_shares[1]) << std::endl;
    std::cout << "RHS: " << (R0_shares[0] ^ R1_shares[0] ^ R0_shares[1] ^ R1_shares[1]) << std::endl;
    block left_out[2]; block right_out[2] ;
     
    P0_message.PB_middle.L_shares_sent[cur_depth] = L0_shares[0];
    P0_message.PB_middle.R_shares_sent[cur_depth] = R0_shares[0];
    P0_message.PB_middle.L_shares_recv[cur_depth] = L1_shares[0];
    P0_message.PB_middle.R_shares_recv[cur_depth] = R1_shares[0];

    P1_message.PB_middle.L_shares_sent[cur_depth] = L1_shares[0];
    P1_message.PB_middle.R_shares_sent[cur_depth]= R1_shares[0];
    P1_message.PB_middle.L_shares_recv[cur_depth] = L0_shares[0];
    P1_message.PB_middle.R_shares_recv[cur_depth] = R0_shares[0];



    conditional_swap(key, L0_shares, R0_shares, L1_shares, R1_shares, left_out, right_out, cur_depth);
    
    
     bool bitL0_shares[2]; 
     bitL0_shares[0] = t1[0];
     bitL0_shares[1] = t0[0] ^ dpfkey[0].t[cur_depth][L] & (bit0_new[cur_depth][0] ^ bit0_new[cur_depth][1]);
    
     bool bitR0_shares[2]; 
     bitR0_shares[0] = t1[1]; 
     bitR0_shares[1] = t0[1] ^ dpfkey[0].t[cur_depth][R] & (bit0_new[cur_depth][0] ^ bit0_new[cur_depth][1]);
    
     bool bitL1_shares[2];
     bitL1_shares[0] = tt1[0];
     bitL1_shares[1] = tt0[0] ^ dpfkey[1].t[cur_depth][L] & (bit1_new[cur_depth][0] ^ bit1_new[cur_depth][1]);
    
     bool bitR1_shares[2];
     bitR1_shares[0] = tt0[1];
     bitR1_shares[1] = tt1[1] ^ dpfkey[1].t[cur_depth][R] & (bit1_new[cur_depth][0] ^ bit1_new[cur_depth][1]);
     
    P0_message.PB_middle.bit_L_shares_sent[cur_depth] = bitL0_shares[0];
    P0_message.PB_middle.bit_R_shares_sent[cur_depth] = bitR0_shares[0];
    P0_message.PB_middle.bit_L_shares_recv[cur_depth] = bitL1_shares[0];
    P0_message.PB_middle.bit_R_shares_recv[cur_depth] = bitR1_shares[0];

    P1_message.PB_middle.bit_L_shares_sent[cur_depth] = bitL1_shares[0];
    P1_message.PB_middle.bit_R_shares_sent[cur_depth] = bitR1_shares[0];
    P1_message.PB_middle.bit_L_shares_recv[cur_depth] = bitL0_shares[0];
    P1_message.PB_middle.bit_R_shares_recv[cur_depth] = bitR0_shares[0];
   
    get_next_bits(bitL0_shares, bitR0_shares, bitL1_shares,  bitR1_shares, cur_depth);
    
    std::cout << "reconstructed L = " << (left_out[0]  ^ left_out[1]) << std::endl;
    std::cout << "reconstructed R = " << (right_out[0] ^ right_out[1]) << std::endl << std::endl;

    
}

  void Simulator::root_layer(LowMC& key, dpf_key<__mX, nitems> dpfkey[2])
  {
    std::cout << "root_layer: " << std::endl;
    block seed0 = dpfkey[0].root; 
    block seed1 = dpfkey[1].root;

    size_t cur_depth = 0;
   
    block CW = dpfkey[0].cw[cur_depth];
   
    bool b0 = get_lsb(dpfkey[0].root);
    
    bool b1 = get_lsb(dpfkey[1].root);

    block s0[2], s1[2];    bool t0[2], t1[2];

    block child[2];

    expand_nonmpc(key, seed0, child, t0);
   
    s0[L] = xor_if(child[L], CW, !b0);
    s0[R] = xor_if(child[R], CW, !b0);

    expand_nonmpc(key, seed1, child, t1);
   
    s1[L] = xor_if(child[L], CW, !b1);
    s1[R] = xor_if(child[R], CW, !b1);

    block L0_shares[2]; block R0_shares[2];
   
    P0rand.next(L0_shares[0]); 
    
    L0_shares[1] = L0_shares[0] ^ s0[L];P0rand.next(R0_shares[0]);
    R0_shares[1] = R0_shares[0] ^ s0[R];
 
    block L1_shares[2]; block R1_shares[2];
   
    P1rand.next(L1_shares[0]);
    L1_shares[1] = L1_shares[0] ^ s1[L];
    P1rand.next(R1_shares[0]); 
    R1_shares[1] = R1_shares[0] ^ s1[R];

    
    // P0_message.PB_root.L_shares_sent = L0_shares[1];
    // P0_message.PB_root.R_shares_sent = R0_shares[1];
    P0_message.PB_root.L_shares_recv = L1_shares[0];
    P0_message.PB_root.R_shares_recv = R1_shares[0];

    // P1_message.PB_root.L_shares_sent = L1_shares[0];
    // P1_message.PB_root.R_shares_sent = R1_shares[0];
    P1_message.PB_root.L_shares_recv = L0_shares[1];
    P1_message.PB_root.R_shares_recv = R0_shares[1];

    block left_out[2]; block right_out[2];
    
    conditional_swap(key, L0_shares, R0_shares, L1_shares, R1_shares, left_out, right_out, cur_depth);

    std::cout << "reconstructed L = " << (left_out[0]  ^ left_out[1]) << std::endl;
   
    bool bitL0_shares[2]; 
    int dd; 
    bitL0_shares[0] = P0rand.next(dd); //rand();
    bitL0_shares[1] =  bitL0_shares[0] ^ (t0[L] ^ dpfkey[0].t[cur_depth][L] & (b0));
    
    bool bitR0_shares[2]; 
    bitR0_shares[0] = P0rand.next(dd); //rand();
    bitR0_shares[1] = bitR0_shares[0] ^ (t0[R] ^ dpfkey[0].t[cur_depth][R] & (b0));
    
    bool bitL1_shares[2];
    bitL1_shares[0] = P1rand.next(dd); //rand();
    bitL1_shares[1] = bitL1_shares[0] ^ (t1[L] ^ dpfkey[1].t[cur_depth][L] & (b1));
    
    bool bitR1_shares[2];
    bitR1_shares[0] = P1rand.next(dd); //rand();
    bitR1_shares[1] = bitR1_shares[0] ^ (t1[R] ^ dpfkey[1].t[cur_depth][R] & (b1));




    // P0_message.PB_root.bit_L_shares_sent = bitL0_shares[0];
    // P0_message.PB_root.bit_R_shares_sent = bitR0_shares[0];
    P0_message.PB_root.bit_L_shares_recv = bitL1_shares[0];
    P0_message.PB_root.bit_R_shares_recv = bitR1_shares[0];

    // P1_message.PB_root.bit_L_shares_sent = bitL1_shares[0];
    // P1_message.PB_root.bit_R_shares_sent = bitR1_shares[0];
    P1_message.PB_root.bit_L_shares_recv = bitL0_shares[1];
    P1_message.PB_root.bit_R_shares_recv = bitR0_shares[1];


    get_next_bits(bitL0_shares, bitR0_shares, bitL1_shares,  bitR1_shares, cur_depth);

  }

  void Simulator::multiply_mpc_b(const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth)
  {
    int next_b;

    bool D0 = P0rand.next(next_b); 
    bool D1 = P1rand.next(next_b); 
     
    bool d0 = P0rand.next(next_b);     
    bool d1 = P1rand.next(next_b);
      
    const bool alpha = P2rand.next(next_b); 

    bool c0 = xor_if(alpha, D0, d1);
    bool c1 = xor_if(alpha, D1, d0);

    P0_message.PB_middle.next_bit_L_sent[cur_depth][mul] = x0 ^ D0;
    P0_message.PB_middle.next_bit_R_sent[cur_depth][mul] = b0 ^ d0;

    P1_message.PB_middle.next_bit_L_sent[cur_depth][mul] = x1 ^ D1;
    P1_message.PB_middle.next_bit_R_sent[cur_depth][mul] = b1 ^ d1;


    const bool gamma0 = xor_if(c0, P1_message.PB_middle.next_bit_L_sent[cur_depth][mul], d0);
  
    bool xx0 = x0;   
    xor_if(gamma0, xx0, (b0 ^ P1_message.PB_middle.next_bit_R_sent[cur_depth][mul])); 

    const bool gamma1 = xor_if(c1,P0_message.PB_middle.next_bit_L_sent[cur_depth][mul], d1);
    bool xx1 = x1;  
    xor_if(gamma1, xx1, (b1 ^ P0_message.PB_middle.next_bit_R_sent[cur_depth][mul])); 

    results[0] = xor_if(gamma0, xx0, (b0 ^ P1_message.PB_middle.next_bit_R_sent[cur_depth][mul]));
  
   results[1] = xor_if(gamma1, xx1, (b1 ^ P0_message.PB_middle.next_bit_R_sent[cur_depth][mul])); 

   if(cur_depth == 0)
    {
      // P0_message.PB_root.next_bit_L_sent[mul] = P0_message.PB_middle.next_bit_L_sent[cur_depth][mul];
      // P1_message.PB_root.next_bit_L_recv[mul] = P0_message.PB_root.next_bit_L_sent[mul];
      // P0_message.PB_root.next_bit_R_sent[mul] = P1_message.PB_middle.bit_blinds_recv[cur_depth][mul];
      // P1_message.PB_root.next_bit_R_recv[mul] = P0_message.PB_root.next_bit_R_sent[mul];
      
      P0_message.PB_root.gamma_bit[mul] = gamma0; 
      P1_message.PB_root.gamma_bit[mul] = gamma1; 

      P1_message.PB_root.next_bit_L_sent[mul] = P0_message.PB_middle.next_bit_L_sent[cur_depth][mul];
      P0_message.PB_root.next_bit_L_recv[mul] = P1_message.PB_root.next_bit_L_sent[mul];
      P1_message.PB_root.next_bit_L_recv[mul] = P0_message.PB_middle.next_bit_L_sent[cur_depth][mul];
      P1_message.PB_root.next_bit_R_sent[mul] = P0_message.PB_middle.next_bit_R_sent[cur_depth][mul];
      P0_message.PB_root.next_bit_R_recv[mul] = P1_message.PB_root.next_bit_R_sent[mul];
      P1_message.PB_root.next_bit_R_recv[mul] = P0_message.PB_middle.next_bit_R_sent[cur_depth][mul];
    }
    else
    {
      P1_message.PB_middle.next_bit_L_recv[cur_depth][mul] = P0_message.PB_middle.next_bit_L_sent[cur_depth][mul];
      P1_message.PB_middle.next_bit_R_recv[cur_depth][mul] = P0_message.PB_middle.next_bit_R_sent[cur_depth][mul];
      P0_message.PB_middle.next_bit_L_recv[cur_depth][mul] = P1_message.PB_middle.next_bit_L_sent[cur_depth][mul];
      P0_message.PB_middle.next_bit_R_recv[cur_depth][mul] = P1_message.PB_middle.next_bit_R_sent[cur_depth][mul];
      P0_message.PB_middle.gamma_bit[cur_depth][mul]       = gamma0;
      P1_message.PB_middle.gamma_bit[cur_depth][mul]       = gamma1;
    }

  }






  template<typename T>
  void Simulator::multiply_mpc(const T& x0, bool b0, const T& x1, bool b1, T results[2], size_t mul, size_t cur_depth)
  {
   
    T D0; P0rand.next(D0);
    
    T D1; P1rand.next(D1);
  //  std::cout << "D1 = " << D1 << std::endl;
    int dd0; 
    
    bool d0 = P0rand.next(dd0);

    int dd1; 
    
    bool d1 = P1rand.next(dd1); 

    T rand;
   
    P2rand.next(rand); 

    const T alpha = rand;

    T c0 = xor_if(alpha, D0, d1);
    T c1 = xor_if(alpha, D1, d0);

    //P0_message.PB_middle.blinds_sent[cur_depth][mul]    = x0 ^ D0;
    
    P1_message.PB_middle.blinds_recv[cur_depth][mul]    = x0 ^ D0;

//    P0_message.PB_middle.bit_blinds_sent[cur_depth][mul] 

    P1_message.PB_middle.bit_blinds_recv[cur_depth][mul] = b0 ^ d0;

    //P1_message.PB_middle.blinds_sent[cur_depth][mul] = x1 ^ D1;
    
    P0_message.PB_middle.blinds_recv[cur_depth][mul]= x1 ^ D1;
    
    P0_message.PB_middle.bit_blinds_recv[cur_depth][mul] = b1 ^ d1;



    const T gamma0 = xor_if(c0, P0_message.PB_middle.blinds_recv[cur_depth][mul], d0);



    T xx0 = x0;

   
    const T gamma1 = xor_if(c1, P1_message.PB_middle.blinds_recv[cur_depth][mul], d1);
    
    
    T xx1 = x1;

    results[0] = xor_if(gamma0, xx0, (b0 ^ P0_message.PB_middle.bit_blinds_recv[cur_depth][mul]));
    
    // if(mul == 0)
    // {
    //     std::cout << "gamma1 =   " << gamma1 << std::endl;
    //     std::cout << "xx1 =  " << xx1 << std::endl;
    //     std::cout << "b1 =  " << b1 << std::endl;
    //     std::cout << "Y1 =    " << P1_message.PB_middle.bit_blinds_recv[cur_depth][mul] << std::endl;
    // }
    results[1] = xor_if(gamma1, xx1, (b1 ^ P1_message.PB_middle.bit_blinds_recv[cur_depth][mul])); 

   if(cur_depth == 0)
    {
      //P0_message.PB_root.blinds_sent[mul] = P1_message.PB_middle.blinds_recv[cur_depth][mul]; //P0_message.PB_middle.blinds_sent[cur_depth][mul];
      P1_message.PB_root.blinds_recv[mul] = P1_message.PB_middle.blinds_recv[cur_depth][mul];//P0_message.PB_root.blinds_sent[mul];
   //  P0_message.PB_root.bit_blinds_sent[mul] = P1_message.PB_middle.bit_blinds_recv[cur_depth][mul];
      P1_message.PB_root.bit_blinds_recv[mul] = P1_message.PB_middle.bit_blinds_recv[cur_depth][mul];// P0_message.PB_root.bit_blinds_sent[mul];
      


  //    P1_message.PB_root.blinds_sent[mul] = P0_message.PB_middle.blinds_recv[cur_depth][mul]; //P1_message.PB_middle.blinds_sent[cur_depth][mul];
      P0_message.PB_root.blinds_recv[mul] = P0_message.PB_middle.blinds_recv[cur_depth][mul];//P1_message.PB_root.blinds_sent[mul];
     // P1_message.PB_root.bit_blinds_sent[mul] = P0_message.PB_middle.bit_blinds_recv[cur_depth][mul];
       P0_message.PB_root.bit_blinds_recv[mul] = P0_message.PB_middle.bit_blinds_recv[cur_depth][mul] ;//P1_message.PB_root.bit_blinds_sent[mul];

      P2_message[0].P2_root.blinds0[mul] = D0;
      P2_message[1].P2_root.blinds1[mul] = D1;
      P2_message[0].P2_root.bit_blinds0[mul] = d0;
      P2_message[1].P2_root.bit_blinds1[mul] = d1;

      P2_message[0].P2_root.gamma[mul] = gamma0;
      P0_message.PB_root.gamma[mul]    = gamma0; 
      P1_message.PB_root.gamma[mul]    = gamma1;
      P2_message[1].P2_root.gamma[mul] = gamma1;
    }
    else
    {

      //P0_message.PB_middle.blinds_sent[cur_depth][mul] = X0;
      //P1_message.PB_middle.blinds_recv[cur_depth][mul] = P0_message.PB_middle.blinds_sent[cur_depth][mul];
      //P0_message.PB_middle.bit_blinds_sent[cur_depth][mul] = Y0;

     // P1_message.PB_middle.bit_blinds_recv[cur_depth][mul] = P0_message.PB_middle.bit_blinds_sent[cur_depth][mul];
      
      //P1_message.PB_middle.blinds_sent[cur_depth][mul] = X1;
     // P0_message.PB_middle.blinds_recv[cur_depth][mul] = P1_message.PB_middle.blinds_sent[cur_depth][mul];
      //P1_message.PB_middle.bit_blinds_sent[cur_depth][mul] = Y1;

      //P0_message.PB_middle.bit_blinds_recv[cur_depth][mul] = P1_message.PB_middle.bit_blinds_sent[cur_depth][mul];

      P0_message.PB_middle.gamma_middle[cur_depth][mul]    = gamma0; 

      P1_message.PB_middle.gamma_middle[cur_depth][mul]    = gamma1;
     
      P2_message[0].P2_middle.blinds0[cur_depth][mul] = D0;
      P2_message[1].P2_middle.blinds1[cur_depth][mul] = D1;
      P2_message[0].P2_middle.bit_blinds0[cur_depth][mul] = d0;
      P2_message[1].P2_middle.bit_blinds1[cur_depth][mul] = d1;

      P2_message[0].P2_middle.gamma[cur_depth][mul] = gamma0;
      P2_message[1].P2_middle.gamma[cur_depth][mul] = gamma1;
    }
  }


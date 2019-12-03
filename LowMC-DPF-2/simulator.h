class  Simulator
{

  public:

  std::bitset<depth> P0direction;
  std::bitset<depth> P1direction;
  
  block seed0[depth][2];
  block seed1[depth][2];
  
  bool bit0_new[depth][2];
  bool bit1_new[depth][2];
  
  Simulator(AES_KEY& aeskey, __m128i seed0, __m128i seed1, __m128i seed2,  size_t len)
      : P0rand(aeskey, seed0, len), P1rand(aeskey, seed1, len), P2rand(aeskey, seed2, len)
  {


  }


  void Gen_Blinds(block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds]);

  
  void get_next_bits(bool L0_shares[2], bool R0_shares[2], bool L1_shares[2], bool R1_shares[2], size_t cur_depth);

  // void get_t_bits_for_next_layer(size_t cur_depth, size_t prev_dir, bool t0[2], bool t1[2], bool tt0[2], bool tt1[2]);
   
  void conditional_swap(LowMC& key, block L0_shares[2], block R0_shares[2], block L1_shares[2], block R1_shares[2], block left_out[2], block right_out[2], size_t cur_depth);

  void middle_layers(LowMC& key, dpf_key<__mX, nitems> dpfkey[2], size_t cur_depth);  

  template<typename T>
  void multiply_mpc(const T&, bool, const T&, bool, T results[2], size_t mul, size_t cur_depth);
  void multiply_mpc_b(const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth);
  void dpf_mpc(LowMC&, AES_KEY&, dpf_key<__mX, nitems>& k0, dpf_key<__mX, nitems> k1); 
  void root_layer(LowMC& key, dpf_key<__mX, nitems> dpfkey[2]);
 

  ~ Simulator()
  {

  }

private:
  MPCrandomness P0rand, P1rand, P2rand;
  
};


  void Simulator::conditional_swap(LowMC& key, block L0_shares[2], block R0_shares[2], block L1_shares[2], block R1_shares[2], 
                                 block left_out[2], block right_out[2], size_t cur_depth)
  {
    
     std::cout << "direction = " << P0direction[cur_depth] << " ^ " << P1direction[cur_depth]
                                 << " = " << (P0direction[cur_depth] ^ P1direction[cur_depth]) 
                                 << std::endl; 
     block bL0Shares[2];
     
     multiply_mpc(L0_shares[0], P0direction[cur_depth], L0_shares[1], P1direction[cur_depth], bL0Shares,0, cur_depth); // L_0 * direction
    
    if(cur_depth == 0) std::cout << "result[" << 0 << "] = " << bL0Shares[0] << std::endl << std::endl;

     block bNotL0Shares[2];
     bNotL0Shares[0] = L0_shares[0] ^ bL0Shares[0];
     bNotL0Shares[1] = L0_shares[1] ^ bL0Shares[1];
     
     if(cur_depth == 0) {
      std::cout << "L0_shares[0]    = " << L0_shares[0] << std::endl;
      std::cout << "bNotL0Shares[0] = " << bNotL0Shares[0] << std::endl; 
     }
     block bL1Shares[2];
     multiply_mpc(L1_shares[0], P0direction[cur_depth], L1_shares[1], P1direction[cur_depth], bL1Shares,1, cur_depth); // L_1 * direction
     if(cur_depth == 0) std::cout << "result[" << 1 << "] = " << bL1Shares[0] << std::endl << std::endl;

     block bNotL1Shares[2];
     bNotL1Shares[0] = L1_shares[0] ^ bL1Shares[0];
     bNotL1Shares[1] = L1_shares[1] ^ bL1Shares[1];
     
     block bR0Shares[2];
     multiply_mpc(R0_shares[0], P0direction[cur_depth], R0_shares[1], P1direction[cur_depth], bR0Shares,2, cur_depth); // R_0 * direction
     if(cur_depth == 0) std::cout << "result[" << 2 << "] = " << bR0Shares[0] << std::endl << std::endl;

     block bNotR0Shares[2];
     bNotR0Shares[0] = R0_shares[0] ^ bR0Shares[0];
     bNotR0Shares[1] = R0_shares[1] ^ bR0Shares[1];

    
     block bR1Shares[2];    
     multiply_mpc(R1_shares[0], P0direction[cur_depth], R1_shares[1], P1direction[cur_depth], bR1Shares, 3, cur_depth); // R_1 * direction
     if(cur_depth == 0) std::cout << "result[" << 3 << "] = " << bR1Shares[0] << std::endl << std::endl;

   
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

     if(cur_depth == 0)
     {
      std::cout << "[SIM] seed0[cur_depth+1][0] =  " << seed0[cur_depth+1][0] << std::endl;
      std::cout << "[SIM] seed1[cur_depth+1][0] =  " << seed1[cur_depth+1][0] << std::endl;
     }
  }



  void Simulator::get_next_bits(bool L0_shares[2], bool R0_shares[2], bool L1_shares[2], bool R1_shares[2], size_t cur_depth)
  {
    
    bool bNotL0Shares[2];
     std::cout << "L0Shares: " <<  (L0_shares[0] ^ L0_shares[1]) << std::endl;
    
    multiply_mpc_b(L0_shares[0], P0direction[cur_depth], L0_shares[1], !P1direction[cur_depth], bNotL0Shares, 30, 0); // L_0 * (direction - 1)

    std::cout << "bNotL0Shares: " << (bNotL0Shares[0]) << " ^  " << (bNotL0Shares[1]) << "  " << (L0_shares[0] ^ L0_shares[1]) << std::endl;

    bool bNotL1Shares[2];
    std::cout << "L1Shares: " <<  (L1_shares[0] ^ L1_shares[1]) << std::endl;
    
    multiply_mpc_b(L1_shares[0], P0direction[cur_depth], L1_shares[1], !P1direction[cur_depth], bNotL1Shares, 30, 0); // L_1 * (direction - 1)
    std::cout << "bNotL1Shares: " <<  (bNotL1Shares[0] ^ bNotL1Shares[1]) << "  " << (L1_shares[0] ^ L1_shares[1]) << std::endl;

    bool bR0Shares[2];
    
    multiply_mpc_b(R0_shares[0], P0direction[cur_depth], R0_shares[1], P1direction[cur_depth],  bR0Shares, 30 , 0);
    
    std::cout << "bR0Shares: " <<  (bR0Shares[0] ^ bR0Shares[1]) << "  " << (R0_shares[0] ^ R0_shares[1]) << std::endl;
    
    bool bR1Shares[2];

    multiply_mpc_b(R1_shares[0], P0direction[cur_depth], R1_shares[1], P1direction[cur_depth], bR1Shares, 30 , 0); // R_1 * direction
    
    std::cout << "bR1Shares: " << (bR1Shares[0] ^ bR1Shares[1]) << "  " << (R1_shares[0] ^ R1_shares[1]) << std::endl;
 
    std::cout << bNotL0Shares[0] << " ^  " << bR0Shares[0] << std::endl;
    
    bit0_new[cur_depth+1][0] = bNotL0Shares[0] ^ bR0Shares[0];
    std::cout << " bit0_new[cur_depth+1][0] " << bit0_new[cur_depth+1][0] << std::endl;
    
    bit0_new[cur_depth+1][1] = bNotL0Shares[1] ^ bR0Shares[1];    
    std::cout << " bit0_new[cur_depth+1][1] " << bit0_new[cur_depth+1][1] << std::endl;
    
    bit1_new[cur_depth+1][0] = bNotL1Shares[0] ^ bR1Shares[0];    
    std::cout << " bit1_new[cur_depth+1][0] " << bit1_new[cur_depth+1][0] << std::endl;
    
    bit1_new[cur_depth+1][1] = bNotL1Shares[1] ^ bR1Shares[1];
    std::cout << " bit1_new[cur_depth+1][1] " << bit1_new[cur_depth+1][1] << std::endl;
  }

  void Simulator::Gen_Blinds(block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds])
  {
    block rand[rounds];
    for (unsigned r = 0; r < rounds; ++r)
     {
        P0rand.next(blinds0[r]);
        P1rand.next(blinds1[r]); 
        P2rand.next(rand[r]);

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

    Gen_Blinds(P0blinds0 ,P0blinds1, P0gamma);
  

    block s0[2], s1[2];
    bool t0[2], t1[2]; 
   
    std::pair<std::vector<block>,std::vector<block>> expand_transcript0 = expand_mpc(key, seed0[cur_depth][0], seed0[cur_depth][1],  s0, t0, s1, t1, P0blinds0, P0blinds1, P0gamma);

     
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

    Gen_Blinds(P1blinds0 ,P1blinds1, P1gamma);
    std::pair<std::vector<block>,std::vector<block>> expand_transcript1 = expand_mpc(key, seed1[cur_depth][0], seed1[cur_depth][1],  ss0, tt0, ss1, tt1, P1blinds0, P1blinds1, P1gamma);

    block L1_shares[2];
    block R1_shares[2];

    L1_shares[0] = xor_if(ss0[L], CW, bit1_new[cur_depth][0]);
    L1_shares[1] = xor_if(ss1[L], CW, !bit1_new[cur_depth][1]);
    R1_shares[0] = xor_if(ss0[R], CW, bit1_new[cur_depth][0]);
    R1_shares[1] = xor_if(ss1[R], CW, !bit1_new[cur_depth][1]);

    std::cout << "LHS: " << (L0_shares[0] ^ L1_shares[0] ^ L0_shares[1] ^ L1_shares[1]) << std::endl;
    std::cout << "RHS: " << (R0_shares[0] ^ R1_shares[0] ^ R0_shares[1] ^ R1_shares[1]) << std::endl;
     block left_out[2]; block right_out[2] ;
     
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

    
    std::cout << "seed0 = " << seed0 << std::endl;
    std::cout << "seed1 = " << seed1 << std::endl;  

    block child[2];

    expand_nonmpc(key, seed0, child, t0);
   
    s0[L] = xor_if(child[L], CW, !b0);
    s0[R] = xor_if(child[R], CW, !b0);

    std::cout << "s0[" << L << "] = " << s0[L] << std::endl;
    std::cout << "s0[" << R << "] = " << s0[R] << std::endl;
    expand_nonmpc(key, seed1, child, t1);
   
    s1[L] = xor_if(child[L], CW, !b1);
    s1[R] = xor_if(child[R], CW, !b1);

    block L0_shares[2]; block R0_shares[2];
   
    P0rand.next(L0_shares[0]); 
    L0_shares[1] = L0_shares[0] ^ s0[L];
    P0rand.next(R0_shares[0]);
    R0_shares[1] = R0_shares[0] ^ s0[R];
 

    std::cout << "L0_shares = " << L0_shares[0] << " " << L0_shares[1] << std::endl;
    std::cout << "R0_shares = " << R0_shares[0] << " " << R0_shares[1] << std::endl;
    block L1_shares[2]; block R1_shares[2];
   
    P1rand.next(L1_shares[0]);
    L1_shares[1] = L1_shares[0] ^ s1[L];
    P1rand.next(R1_shares[0]); 
    R1_shares[1] = R1_shares[0] ^ s1[R];

    P0_message.PB_root.L_shares_sent = L0_shares[1];
    P0_message.PB_root.R_shares_sent = R0_shares[1];
    P0_message.PB_root.L_shares_recv = L1_shares[0];
    P0_message.PB_root.R_shares_recv = R1_shares[0];

    P1_message.PB_root.L_shares_sent = L1_shares[0];
    P1_message.PB_root.R_shares_sent = R1_shares[0];
    P1_message.PB_root.L_shares_recv = L0_shares[1];
    P1_message.PB_root.R_shares_recv = R0_shares[1];

    block left_out[2]; block right_out[2];
    
    conditional_swap(key, L0_shares, R0_shares, L1_shares, R1_shares, left_out, right_out, cur_depth);

    std::cout << "reconstructed L = " << (left_out[0]  ^ left_out[1]) << std::endl;
   
    bool bitL0_shares[2]; 
    bitL0_shares[0] = t0[L];
    bitL0_shares[1] = dpfkey[0].t[cur_depth][L] & (b0);
    
    bool bitR0_shares[2]; 
    bitR0_shares[0] = t0[R];
    bitR0_shares[1] = dpfkey[0].t[cur_depth][R] & (b0);
    
    bool bitL1_shares[2];
    bitL1_shares[0] = t1[L];
  
    bitL1_shares[1] = dpfkey[1].t[cur_depth][L] & (b1);
    
    bool bitR1_shares[2];

    bitR1_shares[0] = t1[R];
    bitR1_shares[1] = dpfkey[1].t[cur_depth][R] & (b1);

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

    bool X0 = x0 ^ D0;
    bool Y0 = b0 ^ d0;

    bool X1 = x1 ^ D1;
    bool Y1 = b1 ^ d1;

    const bool gamma0 = xor_if(c0, X1, d0);
    bool xx0 = x0;   
    xor_if(gamma0, xx0, (b0 ^ Y1)); 

    const bool gamma1 = xor_if(c1, X0, d1);
    bool xx1 = x1;  
    xor_if(gamma1, xx1, (b1 ^ Y0)); 

    results[0] = xor_if(gamma0, xx0, (b0 ^ Y1));
    results[1] = xor_if(gamma1, xx1, (b1 ^ Y0)); 
  }






  template<typename T>
  void Simulator::multiply_mpc(const T& x0, bool b0, const T& x1, bool b1, T results[2], size_t mul, size_t cur_depth)
  {

    //std::cout << "x0 = " << x0 << std::endl;
    //std::cout << "b0 = " << b0 << std::endl;

    T D0; P0rand.next(D0);
    
    T D1; P1rand.next(D1);
    
    int dd0; 
    
    bool d0 = P0rand.next(dd0);

    int dd1; 
    
    bool d1 = P1rand.next(dd1); 

    T rand;
   
    P2rand.next(rand); 

    const T alpha = rand;

    T c0 = xor_if(alpha, D0, d1);
    T c1 = xor_if(alpha, D1, d0);

    T X0    = x0 ^ D0;
    bool Y0 = b0 ^ d0;

    T X1 = x1 ^ D1;
    bool Y1 = b1 ^ d1;



    const T gamma0 = xor_if(c0, X1, d0);

    T xx0 = x0;

    //std::cout << "gamma0 = " << gamma0 << std::endl;
  
    const T gamma1 = xor_if(c1, X0, d1);
    T xx1 = x1;

    results[0] = xor_if(gamma0, xx0, (b0 ^ Y1));
    results[1] = xor_if(gamma1, xx1, (b1 ^ Y0)); 

   if(cur_depth == 0)
    {
      P0_message.PB_root.blinds_sent[mul] = X0;
      P1_message.PB_root.blinds_recv[mul] = P0_message.PB_root.blinds_sent[mul];
      P0_message.PB_root.bit_blinds_sent[mul] = Y0;
      P1_message.PB_root.bit_blinds_recv[mul] = P0_message.PB_root.bit_blinds_sent[mul];
      


      P1_message.PB_root.blinds_sent[mul] = X1;
      P0_message.PB_root.blinds_recv[mul] = P1_message.PB_root.blinds_sent[mul];
      P1_message.PB_root.bit_blinds_sent[mul] = Y1;
      P0_message.PB_root.bit_blinds_recv[mul] = P1_message.PB_root.bit_blinds_sent[mul];

      P2_message[0].P2_root.blinds0[mul] = D0;
      P2_message[1].P2_root.blinds1[mul] = D1;
      P2_message[0].P2_root.bit_blinds0[mul] = d0;
      P2_message[1].P2_root.bit_blinds1[mul] = d1;

      P2_message[0].P2_root.gamma[mul] = gamma0;
      P0_message.PB_root.gamma[mul]    = gamma0; 
      P2_message[1].P2_root.gamma[mul] = gamma1;
      //P1_message.PB_root.gamma    = P2_message[1].P2_root.gamma[mul];
    }
    else
    {

      P0_message.PB_middle.blinds_sent[cur_depth][mul] = X0;
      P1_message.PB_middle.blinds_recv[cur_depth][mul] = P0_message.PB_middle.blinds_sent[cur_depth][mul];
      P0_message.PB_middle.bit_blinds_sent[cur_depth][mul] = Y0;
      P1_message.PB_middle.bit_blinds_recv[cur_depth][mul] = P0_message.PB_middle.bit_blinds_sent[cur_depth][mul];
      
      P1_message.PB_middle.blinds_sent[cur_depth][mul] = X1;
      P0_message.PB_middle.blinds_recv[cur_depth][mul] = P1_message.PB_middle.blinds_sent[cur_depth][mul];
      P1_message.PB_middle.bit_blinds_sent[cur_depth][mul] = Y1;
      P0_message.PB_middle.bit_blinds_recv[cur_depth][mul] = P1_message.PB_middle.bit_blinds_sent[cur_depth][mul];

      P2_message[0].P2_middle.blinds0[cur_depth][mul] = D0;
      P2_message[1].P2_middle.blinds1[cur_depth][mul] = D1;
      P2_message[0].P2_middle.bit_blinds0[cur_depth][mul] = d0;
      P2_message[1].P2_middle.bit_blinds1[cur_depth][mul] = d1;

      P2_message[0].P2_middle.gamma[cur_depth][mul] = gamma0;
      P2_message[1].P2_middle.gamma[cur_depth][mul] = gamma1;
    }
  }


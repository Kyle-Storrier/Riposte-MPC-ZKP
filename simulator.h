class  Simulator
{

  public:

  const size_t depth;

  std::vector<bool> P0direction;
  std::vector<bool> P1direction;
  std::vector<std::vector<block_t>> seed0;//[depth][2];
  std::vector<std::vector<block_t>> seed1;//[depth][2];
   
  std::vector<std::vector<bool>> bit0;//[depth][2];
  std::vector<std::vector<bool>> bit1;//[depth][2];
  
  Simulator(AES_KEY& aeskey, __m128i seed0_, __m128i seed1_, __m128i seed2_, size_t len, size_t depth_)
      : depth(depth_), P0rand(aeskey, seed0_, len), P1rand(aeskey, seed1_, len), P2rand(aeskey, seed2_, len)
  { 
    
    P0direction.resize(depth_+1);
    P1direction.resize(depth_+1);
    seed0.resize(depth_+1); for(size_t j = 0; j < depth_ + 1; ++j ) seed0[j].resize(2);
    seed1.resize(depth_+1); for(size_t j = 0; j < depth_ + 1; ++j ) seed1[j].resize(2);
    bit0.resize(depth_+1); for(size_t j = 0; j < depth_ + 1; ++j ) bit0[j].resize(2);
    bit1.resize(depth_+1); for(size_t j = 0; j < depth_ + 1; ++j ) bit1[j].resize(2);
  }


   void Gen_Blinds(block_t blinds0[], block_t blinds1[], block_t gamma[2][rounds]);

  
   void get_next_bits(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  bool L0_shares[2], bool R0_shares[2], bool L1_shares[2], bool R1_shares[2], size_t cur_depth);
  
  template<typename prgkey_t>
   auto prg_mpc(from_PB_ P0_view[], from_PB_ P1_view[], const prgkey_t & key,block_t seedL, block_t seedR,  const block_t blind0[rounds], const block_t blind1[rounds], block_t gamma[2][rounds], size_t cur_depth, size_t from, size_t len, bool LR);
   
   template<typename prgkey_t>
   void expand_mpc(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const prgkey_t & key,   block_t & seed0,  block_t & seed1,  block_t s0[2], uint8_t t0[2], 
                         block_t s1[2], uint8_t t1[2],   const block_t blind0[],  const block_t blind1[], block_t gamma[2][rounds], size_t cur_depth, bool LR);
  
   void conditional_swap_and_next_seed(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const prgkey_t & key, block_t L0_shares[2], block_t R0_shares[2], block_t L1_shares[2], block_t R1_shares[2], block_t left_out[2], block_t right_out[2], size_t cur_depth);
   
   template<typename leaf_t, typename node_t, typename prgkey_t>
   void middle_layers(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const prgkey_t& key, dpf_key<leaf_t, node_t, prgkey_t> dpfkey0, dpf_key<leaf_t, node_t, prgkey_t> dpfkey1, size_t cur_depth);  

   template<typename leaf_t, typename node_t, typename prgkey_t>
   void leaf_layer(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const prgkey_t& key, dpf_key<leaf_t, node_t, prgkey_t> dpfkey0, dpf_key<leaf_t, node_t, prgkey_t> dpfkey1, size_t cur_depth);  


   template<typename T>
   void multiply_mpc(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const T&, bool, const T&, bool, T results[2], size_t mul, size_t cur_depth);
   void multiply_mpc_b(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth);
   void multiply_mpc(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth);
   // void dpf_mpc(LowMC&, AES_KEY&, dpf_key<__mX, nitems>& k0, dpf_key<__mX, nitems> k1); 
   template<typename leaf_t, typename node_t, typename prgkey_t>
   void root_layer(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const prgkey_t& key, dpf_key<leaf_t, node_t, prgkey_t> dpfkey0, dpf_key<leaf_t, node_t, prgkey_t> dpfkey1);
 

  ~ Simulator()
  {

  }

//private:
  MPCrandomness P0rand, P1rand, P2rand;
  
};



  void Simulator::Gen_Blinds(block_t blinds0[rounds], block_t blinds1[rounds], block_t gamma[2][rounds])
  {
    block_t rand[rounds];
    for (unsigned r = 0; r < rounds; ++r)
     {        
        blinds0[r] = P0rand.next_block();
        blinds1[r] = P1rand.next_block();
        rand[r]    = P2rand.next_block();

        const block_t tmp1 = ((blinds0[r] >> 1) & blinds1[r]) ^ ((blinds1[r] >> 1) & blinds0[r]);
        const block_t tmp2 = ((blinds0[r] >> 2) & blinds1[r]) ^ ((blinds1[r] >> 2) & blinds0[r]);
    
        const block_t bc = (tmp1 << 2) & maska;
        const block_t ac = (tmp2 << 1) & maskb;
        const block_t ab = (tmp1 >> 1) & maskc;
    
        gamma[0][r] = (bc | ac | ab) ^ rand[r];
        gamma[1][r] = rand[r];
     }
  }

    void Simulator::multiply_mpc_b(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth)
  {
    bool D0 = P0rand.next_bool(); 
    bool D1 = P1rand.next_bool(); 
     
    bool d0 = P0rand.next_bool();     
    bool d1 = P1rand.next_bool();
      
    const bool alpha = P2rand.next_bool(); 

     bool c0 = xor_if(alpha, D0, d1);
    bool c1 = xor_if(alpha, D1, d0);
 
    from_P2_to_P0.c_bit[cur_depth][mul]       = c0;
    from_P2_to_P1.c_bit[cur_depth][mul]       = c1;
 
    P0_view[cur_depth].next_bit_L_recv[mul] = x0 ^ D0;
    P0_view[cur_depth].next_bit_R_recv[mul] = b0 ^ d0;

    P1_view[cur_depth].next_bit_L_recv[mul] = x1 ^ D1;
    P1_view[cur_depth].next_bit_R_recv[mul] = b1 ^ d1;

    const bool gamma0 = xor_if(c0, P1_view[cur_depth].next_bit_L_recv[mul], d0);
    bool xx0 = x0;   
 

    const bool gamma1 = xor_if(c1, P0_view[cur_depth].next_bit_L_recv[mul], d1);
    bool xx1 = x1;  
  
    results[0] = xor_if(gamma0, xx0, (b0 ^ P1_view[cur_depth].next_bit_R_recv[mul]));
  
    results[1] = xor_if(gamma1, xx1, (b1 ^ P0_view[cur_depth].next_bit_R_recv[mul])); 
  }


  void Simulator::get_next_bits(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  bool bit0_L[2], bool bit0_R[2], bool bit1_L[2], bool bit1_R[2], size_t cur_depth)
  {
    
    bool notb_bit0_L[2];    
    multiply_mpc_b(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, bit0_L[0], P0direction[cur_depth], bit0_L[1], !P1direction[cur_depth], notb_bit0_L, 0, cur_depth); // L_0 * (direction - 1)


    bool notb_bit1_L[2]; 
    multiply_mpc_b(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, bit1_L[0], P0direction[cur_depth], bit1_L[1], !P1direction[cur_depth], notb_bit1_L, 1, cur_depth); // L_1 * (direction - 1)

    bool b_bit0_R[2];    
    multiply_mpc_b(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, bit0_R[0], P0direction[cur_depth], bit0_R[1],  P1direction[cur_depth], b_bit0_R, 2, cur_depth);

    bool b_bit1_R[2];
    multiply_mpc_b(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, bit1_R[0], P0direction[cur_depth], bit1_R[1],  P1direction[cur_depth], b_bit1_R, 3, cur_depth); // R_1 * direction    
    
    
    bit0[cur_depth+1][0] = notb_bit0_L[0] ^ b_bit0_R[0];
    
    bit0[cur_depth+1][1] = notb_bit0_L[1] ^ b_bit0_R[1];    
    
    bit1[cur_depth+1][0] = notb_bit1_L[0] ^ b_bit1_R[0];    
    
    bit1[cur_depth+1][1] = notb_bit1_L[1] ^ b_bit1_R[1];

  }


 void Simulator::multiply_mpc(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth)
  {

    //std::cout << "Simulator:::multiply_mpc, cur_depth = " << cur_depth << std::endl; 

    bool D0 = P0rand.next_bool(); 
    bool D1 = P1rand.next_bool(); 
     
    bool d0 = P0rand.next_bool();     
    bool d1 = P1rand.next_bool();
    


    const bool alpha = P2rand.next_bool(); 

    bool c0 = xor_if(alpha, D0, d1);
    bool c1 = xor_if(alpha, D1, d0);

    from_P2_to_P0.c_bit2[cur_depth][mul] = c0;
    from_P2_to_P1.c_bit2[cur_depth][mul] = c1;


 

    P0_view[cur_depth].next_bit_L_recv2[mul] = x0 ^ D0;
   // std::cout << "sim: x0 ^ D0 = " << x0 << " ^ " << D0 << std::endl;
    P0_view[cur_depth].next_bit_R_recv2[mul] = b0 ^ d0;

 

    P1_view[cur_depth].next_bit_L_recv2[mul] = x1 ^ D1;
   //  std::cout << "sim: x1 ^ D1 = " << x1 << " ^ " << D1 << std::endl;
    P1_view[cur_depth].next_bit_R_recv2[mul] = b1 ^ d1;


      const bool gamma0 = xor_if(c0, P1_view[cur_depth].next_bit_L_recv2[mul], d0);
  
    bool xx0 = x0;   
    
    
    const bool gamma1 = xor_if(c1, P0_view[cur_depth].next_bit_L_recv2[mul], d1);
   
    bool xx1 = x1;  
   
    results[0] = xor_if(gamma0, xx0, (b0 ^ P1_view[cur_depth].next_bit_R_recv2[mul]));
  
    results[1] = xor_if(gamma1, xx1, (b1 ^ P0_view[cur_depth].next_bit_R_recv2[mul])); 

 
  }



  template<typename T>
  void Simulator::multiply_mpc(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1,  const T& x0, bool b0, const T& x1, bool b1, T results[2], size_t mul, size_t cur_depth)
  {   

    T D0 = P0rand.next_block();   
    T D1 = P1rand.next_block();
  
    bool d0 = P0rand.next_bool();   
    bool d1 = P1rand.next_bool(); 

    const T alpha = P2rand.next_block();

    T c0 = xor_if(alpha, D0, d1);
    T c1 = xor_if(alpha, D1, d0);
 
    from_P2_to_P0.c[cur_depth][mul] = c0;
    from_P2_to_P1.c[cur_depth][mul] = c1;

   
 
    P0_view[cur_depth].blinds_recv[mul]     = x0 ^ D0; 
    P0_view[cur_depth].bit_blinds_recv[mul] = b0 ^ d0;
 
    P1_view[cur_depth].blinds_recv[mul]     = x1 ^ D1;    
    P1_view[cur_depth].bit_blinds_recv[mul] = b1 ^ d1;

    const T gamma0 = xor_if(c0, P1_view[cur_depth].blinds_recv[mul], d0);

    T xx0 = x0;

    const T gamma1 = xor_if(c1, P0_view[cur_depth].blinds_recv[mul], d1);
    
    T xx1 = x1;

    results[0] = xor_if(gamma0, xx0, (b0 ^ P1_view[cur_depth].bit_blinds_recv[mul])); 
    results[1] = xor_if(gamma1, xx1, (b1 ^ P0_view[cur_depth].bit_blinds_recv[mul])); 
  }

template<typename prgkey_t>
inline auto Simulator::prg_mpc(from_PB_ P0_view[], from_PB_ P1_view[], const prgkey_t& key, block_t seedL, block_t seedR,  const block_t blind0[rounds], const block_t blind1[rounds], block_t gamma[2][rounds], size_t cur_depth, size_t from, size_t len, bool LR)
{    

  __m256i seed0_ = clear_lsb(seedL, 0b11);
  __m256i seed1_ = clear_lsb(seedR, 0b11);
 
  //for(size_t i = 0; i < len; ++i)
  {

    auto tmp0 = _mm256_xor_si256(seed0_, _mm256_set_epi64x(0, 0, 0, 0));    
    auto tmp1 = seed1_;
        
    block_t encrypt_outL_0[rounds+1];
    block_t encrypt_outL_1[rounds+1];
    
    auto encrypt_outL =  key.encrypt_MPC_proof(tmp0, tmp1, blind0, blind1, gamma, encrypt_outL_0, encrypt_outL_1 );

 
    
    if(!LR)
    {
      for(size_t r = 0; r <= rounds; ++r) P1_view[cur_depth].seed0L_encrypt[r] = encrypt_outL_1[r];
  
      for(size_t r = 0; r <= rounds; ++r)  P0_view[cur_depth].seed0L_encrypt[r] = encrypt_outL_0[r];
     
    }

    if(LR)
    {

     for(size_t r = 0; r <= rounds; ++r) P1_view[cur_depth].seed1L_encrypt[r] = encrypt_outL_1[r]; 
     for(size_t r = 0; r <= rounds; ++r) P0_view[cur_depth].seed1L_encrypt[r] = encrypt_outL_0[r]; 
    }

    encrypt_outL.first[rounds] = _mm256_xor_si256(encrypt_outL.first[rounds] , tmp0);
    encrypt_outL.second[rounds] = _mm256_xor_si256(encrypt_outL.second[rounds] , tmp1);

    auto tmp = _mm256_xor_si256(tmp0, tmp1);
 
    tmp0 = _mm256_xor_si256(seed0_, _mm256_set_epi64x( 0,  0,  0, 1 ));    
    tmp1 = seed1_;


    block_t encrypt_outR_0[rounds+1];
    block_t encrypt_outR_1[rounds+1];
    auto encrypt_outR =  key.encrypt_MPC_proof(tmp0, tmp1, blind0, blind1, gamma, encrypt_outR_0, encrypt_outR_1);

    if(!LR)
    {

    for(size_t r = 0; r <= rounds; ++r)P1_view[cur_depth].seed0R_encrypt[r] = encrypt_outR_1[r];

    for(size_t r = 0; r <= rounds; ++r) P0_view[cur_depth].seed0R_encrypt[r] = encrypt_outR_0[r]; 
    }

    if(LR)
    {
    
     for(size_t r = 0; r <= rounds; ++r) P1_view[cur_depth].seed1R_encrypt[r] = encrypt_outR_1[r];     
     for(size_t r = 0; r <= rounds; ++r) P0_view[cur_depth].seed1R_encrypt[r] = encrypt_outR_0[r]; 
    }

    encrypt_outR.first[rounds] = _mm256_xor_si256(encrypt_outR.first[rounds] , tmp0);
    encrypt_outR.second[rounds] = _mm256_xor_si256(encrypt_outR.second[rounds] , tmp1);

    return std::make_pair(std::move(encrypt_outL), std::move(encrypt_outR));
  }
}

  template<typename prgkey_t>
 inline void Simulator::expand_mpc(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1, const prgkey_t & key,  block_t & seed0,  block_t & seed1,  block_t s0[2], uint8_t t0[2], 
                        block_t s1[2], uint8_t t1[2],   const block_t blind0[rounds],   const block_t blind1[rounds], block_t gamma[2][rounds], size_t cur_depth, bool LR)
 {
 
    auto prg_out = prg_mpc(P0_view, P1_view, key, seed0, seed1, blind0, blind1, gamma,  cur_depth, 0, 2, LR);
    
    auto encrypt_outL = prg_out.first;
    auto encrypt_outR = prg_out.second;
    
    s0[L] = encrypt_outL.first[rounds];
    s0[R] = encrypt_outR.first[rounds];

    s1[L] = encrypt_outL.second[rounds];
    s1[R] = encrypt_outR.second[rounds];


    t0[L] = get_lsb(s0[L]);
    t1[L] = get_lsb(s1[L]);
    s0[L] = clear_lsb(s0[L], 0b11);
    s1[L] = clear_lsb(s1[L], 0b11);

    t0[R] = get_lsb(s0[R]);
    t1[R] = get_lsb(s1[R]);
    s0[R] = clear_lsb(s0[R], 0b11);
    s1[R] = clear_lsb(s1[R], 0b11); 
 }



   void Simulator::conditional_swap_and_next_seed(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1, const prgkey_t& key, block_t s0_L[2], block_t s0_R[2], block_t s1_L[2], block_t s1_R[2], block_t left_out[2], block_t right_out[2], size_t cur_depth)
  {
     
     block_t bs0_L[2]; // shares of b \times s0_L
     block_t notbs0_L[2]; // shares of (1 - b)  \times s0_L
     block_t bs1_L[2]; // shares of b \times s1_L
     block_t notbs1_L[2]; // shares of (1 - b) \times s1_L
     block_t bs0_R[2]; // shares of b \times s0_R
     block_t bs1_R[2]; // shares of b \times s1_R
     block_t notbs1_R[2];  // shares of (1-b) \times s1_R
     block_t notbs0_R[2]; // shares of (1-b) \times s0_R

    std::cout << "current direction: " << P0direction[cur_depth] << std::endl;
    std::cout << "current direction: " << P1direction[cur_depth] << std::endl;
    std::cout << "cur_depth: " << cur_depth << std::endl;
    std::cout << "depth: " << depth << std::endl;
    
     multiply_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, s0_L[0], P0direction[cur_depth], s0_L[1], P1direction[cur_depth], bs0_L, 0, cur_depth); // L_0 * direction
     notbs0_L[0] = s0_L[0] ^ bs0_L[0];
     notbs0_L[1] = s0_L[1] ^ bs0_L[1];
 
     
     multiply_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, s1_L[0], P0direction[cur_depth], s1_L[1], P1direction[cur_depth], bs1_L, 1, cur_depth); // L_1 * direction     
     notbs1_L[0] = s1_L[0] ^ bs1_L[0];
     notbs1_L[1] = s1_L[1] ^ bs1_L[1];
     
     multiply_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, s0_R[0], P0direction[cur_depth], s0_R[1], P1direction[cur_depth], bs0_R, 2, cur_depth); // R_0 * direction      
     notbs0_R[0] = s0_R[0] ^ bs0_R[0];
     notbs0_R[1] = s0_R[1] ^ bs0_R[1];

     multiply_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, s1_R[0], P0direction[cur_depth], s1_R[1], P1direction[cur_depth], bs1_R, 3, cur_depth); // R_1 * direction
     notbs1_R[0] = s1_R[0] ^ bs1_R[0];
     notbs1_R[1] = s1_R[1] ^ bs1_R[1];

     right_out[0] = notbs0_L[0] ^ bs0_R[0] ^ notbs0_L[1] ^ bs0_R[1];
     right_out[1] = notbs1_L[0] ^ bs1_R[0] ^ notbs1_L[1] ^ bs1_R[1];

     left_out[0] = bs0_L[0] ^ notbs0_R[0] ^ bs1_L[0] ^ notbs1_R[0];
     left_out[1] = bs0_L[1] ^ notbs0_R[1] ^ bs1_L[1] ^ notbs1_R[1];
    
     seed0[cur_depth+1][0] = notbs0_L[0] ^ bs0_R[0];
      
     seed0[cur_depth+1][1] = notbs0_L[1] ^ bs0_R[1];
     seed1[cur_depth+1][0] = notbs1_L[0] ^ bs1_R[0];

     seed1[cur_depth+1][1] = notbs1_L[1] ^ bs1_R[1];
  }

  template<typename leaf_t, typename node_t, typename prgkey_t>
  void Simulator::root_layer(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1, const prgkey_t & prgkey, dpf_key<leaf_t, node_t, prgkey_t> dpfkey0, dpf_key<leaf_t, node_t, prgkey_t> dpfkey1)
  {
    std::cout << "root_layer: " << std::endl;
    size_t layer = 0;
    
    const block_t seed0 = dpfkey0.root; 
    const block_t seed1 = dpfkey1.root;



    auto & cw0 = dpfkey0.cw[layer];
    uint8_t cw_t0[2] = { get_lsb(cw0, 0b01), get_lsb(cw0, 0b10) };

    
    auto & cw1 = dpfkey1.cw[layer];
    uint8_t cw_t1[2] = { get_lsb(cw1, 0b01), get_lsb(cw1, 0b10) };
    

    size_t cur_depth = 0;
   
 
   
    bool b0 = get_lsb(dpfkey0.root);
    
    bool b1 = get_lsb(dpfkey1.root);

    block_t s0[2], s1[2];    uint8_t t0[2], t1[2];

    block_t child0[2];

    
    expand(prgkey, seed0, child0, t0);
  
    // s0[L] and s0[R] are P0's left and right children after the correction word is applied
    s0[L] = clear_lsb(xor_if(child0[L], cw0, !b0), 0b11); 
    s0[R] = clear_lsb(xor_if(child0[R], cw0, !b0), 0b11); 
   
    block_t child1[2];
    
    
    expand(prgkey, seed1,  child1 , t1); 

    s1[L] = clear_lsb(xor_if(child1[L], cw0, !b1), 0b11);
    s1[R] = clear_lsb(xor_if(child1[R], cw0, !b1), 0b11);


    block_t s0_L[2]; 
    block_t s0_R[2];  

    s0_L[0] = P0rand.next_block();
    s0_L[1] = s0_L[0] ^ s0[L];

    s0_R[0] = P0rand.next_block(); 
    s0_R[1] = s0_R[0] ^ s0[R];
  

    block_t s1_L[2]; 
    block_t s1_R[2];
   
    
    s1_L[0] = P1rand.next_block();    
    s1_L[1] = s1_L[0] ^ s1[L];
  
    s1_R[0] = P1rand.next_block();
    s1_R[1] = s1_R[0] ^ s1[R];

    
    P1_view[0].L_shares_recv = s1_L[0];
    P1_view[0].R_shares_recv = s1_R[0];

 
    P0_view[0].L_shares_recv = s0_L[1];
    P0_view[0].R_shares_recv = s0_R[1];


    block_t left_out[2], right_out[2];
    
    conditional_swap_and_next_seed(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, prgkey, s0_L, s0_R, s1_L, s1_R, left_out, right_out, cur_depth);


    std::cout << "reconstructed L = " << (left_out[0].bits  ^ left_out[1].bits) << std::endl;
    std::cout << "reconstructed R = " << (right_out[0].bits ^ right_out[1].bits) << std::endl << std::endl;
    

    // b0_L[2] are the shares of the cwbit for L-child of P0's root
    bool bit0_L[2]; 

    // b0_R[2] are the shares of the cwbit for R-child of P0's root
    bool bit0_R[2]; 

    bit0_L[0] = P0rand.next_bool();
    bit0_L[1] = bit0_L[0] ^ (t0[L] ^ cw_t0[L] & (b0));
    
    bit0_R[0] = P0rand.next_bool();
    bit0_R[1] = bit0_R[0] ^ (t0[R] ^ cw_t0[R] & (b0));
    
    // // b1_L[2] are the shares of the cwbit for L-child of P1's root
    bool bit1_L[2];
   
    // // b1_R[2] are the shares of the cwbit for R-child of P1's root
    bool bit1_R[2];

    bit1_L[0] = P1rand.next_bool();
    bit1_L[1] = bit1_L[0] ^ (t1[L] ^ cw_t1[L] & (b1));
    
    bit1_R[0] = P1rand.next_bool();
    bit1_R[1] = bit1_R[0] ^ (t1[R] ^ cw_t1[R] & (b1));

    P1_view[0].bit_L_shares_recv = bit1_L[0];
    P1_view[0].bit_R_shares_recv = bit1_R[0];
 
    P0_view[0].bit_L_shares_recv = bit0_L[1];
    P0_view[0].bit_R_shares_recv = bit0_R[1];

    get_next_bits(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, bit0_L, bit0_R, bit1_L, bit1_R, cur_depth);

  }

  template<typename leaf_t, typename node_t, typename prgkey_t>
  void Simulator::middle_layers(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1, const prgkey_t& key, dpf_key<leaf_t, node_t, prgkey_t> dpfkey0, dpf_key<leaf_t, node_t, prgkey_t> dpfkey1, size_t cur_depth) 
  {
    
    auto & cw0 = dpfkey0.cw[cur_depth];
    uint8_t cw_t0[2] = { get_lsb(cw0, 0b01), get_lsb(cw0, 0b10) };

    auto & cw1 = dpfkey1.cw[cur_depth];
    uint8_t cw_t1[2] = { get_lsb(cw1, 0b01), get_lsb(cw1, 0b10) };

     auto & final0 = dpfkey0.finalizer;
     auto & final1 = dpfkey1.finalizer;

    block_t P0gamma[2][rounds];
    block_t P0blinds0[rounds];
    block_t P0blinds1[rounds];

    Gen_Blinds(P0blinds0, P0blinds1, P0gamma);
    

    for(size_t j = 0; j < rounds; ++j)
    {
      from_P2_to_P0.gamma0[cur_depth][j] = P0gamma[0][j];
      from_P2_to_P1.gamma0[cur_depth][j] = P0gamma[1][j];
    } 
    
    block_t s0[2], s1[2], s_[2];
    uint8_t t0[2], t1[2], t_[2];    

    expand_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, key, seed0[cur_depth][0], seed0[cur_depth][1],  s0, t0, s1, t1, P0blinds0, P0blinds1, P0gamma, cur_depth, false);
        
    block_t s0_L[2];
    block_t s0_R[2];
  
    s0_L[0] = clear_lsb(xor_if(s0[L], cw0, bit0[cur_depth][0]), 0b11);
    s0_L[1] = clear_lsb(xor_if(s1[L], cw1, !bit0[cur_depth][1]), 0b11);
    s0_R[0] = clear_lsb(xor_if(s0[R], cw0, bit0[cur_depth][0]), 0b11);
    s0_R[1] = clear_lsb(xor_if(s1[R], cw1, !bit0[cur_depth][1]), 0b11);


    block_t ss0[2], ss1[2];

    uint8_t tt0[2], tt1[2]; 

    block_t P1gamma[2][rounds];
    block_t P1blinds0[rounds];
    block_t P1blinds1[rounds];

    Gen_Blinds(P1blinds0, P1blinds1, P1gamma);
    
    for(size_t j = 0; j < rounds; ++j)
    {
      from_P2_to_P0.gamma1[cur_depth][j] = P1gamma[0][j];
      from_P2_to_P1.gamma1[cur_depth][j] = P1gamma[1][j]; 
    } 
 
    expand_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, key, seed1[cur_depth][0], seed1[cur_depth][1],  ss0, tt0, ss1, tt1, P1blinds0, P1blinds1, P1gamma, cur_depth, true);


    block_t s1_L[2];
    block_t s1_R[2];

    s1_L[0] = clear_lsb(xor_if(ss0[L], cw0, bit1[cur_depth][0]), 0b11);
    s1_L[1] = clear_lsb(xor_if(ss1[L], cw1, !bit1[cur_depth][1]), 0b11);
    s1_R[0] = clear_lsb(xor_if(ss0[R], cw0, bit1[cur_depth][0]), 0b11);
    s1_R[1] = clear_lsb(xor_if(ss1[R], cw1, !bit1[cur_depth][1]), 0b11);

 
    block_t left_out[2]; 
 
    block_t right_out[2];
 
    conditional_swap_and_next_seed(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, key, s0_L, s0_R, s1_L, s1_R, left_out, right_out, cur_depth);
    
    std::cout << "reconstructed L = " << (left_out[0].bits  ^ left_out[1].bits) << std::endl;
    std::cout << "reconstructed R = " << (right_out[0].bits  ^ right_out[1].bits) << std::endl << std::endl;
    
    bool bit0_L[2];     
    multiply_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, cw_t0[L],  false, false, bit0[cur_depth][1], bit0_L, 0, cur_depth);
    bit0_L[0] = bit0_L[0] ^ (t0[0]) ^ cw_t0[L] & (bit0[cur_depth][0]);
    bit0_L[1] = bit0_L[1] ^ t1[0];

 
    bool bit0_R[2];         
    multiply_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, cw_t0[R],  false, false, bit0[cur_depth][1], bit0_R, 1, cur_depth);
    bit0_R[0] = bit0_R[0] ^ t0[1] ^ (cw_t0[R] & bit0[cur_depth][0]);
    bit0_R[1] = bit0_R[1] ^ t1[1];

 
    bool bit1_L[2]; 
    multiply_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, bit1[cur_depth][0], false, false, cw_t1[L],  bit1_L, 2, cur_depth);
    bit1_L[0] = bit1_L[0] ^ tt0[0]; 
    bit1_L[1] = bit1_L[1] ^  tt1[0] ^ (cw_t1[L] & bit1[cur_depth][1]);


      
    bool bit1_R[2];     
    multiply_mpc(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, bit1[cur_depth][0],  false, false, cw_t1[R],  bit1_R, 3, cur_depth);
    bit1_R[0] = bit1_R[0] ^  tt0[1];
    bit1_R[1] = bit1_R[1] ^ (tt1[1]) ^ (cw_t1[R] & bit1[cur_depth][1]);

    get_next_bits(P0_view, P1_view, from_P2_to_P0, from_P2_to_P1, bit0_L, bit0_R, bit1_L, bit1_R, cur_depth);
  


    if(cur_depth == depth - 1)
    {
 
      for(size_t j = 0; j < dpfkey0.nodes_per_leaf; ++j)
      {
 
         seed0[cur_depth+1][0] = xor_if(seed0[cur_depth+1][0], final0[j], bit0[cur_depth+1][0]);
         seed0[cur_depth+1][1] = xor_if(seed0[cur_depth+1][1], final0[j], bit0[cur_depth+1][1]);

         seed1[cur_depth+1][0] = xor_if(seed1[cur_depth+1][0], final1[j], bit1[cur_depth+1][0]);
         seed1[cur_depth+1][1] = xor_if(seed1[cur_depth+1][1], final1[j], bit1[cur_depth+1][1]);

         P0_view[cur_depth].final = seed0[cur_depth+1][0] ^ seed1[cur_depth+1][0];
         P1_view[cur_depth].final = seed0[cur_depth+1][1] ^ seed1[cur_depth+1][1];
      }

      block_t final_out = (seed0[cur_depth+1][0] ^ seed0[cur_depth+1][1] ^ seed1[cur_depth+1][0] ^ seed1[cur_depth+1][1]);
      std::cout << " -> [Simulator]: " << final_out.bits << std::endl;
    }    
 }



 template<typename leaf_t, typename node_t, typename prgkey_t>
  void Simulator::leaf_layer(from_PB_ P0_view[], from_PB_ P1_view[], from_P2 &from_P2_to_P0, from_P2 &from_P2_to_P1, const prgkey_t& key, dpf_key<leaf_t, node_t, prgkey_t> dpfkey0, dpf_key<leaf_t, node_t, prgkey_t> dpfkey1, size_t cur_depth) 
  {
    
    auto & cw0 = dpfkey0.cw[cur_depth];
    uint8_t cw_t0[2] = { get_lsb(cw0, 0b01), get_lsb(cw0, 0b10) };

    auto & cw1 = dpfkey1.cw[cur_depth];
    uint8_t cw_t1[2] = { get_lsb(cw1, 0b01), get_lsb(cw1, 0b10) };


     
 }
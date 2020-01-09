class Verifier
{

  public:

  std::vector<bool> Pdirection;
  std::vector<bool> Pdirection_other;


  // block_t seed0_[depth][2];
  // block_t seed1_[depth][2];

  std::vector<std::vector<bool>> bit0_next;//[depth][2];
  std::vector<std::vector<bool>> bit1_next;//[depth][2];



  std::vector<block_t> seed0;
  std::vector<block_t> seed1;

  std::vector<bool> bit0;
  std::vector<bool> bit1;
  //size_t rounds;
  Verifier(AES_KEY& aeskey, __m128i seed, size_t len, size_t depth)
      : PBrand(aeskey, seed, len),
        PBrand_other(aeskey, seed, len),
        P2rand(aeskey, seed, len)
  {

     
    Pdirection.resize(depth+1);
    Pdirection_other.resize(depth+1);

    bit0_next.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) bit0_next[j].resize(2);
    bit1_next.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) bit1_next[j].resize(2);
    seed0.resize(depth+1);
    seed1.resize(depth+1);
    bit0.resize(depth+1);
    bit1.resize(depth+1);
  }

  // Verifier(AES_KEY& aeskey, __m128i seed0, __m128i seed1, __m128i seed2 , size_t len, size_t depth)
  //     : PBrand(aeskey, seed0, len),
  //       PBrand_other(aeskey, seed1, len),
  //       P2rand(aeskey, seed2, len)
  // {
 
  //   Pdirection.resize(depth+1);
  //   Pdirection_other.resize(depth+1);
  // }  



 
  /// Dice role i \neq n
  // template<typename KEY_TYPE, typename __mX>  
  // void root_layer(KEY_TYPE& key,  dpf_key<__mX, nitems> dpfkey, dpf_key<__mX, nitems> dpfkey_other, bool party);
  // void conditional_swap_and_next_seed(const from_P2 from_P2_to_PB, const from_PB from_PB_, LowMC& key, block_t L0_shares[2], block_t R0_shares[2], block_t L1_shares[2], block_t R1_shares[2], block_t left_out[2], block_t right_out[2], size_t cur_depth);
  // template<typename T>
  // void multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_,const T&, bool, const T&, bool, T results[2], size_t mul, size_t cur_depth);

  // void get_next_bits(const from_P2 from_P2_to_PB, const from_PB from_PB_, bool L0_shares[2], bool R0_shares[2], bool L1_shares[2], bool R1_shares[2], size_t cur_depth);
  // void Gen_Blinds(block_t blinds0[rounds], block_t blinds1[rounds], block_t gamma[2][rounds]);
  // void multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_,const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth);
  // void multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_,const bool& x0, bool b0, const bool& x1, bool b1, bool results[2]);
  // //void middle_layers(const from_P2 from_P2_to_PB, const from_PB from_PB_, LowMC& key, dpf_key<__mX, nitems> dpfkey, dpf_key<__mX, nitems> dpfkey_other, size_t cur_depth);
  
  // template<typename KEY_TYPE, typename __mX>
  // void expand_mpc(KEY_TYPE & key, const block_t & seed0, const block_t & seed1,  block_t s0[2], bool t0[2], 
  //                       block_t s1[2], bool t1[2], const block_t blind0[rounds], const block_t blind1[rounds], block_t gamma[2][rounds], size_t cur_depth, bool LR);

  // ///// Functions for dice role i = n
   void get_next_bits(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, bool L0_shares, bool R0_shares, bool L1_shares, bool R1_shares,  bool party ,size_t cur_depth);
  
   void Gen_Blinds(block_t blinds0[]);
  
   //template<typename leaf_t, typename node_t, typename prgkey_t>
   void expand_mpc_verify(const LowMC<__m256i>& key, const block_t seed, const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, block_t &s0_L, block_t & s0_R, 
                                         uint8_t & t0_L, uint8_t & t0_R, block_t blind[], block_t gamma[], bool party, size_t cur_depth, bool LR);
  
   template<typename leaf_t, typename node_t, typename prgkey_t>
   void root_layer(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other,  const LowMC<__m256i>& key,  dpf_key<leaf_t, node_t, prgkey_t> dpfkey,  bool party);
  

   void multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, block_t val_mX, bool val_bool, size_t cur_depth, size_t mul, bool party, block_t& out);
  
   void multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, bool val, bool val_bool, bool blind_mX, bool blind_bool,  bool party, size_t mul , bool& result, size_t cur_depth);
  
   void multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, bool val, bool party, size_t mul , bool& result, size_t cur_depth);
  
  template<typename prgkey_t>
  void conditional_swap_and_next_seed(const from_P2 from_P2_to_PB,  const from_PB from_PB_, from_PB & from_PB_other, block_t L_share, block_t R_share, block_t L1_share, 
                                      block_t R1_share, prgkey_t & key,  bool party, size_t cur_depth);
  
  template<typename leaf_t, typename node_t, typename prgkey_t>
  void middle_layers(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, const LowMC<__m256i> & key, const block_t seed0, const  block_t seed1, dpf_key<leaf_t, node_t, prgkey_t> dpfkey, size_t cur_depth, bool party);

  // void verify_P2();

  ~ Verifier()
  {

  }

  //private:
  MPCrandomness PBrand;
  MPCrandomness PBrand_other;
  MPCrandomness P2rand;
};


  void Verifier::Gen_Blinds(block_t blinds0[])
  {
    for (unsigned r = 0; r < rounds; ++r)
     {
        PBrand.next(blinds0[r]); 
     }
  }



  

  void Verifier::get_next_bits(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, bool L0_shares, bool R0_shares, bool L1_shares, bool R1_shares, bool party, size_t cur_depth)
  {    
    bool blind_L,  blind_R;

    blind_L = from_PB_.next_bit_L_recv[cur_depth][0];
    blind_R = from_PB_.next_bit_R_recv[cur_depth][0];
 
    bool bNotL0_share;
   
    if(!party) multiply_mpc(from_P2_to_PB, from_PB_,from_PB_other , L0_shares, Pdirection[cur_depth], blind_L, blind_R,   party, 0, bNotL0_share, cur_depth);
    if(party)  multiply_mpc(from_P2_to_PB, from_PB_,from_PB_other , L0_shares, !Pdirection[cur_depth], blind_L, blind_R,   party, 0, bNotL0_share, cur_depth);
 
    blind_L = from_PB_.next_bit_L_recv[cur_depth][1];
    blind_R = from_PB_.next_bit_R_recv[cur_depth][1];

    bool bNotL1_share;
   
    if(!party) multiply_mpc(from_P2_to_PB, from_PB_, from_PB_other , L1_shares, Pdirection[cur_depth] , blind_L, blind_R,   party, 1, bNotL1_share, cur_depth);
    if(party ) multiply_mpc(from_P2_to_PB, from_PB_, from_PB_other , L1_shares, !Pdirection[cur_depth] , blind_L, blind_R,  party, 1, bNotL1_share, cur_depth);
    
    blind_L = from_PB_.next_bit_L_recv[cur_depth][2];
    blind_R = from_PB_.next_bit_R_recv[cur_depth][2];
 
    bool bR0_share;
    multiply_mpc(from_P2_to_PB, from_PB_, from_PB_other , R0_shares, Pdirection[cur_depth], blind_L, blind_R,  party, 2, bR0_share, cur_depth);

    blind_L = from_PB_.next_bit_L_recv[cur_depth][3];
    blind_R = from_PB_.next_bit_R_recv[cur_depth][3];
 
    bool bR1_share;
    multiply_mpc(from_P2_to_PB, from_PB_, from_PB_other , R1_shares, Pdirection[cur_depth], blind_L, blind_R,   party, 3, bR1_share, cur_depth);

    bit0[cur_depth+1] = bNotL0_share  ^ bR0_share;    
    bit1[cur_depth+1] = bNotL1_share  ^ bR1_share; 
    
  }

//   //multiply_mpc
  void Verifier::multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, bool val_mX,  bool party, size_t mul , bool& result, size_t cur_depth)
  {

    bool val_bool = false;
     
    bool D0 = PBrand.next_bool();
    bool X0 = val_mX ^ D0;  ;
         
    bool d0 = PBrand.next_bool();
    bool Y0 = val_bool ^ d0;   
    
    bool X1 = from_PB_.next_bit_L_recv2[cur_depth][mul];
    bool Y1 = from_PB_.next_bit_R_recv2[cur_depth][mul];
    bool c  = from_P2_to_PB.c_bit2[cur_depth][mul];
 
    from_PB_other.next_bit_L_recv2[cur_depth][mul] = X0;
    from_PB_other.next_bit_R_recv2[cur_depth][mul] = Y0;
    
    if(party)
    {
     from_PB_other.next_bit_L_recv2[cur_depth][mul] = (val_bool ^ D0);
     from_PB_other.next_bit_R_recv2[cur_depth][mul] = val_mX ^ d0;
    }

    bool other_blinded = from_PB_.next_bit_L_recv2[cur_depth][mul];
     
    bool gamma = xor_if(c, other_blinded, d0);
      
    if(!party) result = xor_if(gamma, val_mX, (val_bool ^ Y1));
    if( party) result = xor_if(gamma, val_bool, (val_mX ^ X1));
  }
  

  
  void Verifier::multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, bool val_mX, bool val_bool, bool blind_mX, bool blind_bool,  bool party, size_t mul, bool& result, size_t cur_depth)
  {
    // int dd;
    bool D0 = PBrand.next_bool();
    bool X0 = val_mX ^ D0; 


    // int dd0;     
  
    bool d0 = PBrand.next_bool();
    bool Y0 = val_bool ^ d0;  
    
    
    bool X1 = from_PB_.next_bit_L_recv[cur_depth][mul];
    bool Y1 = from_PB_.next_bit_R_recv[cur_depth][mul];

    from_PB_other.next_bit_L_recv[cur_depth][mul] = X0;
    from_PB_other.next_bit_R_recv[cur_depth][mul] = Y0;
    
    bool c = from_P2_to_PB.c_bit[cur_depth][mul];
    
    bool other_blinded = from_PB_.next_bit_L_recv[cur_depth][mul];
 
    bool gamma = xor_if(c, other_blinded, d0);
   
    result = xor_if(gamma, val_mX, (val_bool ^ Y1));    
  }



// i == n
void Verifier::multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, block_t val_mX, bool val_bool, size_t cur_depth, size_t mul, bool party, block_t& result)
{   

   block_t D0 = PBrand.next_block();   
   block_t X0 = val_mX ^ D0;   
   
   bool d0 = PBrand.next_bool();
   bool Y0 = val_bool ^ d0;   
  
   block_t X1;
   bool Y1;
   block_t gamma;
  
   X1 = from_PB_.blinds_recv[cur_depth][mul]; 
   Y1 = from_PB_.bit_blinds_recv[cur_depth][mul]; 

   block_t c = from_P2_to_PB.c[cur_depth][mul];
   
   block_t other_blinded = from_PB_.blinds_recv[cur_depth][mul];


   gamma = xor_if(c, other_blinded, d0);
    
   from_PB_other.blinds_recv[cur_depth][mul]          = X0; // generating the other transcript 
   from_PB_other.bit_blinds_recv[cur_depth][mul]      = Y0; // generating the other transcript
 
   result = xor_if(gamma, val_mX, (val_bool ^ Y1));
}

template<typename prgkey_t>
void Verifier::conditional_swap_and_next_seed(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, block_t s0_L, block_t s0_R, block_t s1_L, block_t s1_R, prgkey_t & key, bool party, size_t cur_depth)
{

  block_t bs0_L;

  multiply_mpc(from_P2_to_PB, from_PB_,  from_PB_other ,s0_L, Pdirection[cur_depth],  cur_depth, 0, party, bs0_L);

  block_t notbs0_L = s0_L ^ bs0_L;
   
  block_t bs1_L;
 
  multiply_mpc(from_P2_to_PB, from_PB_, from_PB_other , s1_L, Pdirection[cur_depth], cur_depth , 1, party, bs1_L); 

  block_t notbs1_L = s1_L ^ bs1_L;
  
  block_t bs0_R;
 
  multiply_mpc(from_P2_to_PB, from_PB_, from_PB_other , s0_R, Pdirection[cur_depth],   cur_depth , 2, party, bs0_R);

  block_t bNotR0_share = s0_R ^ bs0_R;

  block_t bs1_R;

  multiply_mpc(from_P2_to_PB, from_PB_,from_PB_other , s1_R, Pdirection[cur_depth],  cur_depth , 3, party, bs1_R);
  
  block_t bNotR1_share = bs1_R ^ s1_R;

  if(!party)
  {
  seed0[cur_depth+1]  = notbs0_L  ^ bs0_R;
  seed1[cur_depth+1]  = notbs1_L ^ bs1_R;
  }

  if(party)
  {
  seed0[cur_depth+1] = notbs0_L  ^ bs0_R;
  seed1[cur_depth+1] = notbs1_L  ^ bs1_R;
  }

 
}

template<typename leaf_t, typename node_t, typename prgkey_t>
void Verifier::root_layer(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, const LowMC<__m256i> & key,  dpf_key<leaf_t, node_t, prgkey_t> dpfkey, bool party)
{

   block_t seed = dpfkey.root;
   size_t cur_depth = 0;
   const block_t seedL = clear_lsb(seed);
   const block_t seedR = set_lsb(seed);
   
   auto & cw = dpfkey.cw[cur_depth];
   uint8_t cw_t[2] = { get_lsb(cw, 0b01), get_lsb(cw, 0b10) };

   

   block_t CW = dpfkey.cw[cur_depth];
   bool b = get_lsb(dpfkey.root);
   
   block_t s[2];
   uint8_t  t[2];
   
   block_t child[2];

   expand(key, seed, child, t);
   
   s[L] = clear_lsb(xor_if(child[L], CW, !b), 0b11);
   s[R] = clear_lsb(xor_if(child[R], CW, !b), 0b11);
  
   block_t s_L_myshare, s_L_othershare, s_R_myshare, s_R_othershare;

   block_t L_shares[2];    
  
   L_shares[0] = PBrand.next_block();
   L_shares[1] = s[L] ^ L_shares[0];

   if(!party) s_L_myshare = L_shares[0];
   if(party ) s_L_myshare = L_shares[1];

   block_t R_shares[2]; 
  
   R_shares[0] = PBrand.next_block();
   R_shares[1] = s[R] ^ R_shares[0];
  
   if(!party) s_R_myshare = R_shares[0];
   if(party ) s_R_myshare = R_shares[1];
  
   if(!party) from_PB_other.L_shares_recv  =  s[L] ^ L_shares[0];  
   if(!party) from_PB_other.R_shares_recv  =  s[R] ^ R_shares[0];  

   if( party) from_PB_other.L_shares_recv  =  L_shares[0];  
   if( party) from_PB_other.R_shares_recv  =  R_shares[0];

   s_L_othershare = from_PB_.L_shares_recv;    
   s_R_othershare = from_PB_.R_shares_recv;   

   if(!party) conditional_swap_and_next_seed(from_P2_to_PB, from_PB_, from_PB_other , s_L_myshare, s_R_myshare, s_L_othershare, s_R_othershare, key, party, cur_depth); 
   if( party) conditional_swap_and_next_seed(from_P2_to_PB, from_PB_, from_PB_other , s_L_othershare, s_R_othershare, s_L_myshare, s_R_myshare, key, party, cur_depth);
   
   bool bitL_shares[2];  
   bool bitR_shares[2]; 
 
   bool bit_L_myshare, bit_L_othershare, bit_R_myshare, bit_R_othershare;

   bool bit_L0_share, bit_R0_share, bit_L1_share, bit_R1_share;
   
   bitL_shares[0] = PBrand.next_bool();
   bitR_shares[0] = PBrand.next_bool();
 

   bitL_shares[1] = bitL_shares[0] ^ t[L] ^ cw_t[L] & (b);
   bitR_shares[1] = bitR_shares[0] ^ t[R] ^ cw_t[R] & (b);
   

   if(!party)
   {
    bit_L_myshare = bitL_shares[0];
    bit_R_myshare = bitR_shares[0];

    from_PB_other.bit_L_shares_recv = bitL_shares[1];
    from_PB_other.bit_R_shares_recv = bitR_shares[1];
   }

   if(party)
   {
    bit_L_myshare = bitL_shares[1];
    bit_R_myshare = bitR_shares[1];

    from_PB_other.bit_L_shares_recv = bitL_shares[0];
    from_PB_other.bit_R_shares_recv = bitR_shares[0];
   }
 
   bit_L_othershare =  from_PB_.bit_L_shares_recv;
   bit_R_othershare =  from_PB_.bit_R_shares_recv; 
   
   
   if(!party) get_next_bits(from_P2_to_PB, from_PB_, from_PB_other , bit_L_myshare, bit_R_myshare, bit_L_othershare, bit_R_othershare, party, cur_depth);
   if( party) get_next_bits(from_P2_to_PB, from_PB_, from_PB_other , bit_L_othershare, bit_R_othershare, bit_L_myshare, bit_R_myshare, party, cur_depth);
}

template<typename leaf_t, typename node_t, typename prgkey_t>
void Verifier::middle_layers(const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other, const LowMC<__m256i> & key, const block_t seed0, const block_t seed1, dpf_key<leaf_t, node_t, prgkey_t> dpfkey, size_t cur_depth, bool party)
{
  
  auto & cw = dpfkey.cw[cur_depth];
  uint8_t cw_t[2] = { get_lsb(cw, 0b01), get_lsb(cw, 0b10) };



  block_t blind[rounds];
  block_t gamma[rounds];
  
  block_t P0gamma[2][rounds];
  block_t P0blinds0[rounds]; 

  Gen_Blinds(P0blinds0);
  
  for(size_t j = 0; j < rounds; ++j)
  {
    P0gamma[0][j] = from_P2_to_PB.gamma0[cur_depth][j];  
    P0gamma[1][j] = from_P2_to_PB.gamma0[cur_depth][j];  
  } 
  
  block_t s0_L,  s0_R; 
  uint8_t  t0_L, t0_R;
  
  if(!party) expand_mpc_verify(key, seed0, from_P2_to_PB, from_PB_, from_PB_other , s0_L, s0_R, t0_L, t0_R, P0blinds0, P0gamma[0], party, cur_depth, false);  
  if(party ) expand_mpc_verify(key, seed0, from_P2_to_PB, from_PB_, from_PB_other , s0_L, s0_R, t0_L, t0_R, P0blinds0, P0gamma[1], party, cur_depth, false);
  
   block_t L0_shares_2, R0_shares_2;
  
  if(!party)
  {
   L0_shares_2 = clear_lsb(xor_if(s0_L, cw, bit0[cur_depth]), 0b11);
   R0_shares_2 = clear_lsb(xor_if(s0_R, cw, bit0[cur_depth]), 0b11);
  }
 
  if(party)
  {
   L0_shares_2 = clear_lsb(xor_if(s0_L, cw, !bit0[cur_depth]), 0b11);
   R0_shares_2 = clear_lsb(xor_if(s0_R, cw, !bit0[cur_depth]), 0b11);
  }

  block_t P1gamma[2][rounds];
  block_t P1blinds0[rounds];
  block_t P1blinds1[rounds];

  Gen_Blinds(P1blinds0); 

  for(size_t j = 0; j < rounds; ++j)
  {
  P1gamma[0][j] = from_P2_to_PB.gamma1[cur_depth][j];  
  P1gamma[1][j] = from_P2_to_PB.gamma1[cur_depth][j];  
  } 


  block_t s1_L,  s1_R; 
  uint8_t t1_L, t1_R;
  if(!party) expand_mpc_verify(key, seed1,  from_P2_to_PB,  from_PB_,from_PB_other ,  s1_L, s1_R, t1_L, t1_R,  P1blinds0, P1gamma[0], party, cur_depth, true);
  if(party)  expand_mpc_verify(key, seed1,  from_P2_to_PB,  from_PB_,from_PB_other ,  s1_L, s1_R, t1_L, t1_R,  P1blinds0, P1gamma[1], party, cur_depth, true);

  block_t L0_shares, R0_shares;
  block_t L1_shares, R1_shares;
  
  block_t L1_shares_2, R1_shares_2;

  if(!party)
  {
   L1_shares_2 = clear_lsb(xor_if(s1_L, cw, bit1[cur_depth]), 0b11);
   R1_shares_2 = clear_lsb(xor_if(s1_R, cw, bit1[cur_depth]), 0b11);
  }

  if(party)
  {
   L1_shares_2 = clear_lsb(xor_if(s1_L, cw, !bit1[cur_depth]), 0b11);
   R1_shares_2 = clear_lsb(xor_if(s1_R, cw, !bit1[cur_depth]), 0b11);
  }


   if(!party) conditional_swap_and_next_seed(from_P2_to_PB, from_PB_,from_PB_other , L0_shares_2, R0_shares_2, L1_shares_2, R1_shares_2,  key, party, cur_depth);
   if( party) conditional_swap_and_next_seed(from_P2_to_PB, from_PB_,from_PB_other , L0_shares_2, R0_shares_2, L1_shares_2, R1_shares_2,  key, party, cur_depth);

   bool bitL_shares, bitR_shares;
  
   bool bit0_L, bit0_R, bit1_L, bit1_R;
  
   bool my_val = cw_t[L];

   if(party) my_val = bit0[cur_depth];
  
   multiply_mpc(from_P2_to_PB, from_PB_,from_PB_other , my_val,  party, 0, bit0_L, cur_depth);

   if(!party) bit0_L = bit0_L ^ t0_L ^ cw_t[L]  &  (bit0[cur_depth]);
  
  if(party)
  {
  
    bit0_L = bit0_L ^ t0_L;
  } 
  
  my_val = cw_t[R];


  if(party) my_val = bit0[cur_depth];

  multiply_mpc(from_P2_to_PB, from_PB_,from_PB_other , my_val, party, 1, bit0_R, cur_depth);

  if(!party) bit0_R = bit0_R ^ t0_R ^ (cw_t[R] & bit0[cur_depth]);

  if(party)
  {
    bit0_R = bit0_R ^ t0_R;
  } 
  
  my_val = bit1[cur_depth];
 
 
  if(party) my_val = cw_t[L];

  multiply_mpc(from_P2_to_PB, from_PB_,from_PB_other , my_val,  party, 2, bit1_L, cur_depth);
 
   if( party)  bit1_L = bit1_L ^  t1_L ^ (cw_t[L] & bit1[cur_depth]);

  if(!party) bit1_L = bit1_L ^ t1_L;
 
  if(party) my_val = cw_t[R];

  multiply_mpc(from_P2_to_PB, from_PB_,from_PB_other , my_val,  party, 3, bit1_R, cur_depth);

  if(!party) bit1_R = bit1_R ^ t1_R;
 
  if(party)  bit1_R = bit1_R ^ (t1_R) ^ (cw_t[R] & bit1[cur_depth]);
  
  get_next_bits(from_P2_to_PB, from_PB_, from_PB_other , bit0_L, bit0_R, bit1_L, bit1_R,  party, cur_depth);
   
 }



// void Verifier::multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_,const bool& x0, bool b0, const bool& x1, bool b1, bool results[2])
  // {
  //   int next_b;

  //   bool D0 = rand(); 
  //   bool D1 = rand(); 
     
  //   bool d0 = rand();     
  //   bool d1 = rand();
      
  //   const bool alpha = rand(); 

  //   bool c0 = xor_if(alpha, D0, d1);
  //   bool c1 = xor_if(alpha, D1, d0);

  //   bool DD0 = x0 ^ D0;
  //   bool dd0 = b0 ^ d0;

  //   bool DD1 = x1 ^ D1;
  //   bool dd1 = b1 ^ d1;

    

  //  const bool gamma0 = xor_if(c0, DD1, d0);
  
  //  bool xx0 = x0;   
   
  //  xor_if(gamma0, xx0, (b0 ^ dd1)); 

  //  const bool gamma1 = xor_if(c1, DD0, d1);
  //  bool xx1 = x1;  
  //  xor_if(gamma1, xx1, (b1 ^ dd0)); 

  //  results[0] = xor_if(gamma0, xx0, (b0 ^ dd1));
  
  //  results[1] = xor_if(gamma1, xx1, (b1 ^ dd0)); 

  // }


//   void Verifier::multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_, const bool& x0, bool b0, const bool& x1, bool b1, bool results[2], size_t mul, size_t cur_depth)
//   {
//     // bool D0 = PBrand.next_bool(); 
//     // bool D1 = PBrand_other.next_bool(); 
     
//     // bool d0 = PBrand.next_bool();     
//     // bool d1 = PBrand_other.next_bool();
      
//     // const bool alpha = P2rand.next_bool(); 

//     // bool c0 = xor_if(alpha, D0, d1);
//     // bool c1 = xor_if(alpha, D1, d0);

//     // P0_message.PB_middle.next_bit_L_sent[cur_depth][mul] = x0 ^ D0;
//     // P0_message.PB_middle.next_bit_R_sent[cur_depth][mul] = b0 ^ d0;

//     // P1_message.PB_middle.next_bit_L_sent[cur_depth][mul] = x1 ^ D1;
//     // P1_message.PB_middle.next_bit_R_sent[cur_depth][mul] = b1 ^ d1;


//     // const bool gamma0 = xor_if(c0, P1_message.PB_middle.next_bit_L_sent[cur_depth][mul], d0);
  
//     // bool xx0 = x0;   
//     // xor_if(gamma0, xx0, (b0 ^ P1_message.PB_middle.next_bit_R_sent[cur_depth][mul])); 

//     // const bool gamma1 = xor_if(c1,P0_message.PB_middle.next_bit_L_sent[cur_depth][mul], d1);
//     // bool xx1 = x1;  
//     // xor_if(gamma1, xx1, (b1 ^ P0_message.PB_middle.next_bit_R_sent[cur_depth][mul])); 

//     // results[0] = xor_if(gamma0, xx0, (b0 ^ P1_message.PB_middle.next_bit_R_sent[cur_depth][mul]));
  
//     // results[1] = xor_if(gamma1, xx1, (b1 ^ P0_message.PB_middle.next_bit_R_sent[cur_depth][mul])); 

//    // if(cur_depth == 0)
//    //  {
      
//    //    P0_message.PB_root.gamma_bit[mul] = gamma0; 
//    //    P1_message.PB_root.gamma_bit[mul] = gamma1; 

//    //    P0_message.PB_root.next_bit_L_sent[mul] = P0_message.PB_middle.next_bit_L_sent[cur_depth][mul];
//    //    P0_message.PB_root.next_bit_L_recv[mul] = P1_message.PB_middle.next_bit_L_sent[cur_depth][mul];
//    //    P1_message.PB_root.next_bit_L_recv[mul] = P0_message.PB_middle.next_bit_L_sent[cur_depth][mul];
//    //    P1_message.PB_root.next_bit_R_sent[mul] = P1_message.PB_middle.next_bit_R_sent[cur_depth][mul];
//    //    P0_message.PB_root.next_bit_R_recv[mul] = P1_message.PB_middle.next_bit_R_sent[cur_depth][mul];
//    //    P1_message.PB_root.next_bit_R_recv[mul] = P0_message.PB_middle.next_bit_R_sent[cur_depth][mul];
//    //  }
//    //  else
//    //  {
//    //    P1_message.PB_middle.next_bit_L_recv[cur_depth][mul] = P0_message.PB_middle.next_bit_L_sent[cur_depth][mul];
//    //    P1_message.PB_middle.next_bit_R_recv[cur_depth][mul] = P0_message.PB_middle.next_bit_R_sent[cur_depth][mul];
//    //    P0_message.PB_middle.next_bit_L_recv[cur_depth][mul] = P1_message.PB_middle.next_bit_L_sent[cur_depth][mul];
//    //    P0_message.PB_middle.next_bit_R_recv[cur_depth][mul] = P1_message.PB_middle.next_bit_R_sent[cur_depth][mul];
//    //    P0_message.PB_middle.gamma_bit[cur_depth][mul]       = gamma0;
//    //    P1_message.PB_middle.gamma_bit[cur_depth][mul]       = gamma1;
//    //  }

//   }


  // template<typename T>
  // void Verifier::multiply_mpc(const from_P2 from_P2_to_PB, const from_PB from_PB_, const T& x0, bool b0, const T& x1, bool b1, T results[2], size_t mul, size_t cur_depth)
  // {
   
  //   T D0 = PBrand.next_block();
    
  //   T D1 = PBrand_other.next_block();
    
  //   bool d0 = PBrand.next_bool(); 
    
  //   bool d1 = PBrand_other.next_bool(); 

  //   T rand;
   
  //   P2rand.next_bool(); 

  //   const T alpha = rand;

  //   T c0 = xor_if(alpha, D0, d1);
  //   T c1 = xor_if(alpha, D1, d0);
    
  //   T x0_blind    = x0 ^ D0;

  //   bool b0_blind = b0 ^ d0;
    
  //   T x1_blind    = x1 ^ D1;
    
  //   bool b1_blind = b1 ^ d1;

  //   const T gamma0 = xor_if(c0, x1_blind, d0);

  //   T xx0 = x0;

  //   const T gamma1 = xor_if(c1, x0_blind, d1);
    
  //   T xx1 = x1;

  //   results[0] = xor_if(gamma0, xx0, (b0 ^ b1_blind));
    
  //   results[1] = xor_if(gamma1, xx1, (b1 ^ b0_blind)); 

  // }

  // void Verifier::Gen_Blinds(block_t blinds0[rounds], block_t blinds1[rounds], block_t gamma[2][rounds])
  // {
  //   block_t rand[rounds];
  //   for (unsigned r = 0; r < rounds; ++r)
  //    {        
  //       blinds0[r] = PBrand.next_block();
  //       blinds1[r] = PBrand_other.next_block();
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




//  template<typename KEY_TYPE, typename __mX>
//  inline void Verifier::expand_mpc(KEY_TYPE & key, const blocks<__mX> & seed0, const blocks<__mX> & seed1,  blocks<__mX> s0[2], bool t0[2], 
//                         blocks<__mX> s1[2], bool t1[2], const block blind0[rounds], const block blind1[rounds], block gamma[2][rounds], size_t cur_depth, bool LR)
//  {
  
//    std::vector<block> c0(2 + rounds + rounds);
//    std::vector<block> c1(2 + rounds + rounds);
  

//    const blocks<__mX> seed0L = clear_lsb(seed0); // _mm_clearlsb_si128(seed);
//    const blocks<__mX> seed0R = set_lsb(seed0); //  _mm_setlsb_si128(seed);

//    const blocks<__mX> seed1L = clear_lsb(seed1); // _mm_clearlsb_si128(seed);
//    const blocks<__mX> seed1R = clear_lsb(seed1); //  _mm_setlsb_si128(seed);

//    s0[L] = seed0L;
//    s0[R] = seed0R;
 
//    s1[L] = seed1L;
//    s1[R] = seed1R;

   
//    std::pair<std::vector<block>,std::vector<block>> encrypt_outL = key.encrypt_MPC_proof(seed0L, seed1L, blind0, blind1, gamma);
      
//    for(size_t j = 0; j < rounds; ++j)
//     {
//       c0[2 + 2*j]     = encrypt_outL.first[j+1];
//       c0[2 + 2*j + 1] = encrypt_outL.second[j+1];
//     }

//    std::pair<std::vector<block>,std::vector<block>> encrypt_outR = key.encrypt_MPC_proof(seed0R, seed1R, blind0, blind1, gamma);
    
 
//    for(size_t j = 0; j < rounds; ++j)
//    {
//      c1[2 + 2*j]     = encrypt_outR.first[j+1];
//      c1[2 + 2*j + 1] = encrypt_outR.second[j+1];
//    }

//    s0[L] = encrypt_outL.first[rounds];
//    s0[R] = encrypt_outR.first[rounds];

//    s1[L] = encrypt_outL.second[rounds];
//    s1[R] = encrypt_outR.second[rounds];

//    s0[L] ^= seed0L; 
//    s1[L] ^= seed1L;
    
//    t0[L] = get_lsb(s0[L]); 
//    t1[L] = get_lsb(s1[L]);
     
//    s0[L] = clear_lsb(s0[L]); 
//    s1[L] = clear_lsb(s1[L]);

//    s0[R] ^= seed0R; s1[R] ^= seed1R;
//    t0[R] = get_lsb(s0[R]); t1[R] = get_lsb(s1[R]);
//    s0[R] = clear_lsb(s0[R]); s1[R] = clear_lsb(s1[R]);   

//    c0[L] = s0[L];
//    c1[L] = s1[L];

//    c0[R] = s0[R];
//    c1[R] = s1[R];


//    std::pair<std::pair<std::vector<block>,std::vector<block>>, std::pair<std::vector<block>,std::vector<block>>> encrypt_out = std::pair(encrypt_outL, encrypt_outR);
  
//  }








//   void Verifier::conditional_swap_and_next_seed(const from_P2 from_P2_to_PB, const from_PB from_PB_ , LowMC& key, block s0_L[2], block s0_R[2], block s1_L[2], block s1_R[2], 
//                                  block left_out[2], block right_out[2], size_t cur_depth)
//   {

//      block bs0_L[2];
     
//      multiply_mpc(from_P2_to_PB, from_PB_, s0_L[0], Pdirection[cur_depth], s0_L[1], Pdirection_other[cur_depth], bs0_L, 0, cur_depth); // L_0 * direction

//      block notbs0_L[2];
//      notbs0_L[0] = s0_L[0] ^ bs0_L[0];
//      notbs0_L[1] = s0_L[1] ^ bs0_L[1];
     
//      block bs1_L[2];
//      multiply_mpc(from_P2_to_PB, from_PB_, s1_L[0], Pdirection[cur_depth], s1_L[1], Pdirection_other[cur_depth], bs1_L, 1, cur_depth); // L_1 * direction

//      block notbs1_L[2];
//      notbs1_L[0] = s1_L[0] ^ bs1_L[0];
//      notbs1_L[1] = s1_L[1] ^ bs1_L[1];
     
     
//      block bs0_R[2];
//      multiply_mpc(from_P2_to_PB, from_PB_, s0_R[0], Pdirection[cur_depth], s0_R[1], Pdirection_other[cur_depth], bs0_R, 2, cur_depth); // R_0 * direction
 
//      block notbs0_R[2];
//      notbs0_R[0] = s0_R[0] ^ bs0_R[0];
//      notbs0_R[1] = s0_R[1] ^ bs0_R[1];

    
//      block bs1_R[2];    
//      multiply_mpc(from_P2_to_PB, from_PB_, s1_R[0], Pdirection[cur_depth], s1_R[1], Pdirection_other[cur_depth], bs1_R, 3, cur_depth); // R_1 * direction
   
//      block notbs1_R[2]; 
//      notbs1_R[0] = s1_R[0] ^ bs1_R[0];
//      notbs1_R[1] = s1_R[1] ^ bs1_R[1];

//      right_out[0] = notbs0_L[0] ^ bs0_R[0] ^ notbs0_L[1] ^ bs0_R[1];
//      right_out[1] = notbs1_L[0] ^ bs1_R[0] ^ notbs1_L[1] ^ bs1_R[1];

//      left_out[0] = bs0_L[0] ^ notbs0_R[0] ^ bs1_L[0] ^ notbs1_R[0];
//      left_out[1] = bs0_L[1] ^ notbs0_R[1] ^ bs1_L[1] ^ notbs1_R[1];
    
//      seed0_[cur_depth+1][0] = notbs0_L[0] ^ bs0_R[0];
//      seed0_[cur_depth+1][1] = notbs0_L[1] ^ bs0_R[1];
//      seed1_[cur_depth+1][0] = notbs1_L[0] ^ bs1_R[0];
//      seed1_[cur_depth+1][1] = notbs1_L[1] ^ bs1_R[1];
 
//   }

//  template<typename KEY_TYPE, typename __mX>
//  inline void Verifier::root_layer(KEY_TYPE& key,  dpf_key<__mX, nitems> dpfkey, dpf_key<__mX, nitems> dpfkey_other, bool party)
//  {

//     // std::cout << "root_layer: " << std::endl;
//     // block seed0 = dpfkey.root; 
//     // block seed1 = dpfkey_other.root;

//     // size_t cur_depth = 0;
   
//     // block CW = dpfkey.cw[cur_depth];
   
//     // bool b0 = get_lsb(dpfkey.root);
    
//     // bool b1 = get_lsb(dpfkey_other.root);

//     // block s0[2], s1[2];    bool t0[2], t1[2];

//     // block child[2];

//     // expand(key, seed0, child, t0);
   
//     // s0[L] = xor_if(child[L], CW, !b0);
//     // s0[R] = xor_if(child[R], CW, !b0);

//     // expand(key, seed1, child, t1);
   
//     // s1[L] = xor_if(child[L], CW, !b1);
//     // s1[R] = xor_if(child[R], CW, !b1);

//     // block s0_L[2]; block s0_R[2];
   
//     // PBrand.next(s0_L[0]); 
    
//     // s0_L[1] = s0_L[0] ^ s0[L];

//     // PBrand.next(s0_R[0]);
    
//     // s0_R[1] = s0_R[0] ^ s0[R];
 
//     // block s1_L[2]; block s1_R[2];
   
//     // PBrand_other.next(s1_L[0]);
    
//     // s1_L[1] = s1_L[0] ^ s1[L];
    
//     // PBrand_other.next(s1_R[0]); 

//     // s1_R[1] = s1_R[0] ^ s1[R];

//     // block left_out[2], right_out[2];

//     // conditional_swap_and_next_seed(key, s0_L, s0_R, s1_L, s1_R, left_out, right_out, cur_depth);
//     // std::cout << "reconstructed L = " << (left_out[0]  ^ left_out[1]) << std::endl;
   
//     // bool bit0_L[2]; 

//     // bit0_L[0] = PBrand.next_bool(); 
//     // bit0_L[1] =  bit0_L[0] ^ (t0[L] ^ dpfkey.t[cur_depth][L] & (b0));
    
//     // bool bit0_R[2]; 
//     // bit0_R[0] = PBrand.next_bool(); 
//     // bit0_R[1] = bit0_R[0] ^ (t0[R] ^ dpfkey.t[cur_depth][R] & (b0));
    
//     // bool bit1_L[2];
//     // bit1_L[0] = PBrand_other.next_bool(); 
//     // bit1_L[1] = bit1_L[0] ^ (t1[L] ^ dpfkey_other.t[cur_depth][L] & (b1));
    
//     // bool bit1_R[2];
//     // bit1_R[0] = PBrand_other.next_bool(); 
//     // bit1_R[1] = bit1_R[0] ^ (t1[R] ^ dpfkey_other.t[cur_depth][R] & (b1));

//  }







//   void Verifier::get_next_bits(const from_P2 from_P2_to_PB, const from_PB from_PB_, bool L0_shares[2], bool R0_shares[2], bool L1_shares[2], bool R1_shares[2], size_t cur_depth)
//   {
    
//     bool notbs0_L[2];    

//     multiply_mpc(from_P2_to_PB, from_PB_, L0_shares[0], Pdirection[cur_depth], L0_shares[1], !Pdirection_other[cur_depth], notbs0_L,  0, cur_depth); // L_0 * (direction - 1)

//     bool notbs1_L[2];
    
//     multiply_mpc(from_P2_to_PB, from_PB_, L1_shares[0], Pdirection[cur_depth], L1_shares[1], !Pdirection_other[cur_depth], notbs1_L,  1, cur_depth); // L_1 * (direction - 1)

//     bool bs0_R[2];
    
//     multiply_mpc(from_P2_to_PB, from_PB_, R0_shares[0], Pdirection[cur_depth], R0_shares[1], Pdirection_other[cur_depth],  bs0_R, 2  , cur_depth);

//     bool bs1_R[2];

//     multiply_mpc(from_P2_to_PB, from_PB_, R1_shares[0], Pdirection[cur_depth], R1_shares[1], Pdirection_other[cur_depth], bs1_R, 3  , cur_depth); // R_1 * direction    
    
    
//     bit0_next[cur_depth+1][0] = notbs0_L[0] ^ bs0_R[0];
    
//     bit0_next[cur_depth+1][1] = notbs0_L[1] ^ bs0_R[1];    
    
//     bit1_next[cur_depth+1][0] = notbs1_L[0] ^ bs1_R[0];    
    
//     bit1_next[cur_depth+1][1] = notbs1_L[1] ^ bs1_R[1];

//   }


//   //void Verifier::middle_layers(const from_P2 from_P2_to_PB, const from_PB from_PB_, LowMC& key, dpf_key<__mX, nitems> dpfkey, dpf_key<__mX, nitems> dpfkey_other , const size_t cur_depth) 
//  // {
    
//  //    std::cout << "middle layer: " << cur_depth << std::endl;
//  //    block CW = dpfkey.cw[cur_depth];
     

//  //    block P0gamma[2][rounds];
//  //    block P0blinds0[rounds];
//  //    block P0blinds1[rounds];

//  //    Gen_Blinds(P0blinds0, P0blinds1, P0gamma);
     
    
//  //    block s0[2], s1[2];
//  //    bool  t0[2], t1[2]; 



//  //    expand_mpc(key, seed0_[cur_depth][0], seed0_[cur_depth][1],  s0, t0, s1, t1, P0blinds0, P0blinds1, P0gamma, cur_depth, false);



//  //    block L0_shares[2]; 
//  //    block R0_shares[2];
  
//  //    L0_shares[0] = xor_if(s0[L], CW, bit0_next[cur_depth][0]);
//  //    L0_shares[1] = xor_if(s1[L], CW, !bit0_next[cur_depth][1]);
//  //    R0_shares[0] = xor_if(s0[R], CW, bit0_next[cur_depth][0]);
//  //    R0_shares[1] = xor_if(s1[R], CW, !bit0_next[cur_depth][1]);


//  //    block ss0[2], ss1[2];

//  //    bool tt0[2], tt1[2]; 

//  //    block P1gamma[2][rounds];
//  //    block P1blinds0[rounds];
//  //    block P1blinds1[rounds];

//  //    Gen_Blinds(P1blinds0, P1blinds1, P1gamma);
    
//  // //    for(size_t j = 0; j < rounds; ++j)
//  // //    {
//  // //    P0_message.PB_middle.gamma1[cur_depth][j] = P1gamma[0][j];
//  // //    P1_message.PB_middle.gamma1[cur_depth][j] = P1gamma[1][j];
//  // //    } 

   
//  //     expand_mpc(key, seed1_[cur_depth][0], seed1_[cur_depth][1],  ss0, tt0, ss1, tt1, P1blinds0, P1blinds1, P1gamma, cur_depth, true);


//  //    block L1_shares[2];
//  //    block R1_shares[2];

//  //    L1_shares[0] = xor_if(ss0[L], CW, bit1_next[cur_depth][0]);
//  //    L1_shares[1] = xor_if(ss1[L], CW, !bit1_next[cur_depth][1]);
//  //    R1_shares[0] = xor_if(ss0[R], CW, bit1_next[cur_depth][0]);
//  //    R1_shares[1] = xor_if(ss1[R], CW, !bit1_next[cur_depth][1]);

 


//  //    block left_out[2]; 
//  //    block right_out[2];
 
    
//  //}

// inline void Verifier::verify_P2()
// {

// }

//template<typename leaf_t, typename node_t, typename prgkey_t>
inline void Verifier::expand_mpc_verify(const LowMC<__m256i> & key, const block_t seed, const from_P2 from_P2_to_PB, const from_PB from_PB_, from_PB & from_PB_other,  block_t &s0_L, block_t & s0_R, 
                                        uint8_t & t0_L, uint8_t & t0_R,  block_t blind[], block_t gamma[], bool party, size_t cur_depth, bool LR)
{

   const block_t seedL = clear_lsb(seed); // _mm_clearlsb_si128(seed);
   
   block_t seedR ;
   
   if(party)  seedR = clear_lsb(seed); //  _mm_setlsb_si128(seed);
   
   if(!party) seedR = set_lsb(seed); //  _mm_setlsb_si128(seed);
  
 
   std::vector<block_t> c2L;
   
   if(!LR) c2L = from_PB_.seed0L_encrypt[cur_depth];   
   if(LR) c2L =  from_PB_.seed1L_encrypt[cur_depth];   
 
  
 
   std::vector<block_t> outL  = key.encrypt_MPC_verify(seedL, c2L, blind, gamma, party);

   if(!LR) from_PB_other.seed0L_encrypt[cur_depth] = outL;
   if( LR) from_PB_other.seed1L_encrypt[cur_depth] = outL;

   // for(size_t j = 1; j < rounds; ++j)
   // {    
   //  if(!LR) 
   //  {
      
   //    std::cout << "outL[" << j << "] = " << ( outL[j] ^ from_PB_.seed0L_encrypt_other[cur_depth][j]) << std::endl;
   //  }
   //  if(LR) 
   //  { 
   //   //from_PB_other.seed1L_encrypt[j] = outL[j];
   //   std::cout << "outL[" << j << "] = " << ( outL[j] ^ from_PB_.seed1L_encrypt_other[cur_depth][j]) << std::endl;
   //  }
   // }

   std::vector<block_t> c2R; 
   
   if(!LR) 
    {

      c2R = from_PB_.seed0R_encrypt[cur_depth];  
    }

   if(LR) 
    {

      c2R =  from_PB_.seed1R_encrypt[cur_depth];  
    }

   std::vector<block_t> outR = key.encrypt_MPC_verify(seedR, c2R, blind, gamma, party);
  
   //from_PB_other.seed0R_encrypt[cur_depth] = outR;
  


   if(!LR) from_PB_other.seed0R_encrypt[cur_depth] = outR;
   if( LR) from_PB_other.seed1R_encrypt[cur_depth] = outR;
   
   // for(size_t j = 0; j <= rounds; ++j)
   // {    
   //  if(!LR)
   //  {
   //   //from_PB_other.seed0R_encrypt[j] = outR[j];  
   //   std::cout << "outR[" << j << "] = " << (outR[j] ^ from_PB_.seed0R_encrypt_other[cur_depth][j])<< std::endl;
   //  } 
   //  if(LR)
   //  {
   //   //from_PB_other.seed1R_encrypt[j] = outR[j];
   //   std::cout << "outR[" << j << "] = " <<  (outR[j] ^ from_PB_.seed1R_encrypt_other[cur_depth][j])<< std::endl;
   //  }
   // }

    s0_L = outL[rounds];
    s0_R = outR[rounds];

  
    s0_L ^= seedL; 
     
    t0_L = get_lsb(s0_L); 

    s0_L = clear_lsb(s0_L); 
    
    s0_R ^= seedR; 
    
    t0_R = get_lsb(s0_R); 
    
    s0_R = clear_lsb(s0_R);     
}
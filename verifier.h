class Verifier
{

  public:
const size_t depth;
  std::vector<bool> Pdirection;
   
  std::vector<block_t> seed0;
  std::vector<block_t> seed1;

  std::vector<bool> bit0;
  std::vector<bool> bit1;
   
  Verifier(AES_KEY& aeskey, __m128i seed, size_t len, size_t depth_)
      : depth(depth_), PBrand(aeskey, seed, len) 
  {  
   Pdirection.resize(depth_ +1); 
   seed0.resize(depth_ +1);
   seed1.resize(depth_+1);
   bit0.resize(depth_ +1);
   bit1.resize(depth_ +1);
  }

   


  // ///// Functions for dice role i = n
   void get_next_bits(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], bool L0_shares, bool R0_shares, bool L1_shares, bool R1_shares,  bool party ,size_t cur_depth);
  
   void Gen_Blinds(block_t blinds0[]);
  

   template<typename prgkey_t>
   auto prg_mpc_verify(const prgkey_t & key, const block_t seed, const from_P2 from_P2_to_PB, 
                               const from_PB_ received_view[], from_PB_ generated_view[], block_t blind[], block_t gamma[], size_t cur_depth ,bool party, bool LR , size_t len);

   template<  typename prgkey_t>
   void expand_mpc_verify(const prgkey_t& key, const block_t seed, const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], block_t &s0_L, block_t & s0_R, 
                                         uint8_t & t0_L, uint8_t & t0_R, block_t blind[], block_t gamma[], bool party, size_t cur_depth, bool LR);
  
   template<typename leaf_t, typename node_t, typename prgkey_t>
   void root_layer(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[],  const prgkey_t & key,  dpf_key<leaf_t, node_t, prgkey_t> dpfkey,  bool party);
  

   void multiply_mpc(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], block_t val_mX, bool val_bool, size_t cur_depth, size_t mul, bool party, block_t& out);
  
   void multiply_mpc(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], bool val, bool val_bool, bool blind_mX, bool blind_bool,  bool party, size_t mul , bool& result, size_t cur_depth);
  
   void multiply_mpc(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], bool val, bool party, size_t mul , bool& result, size_t cur_depth);
  
  template<typename prgkey_t>
  void conditional_swap_and_next_seed(const from_P2 from_P2_to_PB,  const from_PB_ received_view[], from_PB_ generated_view[], block_t L_share, block_t R_share, block_t L1_share, 
                                      block_t R1_share, prgkey_t & key,  bool party, size_t cur_depth);
  
  template<typename leaf_t, typename node_t, typename prgkey_t>
  void middle_layers(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], const prgkey_t & key, const block_t seed0, const  block_t seed1, dpf_key<leaf_t, node_t, prgkey_t> dpfkey, size_t cur_depth, bool party);
 

  ~ Verifier()
  {

  }

  private:
  MPCrandomness PBrand;
};

 
  ////////////////////////////////// The below functions for i == n //////////////////////////////////

 void Verifier::Gen_Blinds(block_t blinds0[])
  {
    for (unsigned r = 0; r < rounds; ++r)
     {
        PBrand.next(blinds0[r]); 
     }
  }


  void Verifier::get_next_bits(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], bool L0_shares, bool R0_shares, bool L1_shares, bool R1_shares, bool party, size_t cur_depth)
  {    
    bool blind_L,  blind_R;

    blind_L = received_view[cur_depth].next_bit_L_recv[0];
    blind_R = received_view[cur_depth].next_bit_R_recv[0];
 
    bool bNotL0_share;
   
    if(!party) multiply_mpc(from_P2_to_PB, received_view,generated_view , L0_shares, Pdirection[cur_depth], blind_L, blind_R,   party, 0, bNotL0_share, cur_depth);
    if(party)  multiply_mpc(from_P2_to_PB, received_view,generated_view , L0_shares, !Pdirection[cur_depth], blind_L, blind_R,   party, 0, bNotL0_share, cur_depth);
 
    blind_L = received_view[cur_depth].next_bit_L_recv[1];
    blind_R = received_view[cur_depth].next_bit_R_recv[1];

    bool bNotL1_share;
   
    if(!party) multiply_mpc(from_P2_to_PB, received_view, generated_view , L1_shares, Pdirection[cur_depth] , blind_L, blind_R,   party, 1, bNotL1_share, cur_depth);
    if(party ) multiply_mpc(from_P2_to_PB, received_view, generated_view , L1_shares, !Pdirection[cur_depth] , blind_L, blind_R,  party, 1, bNotL1_share, cur_depth);
    
    blind_L = received_view[cur_depth].next_bit_L_recv[2];
    blind_R = received_view[cur_depth].next_bit_R_recv[2];
 
    bool bR0_share;
    multiply_mpc(from_P2_to_PB, received_view, generated_view , R0_shares, Pdirection[cur_depth], blind_L, blind_R,  party, 2, bR0_share, cur_depth);

    blind_L = received_view[cur_depth].next_bit_L_recv[3];
    blind_R = received_view[cur_depth].next_bit_R_recv[3];
 
    bool bR1_share;
    multiply_mpc(from_P2_to_PB, received_view, generated_view , R1_shares, Pdirection[cur_depth], blind_L, blind_R,   party, 3, bR1_share, cur_depth);

    bit0[cur_depth+1] = bNotL0_share  ^ bR0_share;    
    bit1[cur_depth+1] = bNotL1_share  ^ bR1_share; 
    
  }

//   //multiply_mpc
  void Verifier::multiply_mpc(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], bool val_mX,  bool party, size_t mul , bool& result, size_t cur_depth)
  {

    bool val_bool = false;
     
    bool D0 = PBrand.next_bool();
    bool X0 = val_mX ^ D0;  ;
         
    bool d0 = PBrand.next_bool();
    bool Y0 = val_bool ^ d0;   
    
    bool X1 = received_view[cur_depth].next_bit_L_recv2[mul];
    bool Y1 = received_view[cur_depth].next_bit_R_recv2[mul];
    bool c  = from_P2_to_PB.c_bit2[cur_depth][mul];
 
    generated_view[cur_depth].next_bit_L_recv2[mul] = X0;
    generated_view[cur_depth].next_bit_R_recv2[mul] = Y0;
    
    if(party)
    {
     generated_view[cur_depth].next_bit_L_recv2[mul] = (val_bool ^ D0);
     generated_view[cur_depth].next_bit_R_recv2[mul] = val_mX ^ d0;
    }

    bool other_blinded = received_view[cur_depth].next_bit_L_recv2[mul];
     
    bool gamma = xor_if(c, other_blinded, d0);
      
    if(!party) result = xor_if(gamma, val_mX, (val_bool ^ Y1));
    if( party) result = xor_if(gamma, val_bool, (val_mX ^ X1));
  }
  

  
  void Verifier::multiply_mpc(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], bool val_mX, bool val_bool, bool blind_mX, bool blind_bool,  bool party, size_t mul, bool& result, size_t cur_depth)
  {
 
    bool D0 = PBrand.next_bool();
    bool X0 = val_mX ^ D0; 
 
  
    bool d0 = PBrand.next_bool();
    bool Y0 = val_bool ^ d0;  
    
    
    bool X1 = received_view[cur_depth].next_bit_L_recv[mul];
    bool Y1 = received_view[cur_depth].next_bit_R_recv[mul];

    generated_view[cur_depth].next_bit_L_recv[mul] = X0;
    generated_view[cur_depth].next_bit_R_recv[mul] = Y0;
    
    bool c = from_P2_to_PB.c_bit[cur_depth][mul];
    
    bool other_blinded = received_view[cur_depth].next_bit_L_recv[mul];
 
    bool gamma = xor_if(c, other_blinded, d0);
   
    result = xor_if(gamma, val_mX, (val_bool ^ Y1));    
  }



// i == n
void Verifier::multiply_mpc(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], block_t val_mX, bool val_bool, size_t cur_depth, size_t mul, bool party, block_t& result)
{   

   block_t D0 = PBrand.next_block();   
   block_t X0 = val_mX ^ D0;   
   
   bool d0 = PBrand.next_bool();
   bool Y0 = val_bool ^ d0;   
  
   block_t X1;
   bool Y1;
   block_t gamma;
  
   X1 = received_view[cur_depth].blinds_recv[mul]; 
   Y1 = received_view[cur_depth].bit_blinds_recv[mul]; 

   block_t c = from_P2_to_PB.c[cur_depth][mul];
   
   block_t other_blinded = received_view[cur_depth].blinds_recv[mul];


   gamma = xor_if(c, other_blinded, d0);
 

   generated_view[cur_depth].blinds_recv[mul]          = X0; // generating the other transcript 
   generated_view[cur_depth].bit_blinds_recv[mul]      = Y0; // generating the other transcript
 
   result = xor_if(gamma, val_mX, (val_bool ^ Y1));
}

template<typename prgkey_t>
void Verifier::conditional_swap_and_next_seed(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], block_t s0_L, block_t s0_R, block_t s1_L, block_t s1_R, prgkey_t & key, bool party, size_t cur_depth)
{

  block_t bs0_L;

  multiply_mpc(from_P2_to_PB, received_view,  generated_view ,s0_L, Pdirection[cur_depth],  cur_depth, 0, party, bs0_L);

  block_t notbs0_L = s0_L ^ bs0_L;
   
  block_t bs1_L;
 
  multiply_mpc(from_P2_to_PB, received_view, generated_view , s1_L, Pdirection[cur_depth], cur_depth , 1, party, bs1_L); 

  block_t notbs1_L = s1_L ^ bs1_L;
  
  block_t bs0_R;
 
  multiply_mpc(from_P2_to_PB, received_view, generated_view , s0_R, Pdirection[cur_depth],   cur_depth , 2, party, bs0_R);

  block_t bNotR0_share = s0_R ^ bs0_R;

  block_t bs1_R;

  multiply_mpc(from_P2_to_PB, received_view,generated_view , s1_R, Pdirection[cur_depth],  cur_depth , 3, party, bs1_R);
  
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
void Verifier::root_layer(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], const prgkey_t & key,  dpf_key<leaf_t, node_t, prgkey_t> dpfkey, bool party)
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
  
   if(!party) generated_view[0].L_shares_recv  =  s[L] ^ L_shares[0];  
   if(!party) generated_view[0].R_shares_recv  =  s[R] ^ R_shares[0];  

   if( party) generated_view[0].L_shares_recv  =  L_shares[0];  
   if( party) generated_view[0].R_shares_recv  =  R_shares[0];

   s_L_othershare = received_view[0].L_shares_recv;    
   s_R_othershare = received_view[0].R_shares_recv;   

   if(!party) conditional_swap_and_next_seed(from_P2_to_PB, received_view, generated_view , s_L_myshare, s_R_myshare, s_L_othershare, s_R_othershare, key, party, cur_depth); 
   if( party) conditional_swap_and_next_seed(from_P2_to_PB, received_view, generated_view , s_L_othershare, s_R_othershare, s_L_myshare, s_R_myshare, key, party, cur_depth);
   
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

    generated_view[0].bit_L_shares_recv = bitL_shares[1];
    generated_view[0].bit_R_shares_recv = bitR_shares[1];
   }

   if(party)
   {
    bit_L_myshare = bitL_shares[1];
    bit_R_myshare = bitR_shares[1];

    generated_view[0].bit_L_shares_recv = bitL_shares[0];
    generated_view[0].bit_R_shares_recv = bitR_shares[0];
   }
 
   bit_L_othershare =  received_view[0].bit_L_shares_recv;
   bit_R_othershare =  received_view[0].bit_R_shares_recv; 
   
   
   if(!party) get_next_bits(from_P2_to_PB, received_view, generated_view , bit_L_myshare, bit_R_myshare, bit_L_othershare, bit_R_othershare, party, cur_depth);
   if( party) get_next_bits(from_P2_to_PB, received_view, generated_view , bit_L_othershare, bit_R_othershare, bit_L_myshare, bit_R_myshare, party, cur_depth);
}

template<typename leaf_t, typename node_t, typename prgkey_t>
void Verifier::middle_layers(const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[], const prgkey_t & key, const block_t seedL, const block_t seedR, dpf_key<leaf_t, node_t, prgkey_t> dpfkey, size_t cur_depth, bool party)
{
  
  auto & cw = dpfkey.cw[cur_depth];
  uint8_t cw_t[2] = { get_lsb(cw, 0b01), get_lsb(cw, 0b10) };
  auto & final = dpfkey.finalizer;


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
  
  if(!party) expand_mpc_verify(key, seedL, from_P2_to_PB, received_view, generated_view , s0_L, s0_R, t0_L, t0_R, P0blinds0, P0gamma[0], party, cur_depth, false);  
  if(party ) expand_mpc_verify(key, seedL, from_P2_to_PB, received_view, generated_view , s0_L, s0_R, t0_L, t0_R, P0blinds0, P0gamma[1], party, cur_depth, false);
  
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
  
  if(!party) expand_mpc_verify(key, seedR,  from_P2_to_PB,  received_view,generated_view ,  s1_L, s1_R, t1_L, t1_R,  P1blinds0, P1gamma[0], party, cur_depth, true);
  if(party)  expand_mpc_verify(key, seedR,  from_P2_to_PB,  received_view,generated_view ,  s1_L, s1_R, t1_L, t1_R,  P1blinds0, P1gamma[1], party, cur_depth, true);
 
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


   if(!party) conditional_swap_and_next_seed(from_P2_to_PB, received_view,generated_view , L0_shares_2, R0_shares_2, L1_shares_2, R1_shares_2,  key, party, cur_depth);
   if( party) conditional_swap_and_next_seed(from_P2_to_PB, received_view,generated_view , L0_shares_2, R0_shares_2, L1_shares_2, R1_shares_2,  key, party, cur_depth);

   bool bitL_shares, bitR_shares;
  
   bool bit0_L, bit0_R, bit1_L, bit1_R;
  
   bool my_val = cw_t[L];

   if(party) my_val = bit0[cur_depth];
  
   multiply_mpc(from_P2_to_PB, received_view,generated_view , my_val,  party, 0, bit0_L, cur_depth);

   if(!party) bit0_L = bit0_L ^ t0_L ^ cw_t[L]  &  (bit0[cur_depth]);
  
  if(party)
  {
  
    bit0_L = bit0_L ^ t0_L;
  } 
  
  my_val = cw_t[R];


  if(party) my_val = bit0[cur_depth];

  multiply_mpc(from_P2_to_PB, received_view,generated_view , my_val, party, 1, bit0_R, cur_depth);

  if(!party) bit0_R = bit0_R ^ t0_R ^ (cw_t[R] & bit0[cur_depth]);

  if(party)
  {
    bit0_R = bit0_R ^ t0_R;
  } 
  
  my_val = bit1[cur_depth];
 
 
  if(party) my_val = cw_t[L];

  multiply_mpc(from_P2_to_PB, received_view,generated_view , my_val,  party, 2, bit1_L, cur_depth);
 
   if( party)  bit1_L = bit1_L ^  t1_L ^ (cw_t[L] & bit1[cur_depth]);

  if(!party) bit1_L = bit1_L ^ t1_L;
 
  if(party) my_val = cw_t[R];

  multiply_mpc(from_P2_to_PB, received_view,generated_view , my_val,  party, 3, bit1_R, cur_depth);

  if(!party) bit1_R = bit1_R ^ t1_R;
 
  if(party)  bit1_R = bit1_R ^ (t1_R) ^ (cw_t[R] & bit1[cur_depth]);
  
  get_next_bits(from_P2_to_PB, received_view, generated_view , bit0_L, bit0_R, bit1_L, bit1_R,  party, cur_depth);


    if(cur_depth == depth - 1)
    {
 
      for(size_t j = 0; j < dpfkey.nodes_per_leaf; ++j)
      {
 
         seed0[cur_depth+1] = xor_if(seed0[cur_depth+1], final[j], bit0[cur_depth+1]);
 
         seed1[cur_depth+1] = xor_if(seed1[cur_depth+1], final[j], bit1[cur_depth+1]);
        
         generated_view[cur_depth].final =  (seed0[cur_depth+1]   ^ seed1[cur_depth+1]);     
      }

      std::cout << "final[j] " << final[0].bits << std::endl;
      block_t final_out = (seed0[cur_depth+1]   ^ seed1[cur_depth+1]) ^ received_view[cur_depth].final;
      std::cout << " -> [Verifier]: " << final_out.bits << std::endl;
    }  
   
 }


 template<typename prgkey_t>
 inline auto Verifier::prg_mpc_verify(const prgkey_t & key, const block_t seed, const from_P2 from_P2_to_PB, 
                               const from_PB_ received_view[], from_PB_ generated_view[], block_t blind[], block_t gamma[], size_t cur_depth, bool party, bool LR , size_t len)
 {
   block_t c2L[rounds+1];

   if(!LR) for(size_t r = 0; r <= rounds; ++r) c2L[r] = received_view[cur_depth].seed0L_encrypt[r];   
   if(LR)  for(size_t r = 0; r <= rounds; ++r) c2L[r] = received_view[cur_depth].seed1L_encrypt[r];   
  
   const block_t seed_ = clear_lsb(seed, 0b11); 
  // std::cout << "seed_ " << seed_.bits << std::endl;
   auto tmp = _mm256_xor_si256(seed_, _mm256_set_epi64x(0, 0, 0, 0));
   
    block_t outL[rounds+1];
    auto ccc = key.encrypt_MPC_verify(tmp, c2L, blind, gamma, party, outL);

   if(!LR) for(size_t r = 0; r <= rounds; ++r)  generated_view[cur_depth].seed0L_encrypt[r] = outL[r];
   if(LR)  for(size_t r = 0; r <= rounds; ++r)  generated_view[cur_depth].seed1L_encrypt[r] = outL[r];

   outL[rounds] = _mm256_xor_si256(outL[rounds] , tmp);
   
   block_t c2R[rounds+1];
   
   if(!LR) 
   {
    for(size_t r = 0; r <= rounds; ++r) c2R[r] = received_view[cur_depth].seed0R_encrypt[r];  
   }

   if(LR) 
   {
    for(size_t r = 0; r <= rounds; ++r) c2R[r] =  received_view[cur_depth].seed1R_encrypt[r];  
   }

   if(!party) tmp = _mm256_xor_si256(seed_, _mm256_set_epi64x(0, 0, 0, 1));
   if(party) tmp = seed_;
   
   block_t outR[rounds+1]; 
   key.encrypt_MPC_verify(tmp, c2R, blind, gamma, party, outR);
   
   if(!LR)for(size_t r = 0; r <= rounds; ++r) generated_view[cur_depth].seed0R_encrypt[r] = outR[r];
   if(LR) for(size_t r = 0; r <= rounds; ++r) generated_view[cur_depth].seed1R_encrypt[r] = outR[r];

   outR[rounds] = _mm256_xor_si256(outR[rounds] , tmp);

   return std::make_pair(std::move(outL[rounds]), std::move(outR[rounds]));
 }

 

template<typename prgkey_t>
inline void Verifier::expand_mpc_verify(const prgkey_t & key, const block_t seed, const from_P2 from_P2_to_PB, const from_PB_ received_view[], from_PB_ generated_view[],  block_t &s0_L, block_t & s0_R, 
                                        uint8_t & t0_L, uint8_t & t0_R,  block_t blind[], block_t gamma[], bool party, size_t cur_depth, bool LR)
{ 
   auto prg_out =  prg_mpc_verify(key, seed, from_P2_to_PB, received_view, generated_view, blind,  gamma,  cur_depth, party, LR, 2);

   auto outL = prg_out.first;
   auto outR = prg_out.second;

   s0_L = outL;
   s0_R = outR; 
     
   t0_L = get_lsb(s0_L); 
   t0_R = get_lsb(s0_R);

   s0_L = clear_lsb(s0_L, 0b11);   
   s0_R = clear_lsb(s0_R, 0b11);     
}
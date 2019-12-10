class Verifier
{

  public:
  Verifier(AES_KEY& aeskey, __m128i seed, size_t len)
      : Prand(aeskey, seed, len)
  {


  }
  std::bitset<depth> Pdirection;
  std::bitset<depth> P1direction;

  block seed0[depth][2], seed1[depth][2];
  bool bit0[depth][2], bit1[depth][2];
  void get_next_bits(bool L0_shares, bool R0_shares, bool L1_shares, bool R1_shares, const PB_Transcripts PB_message, bool party ,size_t cur_depth);
  
  void Gen_Blinds(block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds]);
  
  template<typename KEY_TYPE, typename __mX>
  void expand_mpc_verify(KEY_TYPE & key, const blocks<__mX> seed, const PB_Transcripts PB_message, block &s0_L, block& s0_R, 
                                        bool& t0_L, bool& t0_R, block blind[rounds], block gamma[rounds], bool party, size_t cur_depth, bool LR);
  
  template<typename KEY_TYPE, typename __mX>  
  void root_layer(PB_Transcripts PB_message, KEY_TYPE& key,  dpf_key<__mX, nitems> dpfkey,  bool party);
  

  void multiply_mpc(block val_mX, bool val_bool,  const PB_Transcripts PB_message, size_t cur_depth, size_t mul, bool party, block& out);
  void multiply_mpc_b(bool val_mX, bool val_bool, bool blind_mX, bool blind_bool, const PB_Transcripts PB_message, bool party, size_t mul , bool& result, size_t cur_depth);
  
  template<typename KEY_TYPE>
  void conditional_swap(block L_share, block R_share, block L1_share, block R1_share, const PB_Transcripts PB_message, KEY_TYPE& key,  bool party, size_t cur_depth);
  
  void middle_layers(LowMC& key, block seed0, block seed1, const PB_Transcripts PB_message, dpf_key<__mX, nitems> dpfkey, size_t cur_depth, bool party);

  ~ Verifier()
  {

  }

  //private:
  MPCrandomness Prand;

};


template<typename KEY_TYPE, typename __mX>
inline void Verifier::expand_mpc_verify(KEY_TYPE & key, const blocks<__mX> seed, const PB_Transcripts PB_message, block &s0_L, block& s0_R, 
                                        bool& t0_L, bool& t0_R,  block blind[rounds], block gamma[rounds], bool party, size_t cur_depth, bool LR)
{

   const blocks<__mX> seedL = clear_lsb(seed); // _mm_clearlsb_si128(seed);
   const blocks<__mX> seedR = set_lsb(seed); //  _mm_setlsb_si128(seed);
   

   std::cout << "cur_depth [VER] = " << cur_depth << std::endl;
  std::vector<block> c2L;
   if(!LR) c2L = PB_message.PB_middle.TL[cur_depth];
   if(LR) c2L = PB_message.PB_middle.TL2[cur_depth];
 
  std::vector<block> outL = key.encrypt_MPC_verify(seedL, c2L, blind, gamma, party);


   for(size_t j = 1; j < rounds; ++j)
   {    
   //if(!party) if(!LR)  std::cout << "outL[" << j << "] = " << ( outL[j] ^ PB_message.PB_middle.TL_other[cur_depth][j]) << std::endl;
   //if(!party) if(LR)  std::cout << "outL[" << j << "] = " << ( outL[j] ^ PB_message.PB_middle.TL_other2[cur_depth][j]) << std::endl;
   }

   std::vector<block> c2R; 
   
   if(!LR) c2R = PB_message.PB_middle.TR[cur_depth];
   if(LR) c2R = PB_message.PB_middle.TR2[cur_depth];
   std::vector<block> outR = key.encrypt_MPC_verify(seedR, c2R, blind, gamma, party);

   for(size_t j = 0; j <= rounds; ++j)
   {    
   //if(!party) if(!LR) std::cout << "outR[" << j << "] = " << /*"(" << outR[j]  << " ^ " << PB_message.PB_middle.TR_other[cur_depth][j] << ")" <<*/ (outR[j] ^ PB_message.PB_middle.TR_other[cur_depth][j])<< std::endl;
   //if(!party) if(LR) std::cout << "outR[" << j << "] = " << /*"(" << outR[j]  << " ^ " << PB_message.PB_middle.TR_other[cur_depth][j] << ")" <<*/ (outR[j] ^ PB_message.PB_middle.TR_other2[cur_depth][j])<< std::endl;
   }

    s0_L = outL[rounds];
    s0_R = outR[rounds];

    // block s1_L = outL[rounds];
    // block s1_R = outR[rounds];

    s0_L ^= seedL; 
     
    t0_L = get_lsb(s0_L); 
    // t1[L] = get_lsb(s1[L]);
     
     s0_L = clear_lsb(s0_L); 
    
     s0_R ^= seedR; 
    
    t0_R = get_lsb(s0_R); 
    //t1[R] = get_lsb(s1[R]);
    
     s0_R = clear_lsb(s0_R); 
    
    // c0[L] = s0[L];
    // c1[L] = s1[L];

    // c0[R] = s0[R];
    // c1[R] = s1[R];

}

  void Verifier::get_next_bits(bool L0_shares, bool R0_shares, bool L1_shares, bool R1_shares, const PB_Transcripts PB_message, bool party, size_t cur_depth)
  {
    

    bool blind_L,  blind_R;

    if(cur_depth == 0)
    {
    blind_L = PB_message.PB_root.next_bit_L_recv[0];
    blind_R = PB_message.PB_root.next_bit_R_recv[0];
    }
    else
    {
    blind_L = PB_message.PB_middle.next_bit_L_recv[cur_depth][0];
    blind_R = PB_message.PB_middle.next_bit_R_recv[cur_depth][0];
    }
    bool bNotL0_share;
    if(!party) multiply_mpc_b(L0_shares, Pdirection[cur_depth], blind_L, blind_R, PB_message, party, 0, bNotL0_share, cur_depth);
    if(party) multiply_mpc_b(L0_shares, !Pdirection[cur_depth], blind_L, blind_R, PB_message, party, 0, bNotL0_share, cur_depth);
   
   
    if(cur_depth == 0)
    {
    blind_L = PB_message.PB_root.next_bit_L_recv[1];
    blind_R = PB_message.PB_root.next_bit_R_recv[1];
    }
    else
    {
    blind_L = PB_message.PB_middle.next_bit_L_recv[cur_depth][1];
    blind_R = PB_message.PB_middle.next_bit_R_recv[cur_depth][1];
    }
    bool bNotL1_share;
    if(!party) multiply_mpc_b(L1_shares, Pdirection[cur_depth] , blind_L, blind_R, PB_message, party, 1, bNotL1_share, cur_depth);
    if(party) multiply_mpc_b(L1_shares, !Pdirection[cur_depth] , blind_L, blind_R, PB_message, party, 1, bNotL1_share, cur_depth);
      
    if(cur_depth == 0)
    {
    blind_L = PB_message.PB_root.next_bit_L_recv[2];
    blind_R = PB_message.PB_root.next_bit_R_recv[2];
        }
    else
    {
    blind_L = PB_message.PB_middle.next_bit_L_recv[cur_depth][2];
    blind_R = PB_message.PB_middle.next_bit_R_recv[cur_depth][2];
    }
    bool bR0_share;
    multiply_mpc_b(R0_shares, Pdirection[cur_depth], blind_L, blind_R, PB_message, party, 2, bR0_share, cur_depth);

    if(cur_depth == 0)
    {
    blind_L = PB_message.PB_root.next_bit_L_recv[3];
    blind_R = PB_message.PB_root.next_bit_R_recv[3];
    }
    else
    {
    blind_L = PB_message.PB_middle.next_bit_L_recv[cur_depth][3];
    blind_R = PB_message.PB_middle.next_bit_R_recv[cur_depth][3];
    }
    bool bR1_share;
    multiply_mpc_b(R1_shares, Pdirection[cur_depth], blind_L, blind_R, PB_message, party, 3, bR1_share, cur_depth);


    bit0[cur_depth+1][0] = bNotL0_share  ^ bR0_share;    
    bit1[cur_depth+1][0] = bNotL1_share  ^ bR1_share;

    std::cout << "bit0: [" << cur_depth << "] = " << bit0[cur_depth + 1][0] << std::endl; 
    std::cout << "bit1: [" << cur_depth << "] = " << bit1[cur_depth + 1][0] << std::endl;
    
  }

  void Verifier::Gen_Blinds(block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds])
  {
    block rand[rounds];
    for (unsigned r = 0; r < rounds; ++r)
     {
        Prand.next(blinds0[r]);
        // Prand.next(blinds1[r]); 
        // Prand.next(rand[r]);
    
        // const block tmp1 = ((blinds0[r] >> 1) & blinds1[r]) ^ ((blinds1[r] >> 1) & blinds0[r]); 
        // const block tmp2 = ((blinds0[r] >> 2) & blinds1[r]) ^ ((blinds1[r] >> 2) & blinds0[r]);
    
        // const block bc = (tmp1 << 2) & maska;
        // const block ac = (tmp2 << 1) & maskb;
        // const block ab = (tmp1 >> 1) & maskc;
    
        // gamma[0][r] = (bc | ac | ab) ^ rand[r];
        // gamma[1][r] = rand[r];
     }
  }
void Verifier::middle_layers(LowMC& key, block seed0, block seed1, const PB_Transcripts PB_message, dpf_key<__mX, nitems> dpfkey, size_t cur_depth, bool party)
{
  
  std::cout << "middle layer: " << cur_depth << std::endl;
  block CW = dpfkey.cw[cur_depth];

  block blind[rounds];
  block gamma[rounds];
  
  block P0gamma[2][rounds];
  block P0blinds0[rounds];
  block P0blinds1[rounds];

  Gen_Blinds(P0blinds0 ,P0blinds1, P0gamma);
  
  for(size_t j = 0; j < rounds; ++j)
  {
  P0gamma[0][j] = PB_message.PB_middle.gamma0[cur_depth][j];
  P0gamma[1][j] = PB_message.PB_middle.gamma0[cur_depth][j];
  } 
  

  block s0_L,  s0_R; 
  bool t0_L, t0_R;
  
  if(!party)  expand_mpc_verify(key, seed0, PB_message, s0_L, s0_R, t0_L, t0_R, P0blinds0, P0gamma[0], party, cur_depth, false);  
  if(party)   expand_mpc_verify(key, seed0, PB_message, s0_L, s0_R, t0_L, t0_R, P0blinds1, P0gamma[1], party, cur_depth, false);
  

  std::cout << "bit0 = " << bit0[cur_depth][0] << std::endl;  
  std::cout << "seed0 = " << seed0 << std::endl;  

  block L0_shares_2, R0_shares_2;
  if(!party)
  {
   L0_shares_2 = xor_if(s0_L, CW, bit0[cur_depth][0]);
   R0_shares_2 = xor_if(s0_R, CW, bit0[cur_depth][0]);
  }
  if(party)
  {
   L0_shares_2 = xor_if(s0_L, CW, !bit0[cur_depth][0]);
   R0_shares_2 = xor_if(s0_R, CW, !bit0[cur_depth][0]);
  }

  block P1gamma[2][rounds];
  block P1blinds0[rounds];
  block P1blinds1[rounds];

  Gen_Blinds(P1blinds0 ,P1blinds1, P1gamma); 

  for(size_t j = 0; j < rounds; ++j)
  {
  P1gamma[0][j] = PB_message.PB_middle.gamma1[cur_depth][j];
  P1gamma[1][j] = PB_message.PB_middle.gamma1[cur_depth][j];
  } 

  if(!party) expand_mpc_verify(key, seed1, PB_message, s0_L, s0_R, t0_L, t0_R,  P1blinds0, P1gamma[0], party, cur_depth, true);
  if(party)  expand_mpc_verify(key, seed1, PB_message, s0_L, s0_R, t0_L, t0_R,  P1blinds1, P1gamma[1], party, cur_depth, true);

  block L0_shares, R0_shares;
  block L1_shares, R1_shares;
  
  block L1_shares_2, R1_shares_2;
  if(!party)
  {
   L1_shares_2 = xor_if(s0_L, CW, bit1[cur_depth][0]);
   R1_shares_2 = xor_if(s0_R, CW, bit1[cur_depth][0]);
  }
  if(party)
  {
   L1_shares_2 = xor_if(s0_L, CW, !bit1[cur_depth][0]);
   R1_shares_2 = xor_if(s0_R, CW, !bit1[cur_depth][0]);
  }
  L0_shares = L0_shares_2;

  if(!party) std::cout << "ZERO = " << (L0_shares_2 ^ PB_message.PB_middle.L_shares_sent[cur_depth]) << std::endl;

  R0_shares = R0_shares_2;

  if(!party) std::cout << "ZERO = " << (R0_shares_2 ^ PB_message.PB_middle.R_shares_sent[cur_depth]) << std::endl;


  L1_shares = PB_message.PB_middle.L_shares_recv[cur_depth];
  R1_shares = PB_message.PB_middle.R_shares_recv[cur_depth];

  if(!party) conditional_swap(L0_shares, R0_shares, L1_shares, R1_shares, PB_message, key, party, cur_depth);
  if( party) conditional_swap(L1_shares, R1_shares, L0_shares, R0_shares, PB_message, key, party, cur_depth);

  bool bitL_shares, bitR_shares;
  
  bool L0_share, R0_share, L1_share, R1_share;
  L0_share = PB_message.PB_middle.bit_L_shares_sent[cur_depth];
  R0_share = PB_message.PB_middle.bit_R_shares_sent[cur_depth];
  L1_share = PB_message.PB_middle.bit_L_shares_recv[cur_depth];
  R1_share = PB_message.PB_middle.bit_R_shares_recv[cur_depth];
  
  if(!party) get_next_bits(L0_share, R0_share, L1_share, R1_share, PB_message, party, cur_depth);
  if( party) get_next_bits(L1_share, R1_share, L0_share, R0_share, PB_message, party, cur_depth);
}
  
 
  void Verifier::multiply_mpc_b(bool val_mX, bool val_bool, bool blind_mX, bool blind_bool, const PB_Transcripts PB_message, bool party, size_t mul, bool& result, size_t cur_depth)
  {
    int dd;
    bool D0 = Prand.next(dd);
    bool X0 = val_mX ^ D0;// PB_message.PB_root.blinds_sent[mul];
    


    int dd0;     
    bool d0 = Prand.next(dd);
    bool Y0 = val_bool ^ d0;  //PB_message.PB_root.bit_blinds_sent[mul];
    
    bool X1 = PB_message.PB_root.next_bit_L_recv[mul];
    bool Y1 = PB_message.PB_root.next_bit_R_recv[mul];
    bool gamma = PB_message.PB_root.gamma_bit[mul];
    if(cur_depth >= 1 )
    {
    X1    = PB_message.PB_middle.next_bit_L_recv[cur_depth][mul];
    Y1    = PB_message.PB_middle.next_bit_R_recv[cur_depth][mul];
    gamma = PB_message.PB_middle.gamma_bit[cur_depth][mul];
    }

    result = xor_if(gamma, val_mX, (val_bool ^ Y1));    
  }


void Verifier::multiply_mpc(block val_mX, bool val_bool, const PB_Transcripts PB_message, size_t cur_depth, size_t mul, bool party, block& result)
{
 
  block D0; Prand.next(D0);
  //std::cout << "D0 = " << D0 << std::endl;
  block X0 = val_mX ^ D0;// PB_message.PB_root.blinds_sent[mul];
  
  
  int dd0;     
  bool d0 = Prand.next(dd0);
  bool Y0 = val_bool ^ d0;  //PB_message.PB_root.bit_blinds_sent[mul];
  
  block X1;
  bool Y1;
  block gamma;
  
  if(cur_depth == 0)
  {
    X1 = PB_message.PB_root.blinds_recv[mul];
    Y1  = PB_message.PB_root.bit_blinds_recv[mul];
    gamma = PB_message.PB_root.gamma[mul];  
  }
  else
  {
    X1 = PB_message.PB_middle.blinds_recv[cur_depth][mul];
    Y1  = PB_message.PB_middle.bit_blinds_recv[cur_depth][mul];
   
    gamma = PB_message.PB_middle.gamma_middle[cur_depth][mul];
  }
  
  result = xor_if(gamma, val_mX, (val_bool ^ Y1));
}

template<typename KEY_TYPE>
void Verifier::conditional_swap(block L0_share, block R0_share, block L1_share, block R1_share, const PB_Transcripts PB_message, KEY_TYPE& key, bool party, size_t cur_depth)
{

  block bL0_share;

  multiply_mpc(L0_share, Pdirection[cur_depth], PB_message,cur_depth, 0, party, bL0_share);

  // if(cur_depth == 1) 
  // {
    
  //   std::cout << "L0_share = " << L0_share << std::endl;
  //   std::cout << "L1_share = " << L1_share << std::endl;
  //   std::cout << "R0_share = " << R0_share << std::endl;
  //   std::cout << "R1_share = " << R1_share << std::endl;
  //   std::cout << "bL0_share = " << bL0_share << std::endl;
  // }
  block bNotL0share = L0_share ^ bL0_share;
   
  block bL1_share;
 
  multiply_mpc(L1_share, Pdirection[cur_depth], PB_message,cur_depth , 1, party, bL1_share); 

  
  block bNotL1_share = L1_share ^ bL1_share;
  
  block bR0_share;
 
  multiply_mpc(R0_share, Pdirection[cur_depth],  PB_message,cur_depth , 2, party, bR0_share);
  block bNotR0_share = R0_share ^ bR0_share;


  block bR1_share;
  multiply_mpc(R1_share, Pdirection[cur_depth], PB_message,cur_depth , 3, party, bR1_share);
  
  block bNotR1_share = bR1_share ^ R1_share;

  if(!party)
  {
  seed0[cur_depth+1][0]  = bNotL0share  ^ bR0_share;
  seed1[cur_depth+1][0]  = bNotL1_share ^ bR1_share;
  }

  if(party)
  {
  seed0[cur_depth+1][0] = bNotL0share  ^ bR0_share;
  seed1[cur_depth+1][0] = bNotL1_share ^ bR1_share;
  }


  std::cout << "seed0[cur_depth+1][0] = " << (seed0[cur_depth+1][0] ^ PB_message.PB_middle.seedL[cur_depth+1]) << std::endl;
  std::cout << "seed1[cur_depth+1][0] = " << (seed1[cur_depth+1][0] ^ PB_message.PB_middle.seedR[cur_depth+1]) << std::endl << std::endl;
}

template<typename KEY_TYPE, typename __mX>
void Verifier::root_layer(const PB_Transcripts PB_message, KEY_TYPE& key,  dpf_key<__mX, nitems> dpfkey, bool party)
{

   block seed = dpfkey.root;
   const block seedL = clear_lsb(seed);
   const block seedR = set_lsb(seed);
   
   size_t cur_depth = 0;

   block CW = dpfkey.cw[cur_depth];
   bool b = get_lsb(dpfkey.root);
   
   block s[2];
   bool  t[2];
   
   block child[2];

   expand_nonmpc(key, seed, child, t);
   
   s[L] = xor_if(child[L], CW, !b);
   s[R] = xor_if(child[R], CW, !b);
   
   block L_shares[2]; 
   
   Prand.next(L_shares[0]);
   
   L_shares[1] = s[L] ^ L_shares[0];
  
   
   block R_shares[2]; 
   Prand.next(R_shares[0]);
  
   R_shares[1] = s[R] ^ R_shares[0];
  

   block L1_share = PB_message.PB_root.L_shares_recv;
   block R1_share = PB_message.PB_root.R_shares_recv;

   if(!party) conditional_swap(L_shares[0], R_shares[0], L1_share, R1_share,  PB_message, key, party, cur_depth);
   
   if(party) conditional_swap( L1_share, R1_share, L_shares[1], R_shares[1], PB_message, key, party, cur_depth);
   
   bool bitL_shares[2];  
   bool bitR_shares[2]; 

   bool bit_L0_share, bit_R0_share, bit_L1_share, bit_R1_share;
   
   int dd; 

   bitL_shares[0] = Prand.next(dd);
   bitR_shares[0] = Prand.next(dd);

   bitL_shares[1] = bitL_shares[0] ^ t[L] ^ dpfkey.t[cur_depth][L] & (b);
   bitR_shares[1] = bitR_shares[0] ^ t[R] ^ dpfkey.t[cur_depth][R] & (b);
   
   bit_L1_share = PB_message.PB_root.bit_L_shares_recv;
   bit_R1_share = PB_message.PB_root.bit_R_shares_recv;

   if(!party) get_next_bits(bitL_shares[0], bitR_shares[0], bit_L1_share, bit_R1_share, PB_message, party, cur_depth);
   if( party) get_next_bits(bit_L1_share, bit_R1_share, bitL_shares[1], bitR_shares[1], PB_message, party, cur_depth);
}
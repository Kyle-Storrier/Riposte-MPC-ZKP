class Verifier
{

  public:
  Verifier(AES_KEY& aeskey, __m128i seed0, __m128i seed1, __m128i seed2,  size_t len)
      : P0rand(aeskey, seed0, len), P1rand(aeskey, seed1, len), P2rand(aeskey, seed2, len)
  {


  }
  std::bitset<depth> P0direction;
  std::bitset<depth> P1direction;

  block seed0[depth][2], seed1[depth][2];
  template<typename KEY_TYPE, typename __mX>  
  void root_layer(PB_Transcripts PB_message, KEY_TYPE& key,  dpf_key<__mX, nitems> dpfkey,  bool party);
  

  void multiply_mpc(block val_mX, bool val_bool,  block blind_mX, bool blind_bool, PB_Transcripts PB_message, size_t mul, bool party, block& out);

  template<typename KEY_TYPE>
  void conditional_swap(block L_share, block R_share, bool val_bool, PB_Transcripts PB_message, KEY_TYPE& key,  bool party, size_t cur_depth);
  
  void middle_layers(LowMC& key,  dpf_key<__mX, nitems> dpfkey, size_t cur_depth);

  ~ Verifier()
  {

  }

  private:
  MPCrandomness P0rand, P1rand, P2rand;

};


void Verifier::middle_layers(LowMC& key,  dpf_key<__mX, nitems> dpfkey, size_t cur_depth)
{
  block CW = dpfkey.cw[cur_depth];
}

void Verifier::multiply_mpc(block val_mX, bool val_bool, block blind_mX, bool blind_bool, PB_Transcripts PB_message, size_t mul, bool party, block& result)
{
 
  block X0 = PB_message.PB_root.blinds_sent[mul];
  bool Y0  = PB_message.PB_root.bit_blinds_sent[mul];
  
  block X1 = PB_message.PB_root.blinds_recv[mul];
  bool Y1  = PB_message.PB_root.bit_blinds_recv[mul];
  
  block gamma = PB_message.PB_root.gamma[mul];

 // std::cout << "gamma = " << gamma << std::endl;
  result = xor_if(gamma, val_mX, (val_bool ^ Y1));

//  std::cout << "val_mX = " << val_mX   << std::endl;
//  std::cout << "val_bool = " << val_bool << std::endl;
  std::cout << "result[" << mul << "] = " << result << std::endl << std::endl;

}

template<typename KEY_TYPE>
void Verifier::conditional_swap(block L_share, block R_share, bool val_bool, PB_Transcripts PB_message, KEY_TYPE& key, bool party, size_t cur_depth)
{

  block bL0_share;
  multiply_mpc(L_share, val_bool, PB_message.PB_root.blinds_recv[0], PB_message.PB_root.bit_blinds_recv[0], PB_message, 0, party, bL0_share);
  
  std::cout << "L_share = " << L_share << std::endl;

  block bNotL0share = L_share ^ bL0_share;
  std::cout << "bNotL0share = " << bNotL0share << std::endl;

  block bL1_share;
  multiply_mpc(PB_message.PB_root.L_shares_recv, val_bool, PB_message.PB_root.blinds_recv[1], PB_message.PB_root.bit_blinds_recv[1], PB_message, 1, party, bL1_share);
  block bNotL1_share = PB_message.PB_root.L_shares_recv ^ bL1_share;

  block bR0_share;
  multiply_mpc(R_share, val_bool, PB_message.PB_root.blinds_recv[2], PB_message.PB_root.bit_blinds_recv[2], PB_message, 2, party, bR0_share);
  block bNotR0_share = R_share ^ bR0_share;
  
  block bR1_share;
  multiply_mpc(PB_message.PB_root.R_shares_recv, val_bool, PB_message.PB_root.blinds_recv[3], PB_message.PB_root.bit_blinds_recv[3], PB_message, 3, party, bR1_share);
  block bNotR1_share = bR1_share ^ PB_message.PB_root.R_shares_recv;

  seed0[cur_depth+1][0]  = bNotL0share  ^ bR0_share;
  seed1[cur_depth+1][0]  = bNotL1_share ^ bR1_share;

  std::cout << "seed0[cur_depth+1][0] = " << seed0[cur_depth+1][0] << std::endl;
  std::cout << "seed1[cur_depth+1][0] = " << seed1[cur_depth+1][0] << std::endl;
}

template<typename KEY_TYPE, typename __mX>
void Verifier::root_layer(PB_Transcripts PB_message, KEY_TYPE& key,  dpf_key<__mX, nitems> dpfkey, bool party)
{

   block seed = dpfkey.root;
   const block seedL = clear_lsb(seed);
   const block seedR =   set_lsb(seed);
   

   size_t cur_depth = 0;

   block CW = dpfkey.cw[cur_depth];
   bool b = get_lsb(dpfkey.root);
   block s[2];
   bool  t[2];
   block child[2];

   expand_nonmpc(key, seed, child, t);
   s[L] = xor_if(child[L], CW, !b);
   s[R] = xor_if(child[R], CW, !b);
   
   std::cout << "s[" << L << "] = " << s[L] << std::endl;
   std::cout << "s[" << R << "] = " << s[R] << std::endl;

   block L_shares[2]; 
  // P0rand.next(L_shares[0]);
    
   L_shares[1] = PB_message.PB_root.L_shares_sent;
   L_shares[0] = s[L] ^ L_shares[1];

   std::cout << "L_shares = " << L_shares[0] << " " << L_shares[1] << std::endl;

   block R_shares[2]; 
   R_shares[1] = PB_message.PB_root.R_shares_sent;
 //  P0rand.next(R_shares[0]);
   R_shares[0] = s[R] ^ R_shares[1];
   std::cout << "R_shares = " << R_shares[0] << " " << R_shares[1] << std::endl << std::endl << std::endl;

   conditional_swap(L_shares[0], R_shares[0], P0direction[0],  PB_message, key, true, cur_depth);
   
}
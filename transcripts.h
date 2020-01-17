

struct from_P2
{ 

  std::vector<std::vector<block_t>> gamma0; //[depth][rounds];
  std::vector<std::vector<block_t>> gamma1;//[depth][rounds];
  std::vector<std::vector<block_t>> c;//[depth][4];  
  std::vector<std::vector<bool>> c_bit;//[depth][4];
  std::vector<std::vector<bool>> c_bit2;//[depth][4];

  from_P2(size_t depth, size_t rounds)
  {
    gamma0.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) gamma0[j].resize(rounds + 1);
    gamma1.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) gamma1[j].resize(rounds + 1);
    c.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) c[j].resize(4);
    c_bit.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) c_bit[j].resize(4);
    c_bit2.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) c_bit2[j].resize(4);
  }

};

struct from_P2_
{ 
  block_t gamma0 [rounds]; // Used in encrypt
  block_t gamma1[rounds]; // Used in encrypt
  block_t c[4];  
  bool c_bit[4];
  bool c_bit2[4];
};

 struct from_PB_
 { 
  block_t seed0L_encrypt[rounds+1]; // Used in encrypt
  block_t seed0R_encrypt[rounds+1]; // Used in encrypt
  block_t seed1L_encrypt[rounds+1]; // Used in encrypt
  block_t seed1R_encrypt[rounds+1]; // Used in encrypt
  
  block_t blinds_recv[4];   // To do conditional swap
  bool bit_blinds_recv[4];  // To do conditional swap

  bool next_bit_L_recv[4];  
  bool next_bit_R_recv[4];  

  bool next_bit_L_recv2[4]; // Gen Next Bits
  bool next_bit_R_recv2[4]; // Gen Next Bits

  block_t L_shares_recv; // root layer
  block_t R_shares_recv; // root layer  
  block_t final; // final layer
  bool bit_L_shares_recv , bit_R_shares_recv; // root layer  
 };

//   struct  vals_
//   {
//     std::vector<std::vector<block_t>> blinds_recv;
//     std::vector<std::vector<block_t>> seed0L_encrypt; 
//     std::vector<std::vector<block_t>> seed0R_encrypt; 
//     std::vector<std::vector<block_t>> seed1L_encrypt; 
//     std::vector<std::vector<block_t>> seed1R_encrypt; 
   
//     std::vector<std::vector<bool>> bit_blinds_recv; 
//     std::vector<std::vector<bool>> next_bit_L_recv;   
//     std::vector<std::vector<bool>> next_bit_R_recv;   
//     std::vector<std::vector<bool>> next_bit_L_recv2;
//     std::vector<std::vector<bool>> next_bit_R_recv2;
    
//     block_t L_shares_recv; 
//     block_t R_shares_recv;   
//     block_t final;
//     bool bit_L_shares_recv , bit_R_shares_recv;  
    
//   };

// struct from_PB
// {


 
//    size_t depth_;
  
//    vals_ vals;
 

 
//    from_PB(size_t depth, size_t rounds)
//    {
 
//     depth_ = depth;

//     vals.blinds_recv.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.blinds_recv[j].resize(4);
//     vals.bit_blinds_recv.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.bit_blinds_recv[j].resize(4);
    
//     vals.next_bit_R_recv.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.next_bit_R_recv[j].resize(4);
//     vals.next_bit_L_recv.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.next_bit_L_recv[j].resize(4);
    
//     vals.next_bit_R_recv2.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.next_bit_R_recv2[j].resize(4);
//     vals.next_bit_L_recv2.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.next_bit_L_recv2[j].resize(4);
    

//     vals.seed0L_encrypt.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.seed0L_encrypt[j].resize(rounds + 1);
//     vals.seed0R_encrypt.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.seed0R_encrypt[j].resize(rounds + 1);
//     vals.seed1L_encrypt.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.seed1L_encrypt[j].resize(rounds + 1);
//     vals.seed1R_encrypt.resize(depth); for(size_t j = 0; j < depth; ++j ) vals.seed1R_encrypt[j].resize(rounds + 1);
//   }

//   std::string to_string()
//   {
//     std::string string_val = std::to_string(vals.bit_L_shares_recv);

//     string_val = string_val +  vals.L_shares_recv.bits.to_string()  +  vals.R_shares_recv.bits.to_string() + vals.final.bits.to_string();

//     for(size_t i = 0; i < depth_; ++i)
//     {
//       for (size_t j = 0; j < 4; ++j)
//       {
//         string_val = string_val + vals.blinds_recv[i][j].bits.to_string();

//         string_val = string_val + std::to_string(vals.bit_blinds_recv[i][j]) + std::to_string(vals.next_bit_L_recv[i][j]) 
//                                 + std::to_string(vals.next_bit_R_recv[i][j]) + std::to_string(vals.next_bit_L_recv2[i][j]) + std::to_string(vals.next_bit_R_recv2[i][j]);
//       }

//        for(size_t j = 0; j < rounds; ++j) string_val = string_val + vals.seed0L_encrypt[i][j].bits.to_string();
//        for(size_t j = 0; j < rounds; ++j) string_val = string_val + vals.seed0R_encrypt[i][j].bits.to_string();
//        for(size_t j = 0; j < rounds; ++j) string_val = string_val + vals.seed1L_encrypt[i][j].bits.to_string();
//        for(size_t j = 0; j < rounds; ++j) string_val = string_val + vals.seed1R_encrypt[i][j].bits.to_string();
//       }    

//     return string_val;
//   }
// };


struct proof
 {
   std::string final_hash;
 }; 


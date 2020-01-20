

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

  struct from_PB_compressed_
 { 
  // block_t seed0L_encrypt[rounds+1]; // Used in encrypt
  // block_t seed0R_encrypt[rounds+1]; // Used in encrypt
  // block_t seed1L_encrypt[rounds+1]; // Used in encrypt
  // block_t seed1R_encrypt[rounds+1]; // Used in encrypt
  std::vector<bool> seed0L_encrypt;
  std::vector<bool> seed0R_encrypt;
  std::vector<bool> seed1L_encrypt;
  std::vector<bool> seed1R_encrypt;

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

void addToBitVector(std::vector<bool> &v, block_t data, int n) {
  int lengthToStore = n;
  // std::cout << "Original " << data.bits.to_string() << std::endl;
// std::cout << "Compressed ";
  ssize_t currentIndex = 0;
  block_t current = data;
  for(ssize_t j = 0; j < lengthToStore; j++) {
    v.push_back(get_lsb(current) == 1);
    // std::cout << (get_lsb(current) == 1);
    current >>= 1;
    currentIndex++;
  }
  // std::cout << std::endl;
}

block_t extractFromBitVector(std::vector<bool> &v, int n) {
  block_t data;
  bool bit;
  int lengthToRead = n;
  // std::cout << "Compressed ";
  for(int i = 0; i < lengthToRead; i++) {
    bit = v.back();
    // std::cout << bit;
    v.pop_back();
    data <<= 1;
    if(bit) {
      data.bits.set(0);
    }
  }
  // std::cout << "\nDecompressed " << data.bits.to_string() << std::endl;
  return data;
}

void compressTranscript(from_PB_ originalTranscript,
    from_PB_compressed_ & compressedTranscript) {
  compressedTranscript.L_shares_recv = originalTranscript.L_shares_recv;
  compressedTranscript.R_shares_recv = originalTranscript.R_shares_recv;
  compressedTranscript.final = originalTranscript.final;
  compressedTranscript.bit_L_shares_recv = originalTranscript.bit_L_shares_recv;
  compressedTranscript.bit_R_shares_recv = originalTranscript.bit_R_shares_recv;
  
  for(int i = 0; i < 4; i++) {
    compressedTranscript.blinds_recv[i] = originalTranscript.blinds_recv[i];
    compressedTranscript.bit_blinds_recv[i] = originalTranscript.bit_blinds_recv[i];

    compressedTranscript.next_bit_L_recv[i] = originalTranscript.next_bit_L_recv[i];
    compressedTranscript.next_bit_R_recv[i] = originalTranscript.next_bit_R_recv[i];

    compressedTranscript.next_bit_L_recv2[i] = originalTranscript.next_bit_L_recv2[i];
    compressedTranscript.next_bit_R_recv2[i] = originalTranscript.next_bit_R_recv2[i];
  }

  for(ssize_t i = 0; i < rounds + 1; i++) {
    addToBitVector(
        compressedTranscript.seed0L_encrypt, originalTranscript.seed0L_encrypt[i], 256 - 66);
    addToBitVector(
        compressedTranscript.seed0R_encrypt, originalTranscript.seed0R_encrypt[i], 256 - 66);
    addToBitVector(
        compressedTranscript.seed1L_encrypt, originalTranscript.seed1L_encrypt[i], 256 - 66);
    addToBitVector(
        compressedTranscript.seed1R_encrypt, originalTranscript.seed1R_encrypt[i], 256 - 66);
    // compressedTranscript.seed0L_encrypt[i] = originalTranscript.seed0L_encrypt[i];
    // compressedTranscript.seed0R_encrypt[i] = originalTranscript.seed0R_encrypt[i];
    // compressedTranscript.seed1L_encrypt[i] = originalTranscript.seed1L_encrypt[i];
    // compressedTranscript.seed1R_encrypt[i] = originalTranscript.seed1R_encrypt[i];
  }
}

void decompressTranscript(from_PB_compressed_ compressedTranscript,
    from_PB_ & decompressedTranscript) {
  decompressedTranscript.L_shares_recv = compressedTranscript.L_shares_recv;
  decompressedTranscript.R_shares_recv = compressedTranscript.R_shares_recv;
  decompressedTranscript.final = compressedTranscript.final;
  decompressedTranscript.bit_L_shares_recv = compressedTranscript.bit_L_shares_recv;
  decompressedTranscript.bit_R_shares_recv = compressedTranscript.bit_R_shares_recv;
  
  for(int i = 0; i < 4; i++) {
    decompressedTranscript.blinds_recv[i] = compressedTranscript.blinds_recv[i];
    decompressedTranscript.bit_blinds_recv[i] = compressedTranscript.bit_blinds_recv[i];

    decompressedTranscript.next_bit_L_recv[i] = compressedTranscript.next_bit_L_recv[i];
    decompressedTranscript.next_bit_R_recv[i] = compressedTranscript.next_bit_R_recv[i];

    decompressedTranscript.next_bit_L_recv2[i] = compressedTranscript.next_bit_L_recv2[i];
    decompressedTranscript.next_bit_R_recv2[i] = compressedTranscript.next_bit_R_recv2[i];
  }

  // for(ssize_t i = 0; i < rounds + 1; i++) {
  //   decompressedTranscript.seed0L_encrypt[i] = compressedTranscript.seed0L_encrypt[i];
  //   decompressedTranscript.seed0R_encrypt[i] = compressedTranscript.seed0R_encrypt[i];
  //   decompressedTranscript.seed1L_encrypt[i] = compressedTranscript.seed1L_encrypt[i];
  //   decompressedTranscript.seed1R_encrypt[i] = compressedTranscript.seed1R_encrypt[i];
  // }

  for(ssize_t i = rounds; i >= 0; i--) {
    decompressedTranscript.seed0L_encrypt[i]
        = extractFromBitVector(compressedTranscript.seed0L_encrypt, 256 - 66);
    decompressedTranscript.seed0R_encrypt[i]
        = extractFromBitVector(compressedTranscript.seed0R_encrypt, 256 - 66);
    decompressedTranscript.seed1L_encrypt[i]
        = extractFromBitVector(compressedTranscript.seed1L_encrypt, 256 - 66);
    decompressedTranscript.seed1R_encrypt[i]
        = extractFromBitVector(compressedTranscript.seed1R_encrypt, 256 - 66);
  }
}

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


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


struct from_PB
{
  size_t depth_;

  block_t L_shares_recv; 
 
  block_t R_shares_recv; 
  
  std::vector<std::vector<block_t>> blinds_recv;
   
  std::vector<std::vector<bool>> bit_blinds_recv; 
  
  bool bit_L_shares_recv , bit_R_shares_recv;   
 
  std::vector<std::vector<bool>> next_bit_L_recv; 
  
  std::vector<std::vector<bool>> next_bit_R_recv; 
  
  std::vector<std::vector<bool>> next_bit_L_recv2; 

  std::vector<std::vector<bool>> next_bit_R_recv2; 
 
  std::vector<std::vector<block_t>> seed0L_encrypt; 
  std::vector<std::vector<block_t>> seed0R_encrypt; 
  std::vector<std::vector<block_t>> seed1L_encrypt; 
  std::vector<std::vector<block_t>> seed1R_encrypt; 
 
   from_PB(size_t depth, size_t rounds)
   {
    depth_ = depth;

    std::cout << "reserved " <<  depth << std::endl;
    blinds_recv.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) blinds_recv[j].resize(4);
    bit_blinds_recv.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) bit_blinds_recv[j].resize(4);
    
    next_bit_R_recv.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) next_bit_R_recv[j].resize(4);
    next_bit_L_recv.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) next_bit_L_recv[j].resize(4);
    
    next_bit_R_recv2.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) next_bit_R_recv2[j].resize(4);
    next_bit_L_recv2.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) next_bit_L_recv2[j].resize(4);
    

    seed0L_encrypt.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) seed0L_encrypt[j].resize(rounds + 1);
    seed0R_encrypt.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) seed0R_encrypt[j].resize(rounds + 1);
    seed1L_encrypt.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) seed1L_encrypt[j].resize(rounds + 1);
    seed1R_encrypt.resize(depth+1); for(size_t j = 0; j < depth + 1; ++j ) seed1R_encrypt[j].resize(rounds + 1);
   }

   std::string to_string()
  {
    std::string string_val = std::to_string(bit_L_shares_recv);

    string_val = string_val +  L_shares_recv.bits.to_string()  +  R_shares_recv.bits.to_string();



    for(size_t i = 0; i < depth_; ++i)
    {
      for (size_t j = 0; j < 4; ++j)
      {
        string_val = string_val + blinds_recv[i][j].bits.to_string();

        string_val = string_val + std::to_string(bit_blinds_recv[i][j]) + std::to_string(next_bit_L_recv[i][j]) 
                                + std::to_string(next_bit_R_recv[i][j]) + std::to_string(next_bit_L_recv2[i][j]) + std::to_string(next_bit_R_recv2[i][j]);
      }

      for(size_t j = 0; j < rounds; ++j) string_val = string_val +  seed0L_encrypt[i][j].bits.to_string();
     for(size_t j = 0; j < rounds; ++j) string_val = string_val +  seed0R_encrypt[i][j].bits.to_string();
       for(size_t j = 0; j < rounds; ++j) string_val = string_val +  seed1L_encrypt[i][j].bits.to_string();
       for(size_t j = 0; j < rounds; ++j)  string_val = string_val +  seed1R_encrypt[i][j].bits.to_string();
      }    

    return string_val;
  }
};

struct from_PB_compressed
{
  size_t depth_;

  block_t L_shares_recv; 
 
  block_t R_shares_recv; 
  
  // std::vector<std::vector<block_t>> blinds_recv; // To compress
  std::vector<bool> blinds_recv;
   
  std::vector<std::vector<bool>> bit_blinds_recv; 
  
  bool bit_L_shares_recv , bit_R_shares_recv;   
 
  std::vector<std::vector<bool>> next_bit_L_recv; 
  
  std::vector<std::vector<bool>> next_bit_R_recv; 
  
  std::vector<std::vector<bool>> next_bit_L_recv2; 

  std::vector<std::vector<bool>> next_bit_R_recv2; 
 
  std::vector<std::vector<block_t>> seed0L_encrypt; 
  std::vector<std::vector<block_t>> seed0R_encrypt; 
  std::vector<std::vector<block_t>> seed1L_encrypt; 
  std::vector<std::vector<block_t>> seed1R_encrypt; 
};

void addToBitVector(std::vector<bool> &v, block_t data, int n) {
  int lengthToStore = 3*n;

  ssize_t currentIndex = 0;
  block_t current = data;
  for(ssize_t j = 0; j < lengthToStore; j++) {
    v.push_back(get_lsb(current) == 1);
    current >>= 1;
    currentIndex++;
  }
}

block_t extractFromBitVector(std::vector<bool> &v, int n) {
  block_t data;
  bool bit;
  int lengthToRead = 3*n;
  for(int i = 0; i < lengthToRead; i++) {
    bit = v.back();
    v.pop_back();
    data <<= 1;
    if(bit) {
      data.bits.set(0);
    }
  }
  return data;
}

from_PB_compressed compressTranscript(from_PB originalTranscript) {
  from_PB_compressed compressedTranscript;
  compressedTranscript.depth_ = originalTranscript.depth_;
  compressedTranscript.L_shares_recv = originalTranscript.L_shares_recv;
  compressedTranscript.R_shares_recv = originalTranscript.R_shares_recv;

  compressedTranscript.bit_L_shares_recv = originalTranscript.bit_L_shares_recv;
  compressedTranscript.bit_R_shares_recv = originalTranscript.bit_R_shares_recv;
  
  compressedTranscript.seed0L_encrypt.resize(originalTranscript.seed0L_encrypt.size());
  compressedTranscript.seed0R_encrypt.resize(originalTranscript.seed0R_encrypt.size());
  compressedTranscript.seed1L_encrypt.resize(originalTranscript.seed1L_encrypt.size());
  compressedTranscript.seed1R_encrypt.resize(originalTranscript.seed1R_encrypt.size());
  for(ssize_t i = 0; i < originalTranscript.depth_ + 1; i++) {
    compressedTranscript.seed0L_encrypt[i].resize(originalTranscript.seed0L_encrypt[i].size());
    compressedTranscript.seed0R_encrypt[i].resize(originalTranscript.seed0R_encrypt[i].size());
    compressedTranscript.seed1L_encrypt[i].resize(originalTranscript.seed1L_encrypt[i].size());
    compressedTranscript.seed1R_encrypt[i].resize(originalTranscript.seed1R_encrypt[i].size());
    for(unsigned int j = 0; j < rounds + 1; j++) {
      compressedTranscript.seed0L_encrypt[i][j] = originalTranscript.seed0L_encrypt[i][j];
      compressedTranscript.seed0R_encrypt[i][j] = originalTranscript.seed0R_encrypt[i][j];
      compressedTranscript.seed1L_encrypt[i][j] = originalTranscript.seed1L_encrypt[i][j];
      compressedTranscript.seed1R_encrypt[i][j] = originalTranscript.seed1R_encrypt[i][j];
    }
  }

  for(ssize_t i = 0; i < originalTranscript.depth_ + 1; i++) {
    for(unsigned int j = 0; j < 4; j++) {
      addToBitVector(compressedTranscript.blinds_recv, originalTranscript.blinds_recv[i][j], 1);
    }
  }
  // compressedTranscript.blinds_recv.resize(originalTranscript.blinds_recv.size());
  compressedTranscript.bit_blinds_recv.resize(originalTranscript.bit_blinds_recv.size());
  compressedTranscript.next_bit_L_recv.resize(originalTranscript.next_bit_L_recv.size());
  compressedTranscript.next_bit_R_recv.resize(originalTranscript.next_bit_R_recv.size());
  compressedTranscript.next_bit_L_recv2.resize(originalTranscript.next_bit_L_recv2.size());
  compressedTranscript.next_bit_R_recv2.resize(originalTranscript.next_bit_R_recv2.size());
  for(ssize_t i = 0; i < originalTranscript.depth_ + 1; i++) {
    // compressedTranscript.blinds_recv[i].resize(originalTranscript.blinds_recv[i].size());
    compressedTranscript.bit_blinds_recv[i].resize(originalTranscript.bit_blinds_recv[i].size());
    compressedTranscript.next_bit_L_recv[i].resize(originalTranscript.next_bit_L_recv[i].size());
    compressedTranscript.next_bit_R_recv[i].resize(originalTranscript.next_bit_R_recv[i].size());
    compressedTranscript.next_bit_L_recv2[i].resize(originalTranscript.next_bit_L_recv2[i].size());
    compressedTranscript.next_bit_R_recv2[i].resize(originalTranscript.next_bit_R_recv2[i].size());
    for(unsigned int j = 0; j < 4; j++) {
      // compressedTranscript.blinds_recv[i][j] = originalTranscript.blinds_recv[i][j];
      compressedTranscript.bit_blinds_recv[i][j] = originalTranscript.bit_blinds_recv[i][j];

      compressedTranscript.next_bit_R_recv[i][j] = originalTranscript.next_bit_R_recv[i][j];
      compressedTranscript.next_bit_L_recv[i][j] = originalTranscript.next_bit_L_recv[i][j];

      compressedTranscript.next_bit_R_recv2[i][j] = originalTranscript.next_bit_R_recv2[i][j];
      compressedTranscript.next_bit_L_recv2[i][j] = originalTranscript.next_bit_L_recv2[i][j];
    }
  }
  
  return compressedTranscript;
}

from_PB decompressTranscript(from_PB_compressed compressedTranscript) {
  from_PB decompressedTranscript(compressedTranscript.depth_, rounds);
  // decompressedTranscript.depth_ = compressedTranscript.depth_;
  decompressedTranscript.L_shares_recv = compressedTranscript.L_shares_recv;
  decompressedTranscript.R_shares_recv = compressedTranscript.R_shares_recv;

  decompressedTranscript.bit_L_shares_recv = compressedTranscript.bit_L_shares_recv;
  decompressedTranscript.bit_R_shares_recv = compressedTranscript.bit_R_shares_recv;
  

  for(ssize_t i = 0; i < compressedTranscript.depth_ + 1; i++) {
    for(unsigned int j = 0; j < rounds + 1; j++) {
      decompressedTranscript.seed0L_encrypt[i][j] = compressedTranscript.seed0L_encrypt[i][j];
      decompressedTranscript.seed0R_encrypt[i][j] = compressedTranscript.seed0R_encrypt[i][j];
      decompressedTranscript.seed1L_encrypt[i][j] = compressedTranscript.seed1L_encrypt[i][j];
      decompressedTranscript.seed1R_encrypt[i][j] = compressedTranscript.seed1R_encrypt[i][j];
    }
  }

  decompressedTranscript.blinds_recv.resize(compressedTranscript.depth_ + 1);
  for(ssize_t i = compressedTranscript.depth_; i >= 0; i--) {
    decompressedTranscript.blinds_recv[i].resize(4);
    for(int j = 3; j >= 0; j--) {
      decompressedTranscript.blinds_recv[i][j] = extractFromBitVector(compressedTranscript.blinds_recv, 1);
      // std::cout << "Extraction successful\n" << "i = " << i << '\n' << "j = " << j << '\n';
      // decompressedTranscript.blinds_recv[i][j] = b;
      // std::cout << "Access successful\n";
    }
  }

  for(ssize_t i = 0; i < compressedTranscript.depth_ + 1; i++) {
    for(unsigned int j = 0; j < 4; j++) {
      // decompressedTranscript.blinds_recv[i][j] = compressedTranscript.blinds_recv[i][j];
      decompressedTranscript.bit_blinds_recv[i][j] = compressedTranscript.bit_blinds_recv[i][j];

      decompressedTranscript.next_bit_R_recv[i][j] = compressedTranscript.next_bit_R_recv[i][j];
      decompressedTranscript.next_bit_L_recv[i][j] = compressedTranscript.next_bit_L_recv[i][j];

      decompressedTranscript.next_bit_R_recv2[i][j] = compressedTranscript.next_bit_R_recv2[i][j];
      decompressedTranscript.next_bit_L_recv2[i][j] = compressedTranscript.next_bit_L_recv2[i][j];
    }
  }

  return decompressedTranscript;
}


bool compareTranscripts(from_PB t1, from_PB t2) {
  bool result = (t1.depth_ == t2.depth_)
      && (t1.L_shares_recv.bits == t2.L_shares_recv.bits)
      && (t1.R_shares_recv.bits == t2.R_shares_recv.bits)
      && (t1.bit_L_shares_recv == t2.bit_L_shares_recv)
      && (t1.bit_R_shares_recv == t2.bit_R_shares_recv);
  // int rounds = t1.seed0L_encrypt[0].size();
  for(size_t k = 0; (k < t1.depth_ + 1) && result; k++) {
    for(int i = 0; i < 4 && result; i++) {
      result = result && (t1.blinds_recv[k][i].bits == t2.blinds_recv[k][i].bits);
      result = result && (t1.bit_blinds_recv[k][i] == t2.bit_blinds_recv[k][i]);

      result = result && (t1.next_bit_R_recv[k][i] == t2.next_bit_R_recv[k][i]);
      result = result && (t1.next_bit_L_recv[k][i] == t2.next_bit_L_recv[k][i]);

      result = result && (t1.next_bit_R_recv2[k][i] == t2.next_bit_R_recv2[k][i]);
      result = result && (t1.next_bit_L_recv2[k][i] == t2.next_bit_L_recv2[k][i]);
    }

    for(int i = 0; (i < rounds + 1) && result; i++) {
      result = result && (t1.seed0L_encrypt[k][i].bits == t2.seed0L_encrypt[k][i].bits);
      result = result && (t1.seed0R_encrypt[k][i].bits == t2.seed0R_encrypt[k][i].bits);
      result = result && (t1.seed1L_encrypt[k][i].bits == t2.seed1L_encrypt[k][i].bits);
      result = result && (t1.seed1R_encrypt[k][i].bits == t2.seed1R_encrypt[k][i].bits);
    }
  }
  return result;
}
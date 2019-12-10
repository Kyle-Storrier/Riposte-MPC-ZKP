struct  P2_root_layer
{

 block   blinds0[4], blinds1[4];
 block   gamma[4];
 bool    bit_blinds0[4], bit_blinds1[4];
 bool    next_bit_gamma[4];
};


struct  P2_middle_layer
{
  block blind0[depth][rounds], blind1[depth][rounds]; 
  block gamma0[depth][rounds], gamma1[depth][rounds];
  block blinds0[depth][4], blinds1[depth][4];
  block gamma[depth][4];

  bool bit_blinds0[depth][4], bit_blinds1[depth][4];   
  bool next_bit_gamma[depth][4];
};

struct  PB_root_layer
{
 //block L_shares_sent; // This is temporary
 //block R_shares_sent; // This is temporary
 
 block L_shares_recv;
 block R_shares_recv;

 block gamma[4];
 //block blinds_sent[4]; 
 block blinds_recv[4];

 bool gamma_bit[4];

 //bool bit_blinds_sent[4];
 bool bit_blinds_recv[4];
 
 bool bit_L_shares_sent; 
 bool bit_R_shares_sent;
 bool bit_L_shares_recv; 
 bool bit_R_shares_recv; 

 bool next_bit_L_sent[4], next_bit_L_recv[4];
 bool next_bit_R_sent[4], next_bit_R_recv[4];

};

struct PB_middle_layer
{
  block L_shares_sent[depth];
 
  block L_shares_recv[depth];
 
  block R_shares_sent[depth];
 
  block R_shares_recv[depth];
  
  block seedL[depth], seedR[depth];
  
  block gamma0[depth][rounds], gamma1[depth][rounds];
  //block blinds_sent[depth][4];
  
  block blinds_recv[depth][4];
  block gamma_middle[depth][4];

  //bool bit_blinds_sent[depth][4];
  bool bit_blinds_recv[depth][4];
  bool bit_L_shares_sent[depth], bit_R_shares_sent[depth];
  bool bit_L_shares_recv[depth], bit_R_shares_recv[depth]; 
 
  bool next_bit_L_sent[depth][4], next_bit_L_recv[depth][4];
  bool next_bit_R_sent[depth][4], next_bit_R_recv[depth][4];
  
  bool gamma_bit[depth][4];
  
  std::vector<block> TL[depth+1], TL_other[depth+1];
  std::vector<block> TR[depth+1], TR_other[depth+1];
  std::vector<block> TL2[depth+1], TL_other2[depth+1];
  std::vector<block> TR2[depth+1], TR_other2[depth+1];
};

struct  P2_Transcripts
{
  P2_root_layer   P2_root;  
  P2_middle_layer P2_middle;
};

struct  PB_Transcripts
{
  PB_root_layer    PB_root;  
  PB_middle_layer  PB_middle;
};


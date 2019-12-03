struct  P2_root_layer
{

 block   blinds0[4], blinds1[4];
 block   gamma[4];
 bool    bit_blinds0[4], bit_blinds1[4];
};


struct  P2_middle_layer
{
  block blinds0[depth][4], blinds1[depth][4];
  block gamma[depth][4];
  bool bit_blinds0[depth][4], bit_blinds1[depth][4];   
};

struct  PB_root_layer
{
 block L_shares_sent, L_shares_recv;
 block R_shares_sent, R_shares_recv;

 block   gamma[4];
 block   blinds_sent[4], blinds_recv[4];

 bool    bit_blinds_sent[4], bit_blinds_recv[4];
};

struct PB_middle_layer
{
  std::vector<block> T[depth];
  block   blinds_sent[depth][4],     blinds_recv[depth][4];
  bool    bit_blinds_sent[depth][4], bit_blinds_recv[depth][4];
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


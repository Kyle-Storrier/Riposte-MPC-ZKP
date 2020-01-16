// #include "simulator.h"


//   void Simulator::Gen_Blinds(block_t blinds0[rounds], block_t blinds1[rounds], block_t gamma[2][rounds])
//   {
//     block_t rand[rounds];
//     for (unsigned r = 0; r < rounds; ++r)
//      {        
//         blinds0[r] = P0rand.next_block();
//         blinds1[r] = P1rand.next_block();
//         rand[r]    = P2rand.next_block();

//         const block_t tmp1 = ((blinds0[r] >> 1) & blinds1[r]) ^ ((blinds1[r] >> 1) & blinds0[r]);
//         const block_t tmp2 = ((blinds0[r] >> 2) & blinds1[r]) ^ ((blinds1[r] >> 2) & blinds0[r]);
    
//         const block_t bc = (tmp1 << 2) & maska;
//         const block_t ac = (tmp2 << 1) & maskb;
//         const block_t ab = (tmp1 >> 1) & maskc;
    
//         gamma[0][r] = (bc | ac | ab) ^ rand[r];
//         gamma[1][r] = rand[r];
//      }
//   }

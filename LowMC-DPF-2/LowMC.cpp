#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <bsd/stdlib.h>

#include "LowMC.h"



/////////////////////////////
//     LowMC functions     //
/////////////////////////////

block LowMC::encrypt(const block & message) {
    block c = message ^ roundkeysXORconstants[0];

     for (unsigned r = 1; r <= rounds; ++r) 
     {
         c = Substitution(c);
         c = MultiplyWithGF2Matrix(LinMatrices[r-1], c, roundkeysXORconstants[r]);
     }
    return c;
}

std::vector<block> LowMC::encrypt_MPC_verify(const block & message, const std::vector<block> c2, const block blind[rounds], block gamma[rounds], const bool P) {
    std::vector<block> c(rounds+1);
    block tmp = P ? message : message ^ roundkeysXORconstants[0];

    for (unsigned r = 1; r <= rounds; ++r) {
        c[r-1] = tmp ^ blind[r-1];
        tmp = Substitution_MPC(tmp, c2[r-1], blind[r-1], gamma[r-1]);
        if(P)   tmp = MultiplyWithGF2Matrix(LinMatrices[r-1], tmp);
        if(!P)  tmp = MultiplyWithGF2Matrix(LinMatrices[r-1], tmp, roundkeysXORconstants[r]);
    }
    
    c[rounds] = tmp;

    return c;
}

std::pair<std::vector<block>,std::vector<block>> LowMC::encrypt_MPC_proof(const block & m0, const block & m1, const block blind0[rounds], const block blind1[rounds], block gamma[2][rounds]) {
    std::vector<block> c0(rounds + 1);
    std::vector<block> c1(rounds + 1);
    block tmp0 = m0;
    block tmp1 = m1 ^ roundkeysXORconstants[0];

    for (unsigned r = 1; r <= rounds; ++r) 
    {
        c0[r-1] = tmp0 ^ blind0[r-1];
        c1[r-1] = tmp1 ^ blind1[r-1];

        tmp0 = Substitution_MPC(tmp0, c1[r-1], blind0[r-1], gamma[0][r-1]);
        tmp1 = Substitution_MPC(tmp1, c0[r-1], blind1[r-1], gamma[1][r-1]);
        tmp0 = MultiplyWithGF2Matrix(LinMatrices[r-1], tmp0, roundkeysXORconstants[r]);
        tmp1 = MultiplyWithGF2Matrix(LinMatrices[r-1], tmp1);
    }


    block reconstructed = tmp0 ^ tmp1;
    c0[rounds] = tmp0;
    c1[rounds] = tmp1;
   // std::cout << "m0 ^ m1 = " << (m0 ^ m1) << std::endl;
   // std::cout << "(encrypt_MPC_proof) of (" <<  "--"  <<  ") : " << reconstructed << std::endl << std::endl;
    return std::make_pair(c0, c1);
}

void LowMC::runP2(AES_KEY& aeskey, block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds]) {
  //  GenBlinds(aeskey, blinds0, blinds1, gamma);
}


block LowMC::Substitution(const block & message) {
    const block srli1 = (message >> 1) & maskbc;
    const block srli2 = (message >> 2) & maskc;

    const block tmp = message & srli1;
    const block bc = (tmp << 2) & maska;
    const block ac = (message & srli2) << 1;
    const block ab = (tmp >> 1) & maskc;

    return (bc | ac | ab) ^ message ^ srli1 ^ srli2;
}

block LowMC::Substitution_MPC(const block & message, const block & message2, const block & blind, block  gamma) {
    
    const block srli1 = (message >> 1) & maskbc;
    const block srli2 = (message >> 2) & maskc;


    const block message3 = message ^ message2;
    const block tmp = (message3 & srli1) ^ (blind & (message2 >> 1));
    const block bc = (tmp << 2) & maska;
    const block ac = (((message3 & srli2) ^ (blind & (message2 >> 2))) << 1) & maskb;
    //const block ac = ((message3 ^ (blind & (message2 >> 2))) & srli2) << 1;
    const block ab = (tmp >> 1) & maskc;


    return (bc | ac | ab) ^ message ^ srli1 ^ srli2 ^ gamma;
}

void LowMC::GenBlinds(AES_KEY& aeskey, block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds]) {
   
     block rand[rounds];

    __m128i seed0; 
    __m128i seed1;
    seed0[0] = 21; seed0[1] = 21;
    seed1[0] = 2; seed1[1] = 2;
    __m128i seed2 = seed0 ^ seed1;
    size_t buflen =  rounds * sizeof(block);    
    PRG(aeskey, seed0, (__m128i *)blinds0, buflen / sizeof(__m128i));
    PRG(aeskey, seed1, (__m128i *)blinds1, buflen / sizeof(__m128i));
    PRG(aeskey, seed2, (__m128i *)rand, buflen / sizeof(__m128i));
    

 
    for (unsigned r = 0; r < rounds; ++r)
    {
        const block tmp1 = ((blinds0[r] >> 1) & blinds1[r]) ^ ((blinds1[r] >> 1) & blinds0[r]);
        const block tmp2 = ((blinds0[r] >> 2) & blinds1[r]) ^ ((blinds1[r] >> 2) & blinds0[r]);
    
        const block bc = (tmp1 << 2) & maska;
        const block ac = (tmp2 << 1) & maskb;
        const block ab = (tmp1 >> 1) & maskc;
    
        gamma[0][r] = (bc | ac | ab) ^ rand[r];
        gamma[1][r] = rand[r];// ^ roundkeysXORconstants[r+1];
    }
}

block LowMC::MultiplyWithGF2Matrix
        (const std::vector<block> & matrix, const block & message, const block & initial_value) {
    block temp = initial_value;

    uint64_t bitset;
    for (size_t k = 0; k < 4; ++k) {
        bitset = static_cast<__m256i>(message)[k];
        while (bitset != 0) {
          uint64_t t = bitset & -bitset;
          int i = k * 64 + __builtin_ctzl(bitset);
          temp = _mm256_xor_si256(temp, matrix[i]);
          bitset ^= t;
        }
    }

    return temp;
}

block LowMC::MultiplyWithGF2Matrix_Key
(const std::vector<keyblock> matrix, const keyblock k) {
block temp = 0;
for (unsigned i = 0; i < blocksize; ++i) {
temp[i] = (k & matrix[i]).parity();
}
return temp;
}

void LowMC::keyschedule () {
    for (unsigned r = 0; r <= rounds; ++r) {
        roundkeysXORconstants[r] ^= MultiplyWithGF2Matrix_Key (KeyMatrices[r], key);
    }
    return;
}


void LowMC::instantiate_LowMC () {
    // Create LinMatrices and invLinMatrices
    LinMatrices.clear();
    invLinMatrices.clear();
    for (unsigned r = 0; r < rounds; ++r) {
        //std::cout << "r = " << r << " rounds = " << rounds << std::endl;
        // Create matrix
        std::vector<block> mat;
        // Fill matrix with random bits
        do {
            mat.clear();
            for (unsigned i = 0; i < blocksize; ++i) {
                mat.push_back( getrandblock () );
            }

        // Repeat if matrix is not invertible
        } while ( rank_of_Matrix(mat) != blocksize );
        LinMatrices.push_back(mat);
        invLinMatrices.push_back(invert_Matrix (LinMatrices.back()));
    }

    // Create roundconstants
    roundkeysXORconstants.clear();
    roundkeysXORconstants.push_back(0);
    for (unsigned r = 0; r < rounds; ++r) {
        roundkeysXORconstants.push_back( getrandblock () );
    }
    // Create KeyMatrices
    KeyMatrices.clear();
    for (unsigned r = 0; r <= rounds; ++r) {
        // Create matrix
        std::vector<keyblock> mat;
        // Fill matrix with random bits
        do {
            mat.clear();
            for (unsigned i = 0; i < blocksize; ++i) {
                mat.push_back( getrandkeyblock () );
            }
        // Repeat if matrix is not of maximal rank
        } while ( rank_of_Matrix_Key(mat) < std::min(blocksize, keysize) );
        KeyMatrices.push_back(mat);
    }
    return;
}


/////////////////////////////
// Binary matrix functions //
/////////////////////////////


unsigned LowMC::rank_of_Matrix (const std::vector<block> matrix) {
    std::vector<block> mat; //Copy of the matrix 
    for (auto u : matrix) {
        mat.push_back(u);
    }
    unsigned size = mat[0].size();
    //Transform to upper triangular matrix
    unsigned row = 0;
    for (unsigned col = 1; col <= size; ++col) {
        
        if ( !mat[row][size-col] ) {
            unsigned r = row;
            while (r < mat.size() && !mat[r][size-col]) {
                ++r;
            }


            if (r >= mat.size()) {
                continue;
            } else {
                auto temp = mat[row];
                mat[row] = mat[r];
                mat[r] = temp;
            }
        }
        for (unsigned i = row+1; i < mat.size(); ++i) {
            if ( mat[i][size-col] ) mat[i] ^= mat[row];
        }
        ++row;
        if (row == size) break;
    }
    return row;
}


unsigned LowMC::rank_of_Matrix_Key (const std::vector<keyblock> matrix) {
    std::vector<keyblock> mat; //Copy of the matrix 
    for (auto u : matrix) {
        mat.push_back(u);
    }
    unsigned size = mat[0].size();
    //Transform to upper triangular matrix
    unsigned row = 0;
    for (unsigned col = 1; col <= size; ++col) {
        if ( !mat[row][size-col] ) {
            unsigned r = row;
            while (r < mat.size() && !mat[r][size-col]) {
                ++r;
            }
            if (r >= mat.size()) {
                continue;
            } else {
                auto temp = mat[row];
                mat[row] = mat[r];
                mat[r] = temp;
            }
        }
        for (unsigned i = row+1; i < mat.size(); ++i) {
            if ( mat[i][size-col] ) mat[i] ^= mat[row];
        }
        ++row;
        if (row == size) break;
    }
    return row;
}


std::vector<block> LowMC::invert_Matrix (const std::vector<block> matrix) {
    //std::cout << "invert" << std::endl;
    std::vector<block> mat; //Copy of the matrix 
    for (auto u : matrix) {
        mat.push_back(u);
    }
    
   std::vector<block> invmat(blocksize, 0); //To hold the inverted matrix
   
   // for (unsigned i = 0; i < blocksize; ++i) {
   //      for (unsigned j = 0; j < blocksize; ++j) {
   //      invmat[i][j] = 0;
   //      }
   //  }
        //std::cout << "invert" << std::endl;
    for (unsigned i = 0; i < blocksize; ++i) {
        invmat[i][i] = 1;
    }

    unsigned size = mat[0].size();

    //Transform to upper triangular matrix

    unsigned row = 0;
    for (unsigned col = 0; col < size; ++col) {
        if ( !mat[row][col] ) {
            unsigned r = row+1;
            while (r < mat.size() && !mat[r][col]) {
                ++r;
            }
            if (r >= mat.size()) {
                continue;
            } else {
                auto temp = mat[row];
                mat[row] = mat[r];
                mat[r] = temp;
                temp = invmat[row]; 
                invmat[row] = invmat[r];
                invmat[r] = temp;
            }
        }
        for (unsigned i = row+1; i < mat.size(); ++i) {
            if ( mat[i][col] ) {
                mat[i] ^= mat[row];
                invmat[i] ^= invmat[row]; 
            }
        }

        ++row;
    }

  
    //Transform to identity matrix
    for (unsigned col = size; col > 0; --col) {
        for (unsigned r = 0; r < col-1; ++r) {
            if (mat[r][col-1]) {
                mat[r] ^= mat[col-1];
                invmat[r] ^= invmat[col-1];
            }
        }
    }
 // std::cout << "done: " << std::endl;
    return invmat;
}

///////////////////////
// Pseudorandom bits //
///////////////////////


block LowMC::getrandblock () {
    block tmp = 0;
    for (unsigned i = 0; i < blocksize; ++i)
    {
     tmp[i] = getrandbit ();
    }
    return tmp;
}

keyblock LowMC::getrandkeyblock () {
    keyblock tmp = 0;
    for (unsigned i = 0; i < keysize; ++i) tmp[i] = getrandbit ();
    return tmp;
}


// Uses the Grain LSFR as self-shrinking generator to create pseudorandom bits
// Is initialized with the all 1s state
// The first 160 bits are thrown away
bool LowMC::getrandbit () {
    static std::bitset<80> state; //Keeps the 80 bit LSFR state
    bool tmp = 0;
    //If state has not been initialized yet
    if (state.none ()) {
        state.set (); //Initialize with all bits set
        //Throw the first 160 bits away
        for (unsigned i = 0; i < 160; ++i) {
            //Update the state
            tmp =  state[0] ^ state[13] ^ state[23]
                       ^ state[38] ^ state[51] ^ state[62];
            state >>= 1;
            state[79] = tmp;
        }
    }
    //choice records whether the first bit is 1 or 0.
    //The second bit is produced if the first bit is 1.
    bool choice = false;
    do {
        //Update the state
        tmp =  state[0] ^ state[13] ^ state[23]
                   ^ state[38] ^ state[51] ^ state[62];
        state >>= 1;
        state[79] = tmp;
        choice = tmp;
        tmp =  state[0] ^ state[13] ^ state[23]
                   ^ state[38] ^ state[51] ^ state[62];
        state >>= 1;
        state[79] = tmp;
    } while (! choice);
    return tmp;
}




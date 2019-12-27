#ifndef __LowMC_h__
#define __LowMC_h__

#include "block.h"
#include "prg.h"

typedef blocks<__m128i> keyblock;
typedef blocks<__m256i> block;


// Size of the identity part in the Sbox layer

static const block mask = std::string("0100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100");
static const block maska = mask.shiftr_bits(identitysize-1);
static const block maskb = maska >> 1;
static const block maskc = maska >> 2;
static const block maskbc = maskb | maskc;



class LowMC {
public:
    LowMC() : key(0) {
        instantiate_LowMC();
        keyschedule();
    };

    block encrypt (const block & message);
    block encrypt_MPC (const block & message, const block blind[rounds], const block gamma[rounds], const bool P = 0);
    std::pair<std::vector<block>,std::vector<block>> encrypt_MPC_proof(const block & m0, const block & m1, const block blind0[rounds], const block blind1[rounds], block gamma[2][rounds]);
    std::vector<block> encrypt_MPC_verify(const block & message, const std::vector<block> c2, const block blind[rounds],  block gamma[rounds], const bool P);
    void runP2(AES_KEY& aeskey, block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds]);
    block Substitution (const block & message);
    block Substitution_MPC (const block & message, const block & message2, const block & blind0, block gamma);
 
private:
// LowMC private data members //
    std::vector<std::vector<block>> LinMatrices;
        // Stores the binary matrices for each round
    std::vector<std::vector<block>> invLinMatrices;
        // Stores the inverses of LinMatrices
    std::vector<block> roundkeysXORconstants;
        // Stores the round constants
    const keyblock key;
        //Stores the master key
    std::vector<std::vector<keyblock>> KeyMatrices;
        // Stores the matrices that generate the round keys
    
// LowMC private functions //
    // block Substitution (const block & message);
    // block Substitution_MPC (const block & message, const block & message2, const block & blind, const block & gamma);
    void GenBlinds(AES_KEY& aeskey, block blinds0[rounds], block blinds1[rounds], block gamma[2][rounds]);

    block MultiplyWithGF2Matrix
        (const std::vector<block> & matrix, const block & message, const block & initial_value = 0);
        // For the linear layer
    block MultiplyWithGF2Matrix_Key
        (const std::vector<keyblock> matrix, const keyblock k);
        // For generating the round keys

    void keyschedule ();
        //Creates the round keys from the master key

    void instantiate_LowMC ();
        //Fills the matrices and roundconstants with pseudorandom bits 

// Binary matrix functions //
    unsigned rank_of_Matrix (const std::vector<block> matrix);
    unsigned rank_of_Matrix_Key (const std::vector<keyblock> matrix);
    std::vector<block> invert_Matrix (const std::vector<block> matrix);

// Random bits functions //
    block getrandblock ();
    keyblock getrandkeyblock ();
    bool  getrandbit ();

};

#endif

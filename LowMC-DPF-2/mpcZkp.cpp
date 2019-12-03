#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h> 
#include "dpf.h"


using namespace std;

block multiplicationConstants[2];

#define L 0
#define R 1

struct TranscriptEntry {
    bool generated; // True for data which is being randomly generated.
    int sender; // The party sending the data (same as receiver for generated data).
    int receiver; // The party receiving the data.
    block data; // The data.
    string label; // A label mainly intended for debugging and human readable output.
};

// The transcript of the MPC.
vector<TranscriptEntry> transcripts[3];

void generateData(int party, block data, string label) {
	TranscriptEntry t;
    t.generated = true;
    t.sender = party;
    t.receiver = party;
    t.data = data;
    t.label = label;

    transcripts[party].push_back(t);
}

void sendData(int sender, int receiver, block data, string label) {
    TranscriptEntry t;
    t.generated = false;
    t.sender = sender;
    t.receiver = receiver;
    t.data = data;
    t.label = label;

    transcripts[sender].push_back(t);
    transcripts[receiver].push_back(t);
}

// Perform Du-Attalah multiplication between a bit and an integer.
block* DuAttalahMultiplication(const block& x0, bool b0, const block& x1, bool b1)  {
    block* results = new block[2];

    // P2 generates random values D0, d0, D1, d1 and alpha
    block D0;
    arc4random_buf(&D0, sizeof(block));
    generateData(2, D0, "D_0");

    block D1;
    arc4random_buf(&D1, sizeof(block));
    generateData(2, D1, "D1");

    bool d0 = arc4random();
    generateData(2, d0, "d_0");

    bool d1 = arc4random();
    generateData(2, d1, "d_1");

    block alpha;
    arc4random_buf(&alpha, sizeof(block));
    generateData(2, alpha, "alpha");

    // P2 calculates c0 and c1.
    block c0 = (D0 & multiplicationConstants[d1]) ^ alpha; // = (D0 * d1) ^ alpha
    block c1 = (D1 & multiplicationConstants[d0]) ^ alpha; // = D1 * d0 ^ alpha

    // P2 sends D0, d0, c0 to P0
    sendData(2, 0, D0, "D_0");
    sendData(2, 0, d0, "d_0");
    sendData(2, 0, c0, "c_0 = D_0 * d_1 + alpha");
    // P2 sends D1, d1, and c1 to P1
    sendData(2, 1, D1, "D_1");
    sendData(2, 1, d1, "d_1");
    sendData(2, 1, c1, "c_1 = D_1 * d_0 - alpha");

    // P0 Computes X0 = x0 + D0 and Y0 = b0 + d0 and sends them to P1
    block X0 = x0 ^ D0;
    bool Y0 = b0 ^ d0;

    sendData(0, 1, X0, "X_0 = x_0 + D_0");
    sendData(0, 1, Y0, "Y_0 = b_0 + d_0");

    // P1 computes X1 = x1 + D1 and Y1 = b1 + d1 and sends them to P0
    block X1 = x1 ^ D1;
    bool Y1 = b1 ^ d1;

    sendData(1, 0, X1, "X_1 = x_1 + D_1");
    sendData(1, 0, Y1, "Y_1 = b_1 + d_1");

    results[0] = (x0 & multiplicationConstants[b0 ^ Y1]) ^ (multiplicationConstants[d0] & X1) ^ c0;
    results[1] = (x1 & multiplicationConstants[b1 ^ Y0]) ^ (multiplicationConstants[d1] & X0) ^ c1;
    return results;
}

void printTranscriptEntry(TranscriptEntry t) {
    if(t.generated) {
        cout << "Party " << t.sender << " generated " << t.label << " = " << t.data;
    } else {
        cout << "Party " << t.sender << " sent " << t.label << " = " << t.data << " to Party " << t.receiver;
    }
    cout << endl;
}

void printTranscript(int party) {
    for(unsigned long  i = 0; i < transcripts[party].size(); i++) {
        printTranscriptEntry(transcripts[party].at(i));
    }
}

template <typename KEY_TYPE, size_t nitems, typename __mX>
block mpcFirstStage(KEY_TYPE key, dpf_key<__mX, nitems> dpfkey[2], block seed0, uint8_t b0, block seed1, uint8_t b1, block* newSeed0Shares, uint8_t* newB0Shares, block* newSeed1Shares, uint8_t* newB1Shares) {
    block CW = dpfkey[0].cw[0];

    // P_0 computes L_0 || R_0 = G(seed_0) + b_0 * CW
    block expansion0[2];
    uint8_t exp0Bits[2];
    expand(key, seed0, expansion0, exp0Bits);
    expansion0[L] ^= multiplicationConstants[b0 ^ 1] & CW;
    expansion0[R] ^= multiplicationConstants[b0 ^ 1] & CW;

    // P_1 computes L_1 || R_1 = G(seed_1) + b_1 * CW
    block expansion1[2];
    uint8_t exp1Bits[2];
    expand(key, seed1, expansion1, exp1Bits);
    expansion1[L] ^= multiplicationConstants[b1 ^ 1] & CW;
    expansion1[R] ^= multiplicationConstants[b1 ^ 1] & CW;

    // cout << "Left Side: " << (expansion0[L] ^ expansion1[L]).b.to_string() << endl
    //     << "Right Side: " << (expansion0[R] ^ expansion1[R]).b.to_string() << endl << endl;

    // Shares of the direction value
    uint8_t directionShares[2];
    directionShares[0] = exp0Bits[R];
    directionShares[1] = exp1Bits[R] ^ dpfkey[0].t[0][R];

    // Use Du-Attalah multiplication to compute shares (1 - b) * L_0 + b * R_0 and to compute shares of b * L_0 + (1 - b) * R_0
    // TODO: Do this without the hardcoded 0 shares.
    
    block* bL0Shares    = DuAttalahMultiplication(expansion0[L], directionShares[0], 0, directionShares[1]); // L_0 * direction
    block* bNotL0Shares = DuAttalahMultiplication(expansion0[L], directionShares[0], 0, directionShares[1] ^ 1); // L_0 * (direction - 1)
    block* bL1Shares    = DuAttalahMultiplication(0, directionShares[0], expansion1[L], directionShares[1]); // L_1 * direction
    block* bNotL1Shares = DuAttalahMultiplication(0, directionShares[0], expansion1[L], directionShares[1] ^ 1); // L_1 * (direction - 1)
    block* bR0Shares    = DuAttalahMultiplication(expansion0[R], directionShares[0], 0, directionShares[1]); // R_0 * direction
    block* bNotR0Shares = DuAttalahMultiplication(expansion0[R], directionShares[0], 0, directionShares[1] ^ 1); // R_0 * (direction - 1)
    block* bR1Shares    = DuAttalahMultiplication(0, directionShares[0], expansion1[R], directionShares[1]); // R_1 * direction
    block* bNotR1Shares = DuAttalahMultiplication(0, directionShares[0], expansion1[R], directionShares[1] ^ 1); // R_1 * (direction - 1)

    // Compute shares of the corrected left side.
    block correctedLeftSideShares[2];
    correctedLeftSideShares[0] = bL0Shares[0] ^ bNotR0Shares[0] ^ bL1Shares[0] ^ bNotR1Shares[0];
    correctedLeftSideShares[1] = bL0Shares[1] ^ bNotR0Shares[1] ^ bL1Shares[1] ^ bNotR1Shares[1];

    // P0 and P1 reveal their shares of the corrected left side and calculate the result.
    sendData(0, 1, correctedLeftSideShares[0], "Corrected left side");
    sendData(1, 0, correctedLeftSideShares[1], "Corrected left side");
    // Compute the corrected left side
    block correctedLeftSide = correctedLeftSideShares[0] ^ correctedLeftSideShares[1];
    // = bL0Shares[0] ^ bNotR0Shares[0] ^ bL1Shares[0] ^ bNotR1Shares[0] ^ bL0Shares[1] ^ bNotR0Shares[1] ^ bL1Shares[1] ^ bNotR1Shares[1]
    // = (b * L0) ^ ((b-1) * R0) ^ (b * L1) ^ ((b-1) * R1) = L * b + R * (1 - b)

    // Set the new values for the seeds
    newSeed0Shares[0] = bNotL0Shares[0] ^ bR0Shares[0];
    newSeed0Shares[1] = bNotL0Shares[1] ^ bR0Shares[1];
    newSeed1Shares[0] = bNotL1Shares[0] ^ bR1Shares[0];
    newSeed1Shares[1] = bNotL1Shares[1] ^ bR1Shares[1];

    // Set the new values for the bits.
    // TODO: Do this properly using MPC.
    newB0Shares[0] = (exp0Bits[0] & directionShares[0]) ^ (exp0Bits[0] & (!directionShares[1]));
    newB0Shares[1] = (exp0Bits[1] & directionShares[0]) ^ (exp0Bits[1] & directionShares[1]) ^ (dpfkey[0].t[0][directionShares[0] ^ directionShares[1]] & b0);
    newB1Shares[0] = (exp1Bits[0] & directionShares[0]) ^ (exp1Bits[0] & (!directionShares[1]));
    newB1Shares[1] = (exp1Bits[0] & directionShares[0]) ^ (exp1Bits[0] & directionShares[1]) ^ (dpfkey[1].t[0][directionShares[0] ^ directionShares[1]] & b1);


    //printf("Expected: 0\nActual: %s\n", correctedLeftSide.b.to_string().c_str());
    return correctedLeftSide;
}

template <typename KEY_TYPE, size_t nitems, typename __mX>
block mpcInnerStage(KEY_TYPE key, dpf_key<__mX, nitems> dpfkey[2], ssize_t index, block* seed0Shares, uint8_t* b0Shares, block* seed1Shares, uint8_t* b1Shares) {
    block CW = dpfkey[0].cw[index];

    block expansion0Shares[2][2]; // First dimension is party 0 or 1, second dimension is L or R.
    block expansion1Shares[2][2]; // First dimension is party 0 or 1, second dimension is L or R.
    uint8_t expansion0BitShares[2][2];  // First dimension is party 0 or 1, second dimension is L or R.
    uint8_t expansion1BitShares[2][2];  // First dimension is party 0 or 1, second dimension is L or R.

    //**********************************************************************************************
    // TODO: Replace this code with the MPC expansion code.
    // P_0 computes L_0 || R_0 = G(seed_0) + b_0 * CW
    block expansion0Hidden[2];
    uint8_t exp0Bits[2];
    expand(key, seed0Shares[0] ^ seed0Shares[1], expansion0Hidden, exp0Bits);
    // cout << "Seed 0:\nLeft Expansion: " << (expansion0Hidden[L]).b.to_string() << endl << "Right Expansion: " << (expansion0Hidden[R]).b.to_string() << endl;

    // Secret share the expansion's left side.
    arc4random_buf(&expansion0Shares[0][L], sizeof(block));
    expansion0Shares[1][L] = expansion0Shares[0][L] ^ expansion0Hidden[L];
    expansion0BitShares[0][L] = arc4random() & 0x01;
    expansion0BitShares[1][L] = exp0Bits[L] ^ expansion0BitShares[0][L];

    // Secret share the expansion's right side.
    arc4random_buf(&expansion0Shares[0][R], sizeof(block));
    expansion0Shares[1][R] = expansion0Hidden[R] ^ expansion0Shares[0][R];
    expansion0BitShares[0][R] = arc4random() & 0x01;
    expansion0BitShares[1][R] = exp0Bits[R] ^ expansion0BitShares[0][R];
    //**********************************************************************************************


    expansion0Shares[0][L] ^= multiplicationConstants[b0Shares[0]] & CW;
    expansion0Shares[1][L] ^= multiplicationConstants[b0Shares[1] ^ 1] & CW;
    expansion0Shares[0][R] ^= multiplicationConstants[b0Shares[0]] & CW;
    expansion0Shares[1][R] ^= multiplicationConstants[b0Shares[1] ^ 1] & CW;


    //**********************************************************************************************
    // TODO: Replace this code with the MPC expansion code.
    block expansion1Hidden[2];
    uint8_t exp1BitsHidden[2];
    expand(key, seed1Shares[0] ^ seed1Shares[1], expansion1Hidden, exp1BitsHidden);

    // cout << "Seed 1:\nLeft Expansion: " << (expansion1Hidden[L]).b.to_string() << endl << "Right Expansion: " << (expansion1Hidden[R]).b.to_string() << endl;
    // Secret share the expansion's left side.
    arc4random_buf(&expansion1Shares[0][L], sizeof(block));
    expansion1Shares[1][L] = expansion1Hidden[L] ^ expansion1Shares[0][L];
    expansion1BitShares[0][L] = arc4random() & 0x01;
    expansion1BitShares[1][L] = exp1BitsHidden[L] ^ expansion1BitShares[0][L];

    // Secret share the expansion's right side.
    arc4random_buf(&expansion1Shares[0][R], sizeof(block));
    expansion1Shares[1][R] = expansion1Hidden[R] ^ expansion1Shares[0][R];
    expansion1BitShares[0][R] = arc4random() & 0x01;
    expansion1BitShares[1][R] = exp1BitsHidden[R] ^ expansion1BitShares[0][R];
    //**********************************************************************************************

    expansion1Shares[0][L] ^= multiplicationConstants[b1Shares[0]] & CW;
    expansion1Shares[1][L] ^= multiplicationConstants[b1Shares[1] ^ 1] & CW;
    expansion1Shares[0][R] ^= multiplicationConstants[b1Shares[0]] & CW;
    expansion1Shares[1][R] ^= multiplicationConstants[b1Shares[1] ^ 1] & CW;

    // cout << "Left Side: " << (expansion0Shares[0][L] ^ expansion0Shares[1][L] ^ expansion1Shares[0][L] ^ expansion1Shares[1][L]).b.to_string() << endl
    //     << "Right Side: " << (expansion0Shares[0][R] ^ expansion0Shares[1][R] ^ expansion1Shares[0][R] ^ expansion1Shares[1][R]).b.to_string() << endl << endl;

    // Shares of the direction value
    uint8_t directionShares[2];
    directionShares[0] = expansion0BitShares[0][R] ^ expansion1BitShares[0][R];
    directionShares[1] = expansion0BitShares[1][R] ^ expansion1BitShares[1][R] ^ dpfkey[0].t[index][R];

    cout << "Direction: " << (directionShares[0] ^ directionShares[1]) << endl;

    // Use Du-Attalah multiplication to compute shares (1 - b) * L_0 + b * R_0 and to compute shares of b * L_0 + (1 - b) * R_0
    block* bL0Shares = DuAttalahMultiplication(expansion0Shares[0][L], directionShares[0], expansion0Shares[1][L], directionShares[1]); // L_0 * direction

    // L_0 * (direction - 1) = L_0 * direction - L_0
    block bNotL0Shares[2];
    bNotL0Shares[0] = bL0Shares[0] ^ expansion0Shares[0][L];
    bNotL0Shares[1] = bL0Shares[1] ^ expansion0Shares[1][L];    

    block* bL1Shares = DuAttalahMultiplication(expansion1Shares[0][L], directionShares[0], expansion1Shares[1][L], directionShares[1]); // L_1 * direction

    // L_1 * (direction - 1) = L_1 * direction - L_1
    block bNotL1Shares[2];
    bNotL1Shares[0] = bL1Shares[0] ^ expansion1Shares[0][L];
    bNotL1Shares[1] = bL1Shares[1] ^ expansion1Shares[1][L];

    block* bR0Shares = DuAttalahMultiplication(expansion0Shares[0][R], directionShares[0], expansion0Shares[1][R], directionShares[1]); // R_0 * direction

    // R_0 * (direction - 1) = R_0 * direction - R_0
    block bNotR0Shares[2];
    bNotR0Shares[0] = bR0Shares[0] ^ expansion0Shares[0][R];
    bNotR0Shares[1] = bR0Shares[1] ^ expansion0Shares[1][R];

    block* bR1Shares = DuAttalahMultiplication(expansion1Shares[0][R], directionShares[0], expansion1Shares[1][R], directionShares[1]); // R_1 * direction

    // R_1 * (direction - 1) = R_1 * direction - R_1
    block bNotR1Shares[2];
    bNotR1Shares[0] = bR1Shares[0] ^ expansion1Shares[0][R];
    bNotR1Shares[1] = bR1Shares[1] ^ expansion1Shares[1][R];

    // Compute shares of the corrected left side.
    block correctedLeftSideShares[2];
    correctedLeftSideShares[0] = bL0Shares[0] ^ bNotR0Shares[0] ^ bL1Shares[0] ^ bNotR1Shares[0];
    correctedLeftSideShares[1] = bL0Shares[1] ^ bNotR0Shares[1] ^ bL1Shares[1] ^ bNotR1Shares[1];

    // P0 and P1 reveal their shares of the corrected left side and calculate the result.
    sendData(0, 1, correctedLeftSideShares[0], "Corrected left side");
    sendData(1, 0, correctedLeftSideShares[1], "Corrected left side");

    // Compute the corrected left side
    block correctedLeftSide = correctedLeftSideShares[0] ^ correctedLeftSideShares[1];
    // = bL0Shares[0] ^ bNotR0Shares[0] ^ bL1Shares[0] ^ bNotR1Shares[0] ^ bL0Shares[1] ^ bNotR0Shares[1] ^ bL1Shares[1] ^ bNotR1Shares[1]
    // = (b * L0) ^ ((b-1) * R0) ^ (b * L1) ^ ((b-1) * R1) = L * b + R * (1 - b)

    // Set the new values for the seeds
    seed0Shares[0] = bNotL0Shares[0] ^ bR0Shares[0];
    seed0Shares[1] = bNotL0Shares[1] ^ bR0Shares[1];
    seed1Shares[0] = bNotL1Shares[0] ^ bR1Shares[0];
    seed1Shares[1] = bNotL1Shares[1] ^ bR1Shares[1];

    // // Set the new values for the bits.
    // TODO: Do this properly using MPC.
    uint8_t newB0Shares[2];
    newB0Shares[0] = 0;
    newB0Shares[1] = expansion0BitShares[0][directionShares[0] ^ directionShares[1]]
        ^ expansion0BitShares[1][directionShares[0] ^ directionShares[1]]
        ^ (dpfkey[0].t[index][directionShares[0] ^ directionShares[1]] & (b0Shares[0] ^ b0Shares[1]));

    uint8_t newB1Shares[2];
    newB1Shares[0] = 0;
    newB1Shares[1] = expansion1BitShares[0][directionShares[0] ^ directionShares[1]]
        ^ expansion1BitShares[1][directionShares[0] ^ directionShares[1]]
        ^ (dpfkey[1].t[index][directionShares[0] ^ directionShares[1]] & (b1Shares[0] ^ b1Shares[1]));

    // Update the state for the next round.
    b0Shares[0] = newB0Shares[0];
    b0Shares[1] = newB0Shares[1];
    b1Shares[0] = newB1Shares[0];
    b1Shares[1] = newB1Shares[1];

   // printf("Expected: 0\nActual: %s\n", correctedLeftSide.b.to_string().c_str());

    return correctedLeftSide;
}

int main() {
    typedef __m256i __mX;
    LowMC key;


    // Define the constant all 1s and all 0s blocks used for multiplication by boolean values.
    block *allZeros = new block();
    block *allOnes = new block();
   
    /////////allOnes->set();
   
    multiplicationConstants[0] = *allZeros;
    multiplicationConstants[1] = *allOnes;

    // The prover must generate the DPF keys and then the zero-knowledge proof of its correctness.
    const size_t nitems   =  1ULL << 15; // TODO: Change this value to reflect the actual number of items.
    dpf_key<__mX, nitems> dpfkey[2] = { 0 }; 
    gen(key, 2 , dpfkey);
    constexpr size_t depth = dpf_key<__mX, nitems>::depth;

    // Apply layer one of the proof generation MPC protocol
    block seed0Shares[2];
    block seed1Shares[2];
    uint8_t bit0Shares[2];
    uint8_t bit1Shares[2];
    bool failed = false;
    block result;
    ssize_t index = 0;

    result = mpcFirstStage(key, dpfkey, dpfkey[0].root, get_lsb(dpfkey[0].root), dpfkey[1].root, get_lsb(dpfkey[1].root), seed0Shares, bit0Shares, seed1Shares, bit1Shares);
   
   ////////// failed = result.count() > 0;
   
    for(ssize_t i = 1; i < depth && !failed; i++) {
        // cout << endl
        //     << "Bit 0: " << (bit0Shares[0] ^ bit0Shares[1]) << endl
        //     << "Bit 1: " << (bit1Shares[0] ^ bit1Shares[1]) << endl
        //     << "Seed 0: " << (seed0Shares[0] ^ seed0Shares[1]).b.to_string() << endl
        //     << "Seed 1: " << (seed1Shares[0] ^ seed1Shares[1]).b.to_string() << endl
        //     << endl;
        result = mpcInnerStage(key, dpfkey, i, seed0Shares, bit0Shares, seed1Shares, bit1Shares);
        index = i;
        
        /////failed = result.count() > 0;
    }

    if(failed) {
        cout << "\nFailed at layer: " << index << endl;
    } else {
        cout << "\nValid" << endl;
    }
    
}

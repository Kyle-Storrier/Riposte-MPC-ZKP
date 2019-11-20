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
void mpcFirstStage(KEY_TYPE key, dpf_key<__mX, nitems> dpfkey[2], block seed0, uint8_t b0, block seed1, uint8_t b1, block* newSeed0Shares, uint8_t* newB0Shares, block* newSeed1Shares, uint8_t* newB1Shares) {
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

    cout << "Left Side: " << (expansion0[L] ^ expansion1[L]).b.to_string() << endl
        << "Right Side: " << (expansion0[R] ^ expansion1[R]).b.to_string() << endl << endl;

    // Shares of the direction value
    uint8_t directionShares[2];
    directionShares[0] = exp0Bits[R];
    directionShares[1] = exp1Bits[R] ^ dpfkey[0].t[0][R];

    // Use Du-Attalah multiplication to compute shares (1 - b) * L_0 + b * R_0 and to compute shares of b * L_0 + (1 - b) * R_0
    // TODO: Do this without the hardcoded 0 shares.
    block* bL0Shares = DuAttalahMultiplication(expansion0[L], directionShares[0], 0, directionShares[1]); // L_0 * direction
    block* bNotL0Shares = DuAttalahMultiplication(expansion0[L], directionShares[0], 0, directionShares[1] ^ 1); // L_0 * (direction - 1)
    block* bL1Shares = DuAttalahMultiplication(0, directionShares[0], expansion1[L], directionShares[1]); // L_1 * direction
    block* bNotL1Shares = DuAttalahMultiplication(0, directionShares[0], expansion1[L], directionShares[1] ^ 1); // L_1 * (direction - 1)
    block* bR0Shares = DuAttalahMultiplication(expansion0[R], directionShares[0], 0, directionShares[1]); // R_0 * direction
    block* bNotR0Shares = DuAttalahMultiplication(expansion0[R], directionShares[0], 0, directionShares[1] ^ 1); // R_0 * (direction - 1)
    block* bR1Shares = DuAttalahMultiplication(0, directionShares[0], expansion1[R], directionShares[1]); // R_1 * direction
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

    // // Set the new values for the bits.
    newB0Shares[0] = 0;
    newB0Shares[1] = exp0Bits[directionShares[0] ^ directionShares[1]] ^ (dpfkey[0].t[0][directionShares[0] ^ directionShares[1]] & b0);
    newB1Shares[0] = 0;
    newB1Shares[1] = exp1Bits[directionShares[0] ^ directionShares[1]] ^ (dpfkey[1].t[0][directionShares[0] ^ directionShares[1]] & b1);
    //  bit0 = exp0Bits[direction] ^ (dpfkey[0].t[i][direction] & bit0);
    //  bit1 = exp1Bits[direction] ^ (dpfkey[1].t[i][direction] & bit1);
    // newB0Shares[0] = newSeed0Shares[0].get_lsb();
    // newB0Shares[1] = newSeed0Shares[1].get_lsb();
    // newB1Shares[0] = newSeed1Shares[0].get_lsb();
    // newB1Shares[1] = newSeed1Shares[1].get_lsb();


    printf("Expected: 0\nActual: %s\n", correctedLeftSide.b.to_string().c_str());
}

template <typename KEY_TYPE, size_t nitems, typename __mX>
void mpcInnerStage(KEY_TYPE key, dpf_key<__mX, nitems> dpfkey[2], ssize_t index, block* seed0Shares, uint8_t* b0Shares, block* seed1Shares, uint8_t* b1Shares) {
    block CW = dpfkey[0].cw[index];

    block expansion0Shares[2][2]; // First dimension is party 0 or 1, second dimension is L or R.
    block expansion1Shares[2][2]; // First dimension is party 0 or 1, second dimension is L or R.
    uint8_t expansion0BitShares[2][2];  // First dimension is party 0 or 1, second dimension is L or R.
    uint8_t expansion1BitShares[2][2];  // First dimension is party 0 or 1, second dimension is L or R.

    //***********************************************
    // TODO: Replace this code with the MPC expansion code.
    // P_0 computes L_0 || R_0 = G(seed_0) + b_0 * CW
    block expansion0Hidden[2];
    uint8_t exp0Bits[2];
    expand(key, seed0Shares[0] ^ seed0Shares[1], expansion0Hidden, exp0Bits);
    // cout << "Seed 0:\nLeft Expansion: " << (expansion0Hidden[L]).b.to_string() << endl << "Right Expansion: " << (expansion0Hidden[R]).b.to_string() << endl;
    // Secret share the expansion's left side.
    // arc4random_buf(&expansion0Shares[0][L], sizeof(block));
    // expansion0Shares[1][L] = expansion0Shares[0][L] ^ expansion0Hidden[L];
    expansion0Shares[0][L] = 0;
    expansion0Shares[1][L] = expansion0Hidden[L];

    expansion0BitShares[0][L] = 0; //expansion0Shares[0][L].get_lsb();
    expansion0BitShares[1][L] = exp0Bits[L]; // expansion0Shares[1][L].get_lsb();
    // Secret share the expansion's right side.
    // arc4random_buf(&expansion0Shares[0][R], sizeof(block));
    // expansion0Shares[1][R] = expansion0Shares[0][R] ^ expansion0Hidden[R];
    expansion0Shares[0][R] = 0;
    expansion0Shares[1][R] = expansion0Hidden[R];

    expansion0BitShares[0][R] = 0; //expansion0Shares[0][R].get_lsb();
    expansion0BitShares[1][R] = exp0Bits[R]; //expansion0Shares[1][R].get_lsb();
    //***********************************************


    expansion0Shares[0][L] ^= multiplicationConstants[b0Shares[0]] & CW;
    expansion0Shares[1][L] ^= multiplicationConstants[b0Shares[1] ^ 1] & CW;
    expansion0Shares[0][R] ^= multiplicationConstants[b0Shares[0]] & CW;
    expansion0Shares[1][R] ^= multiplicationConstants[b0Shares[1] ^ 1] & CW;


    //***********************************************
    // TODO: Replace this code with the MPC expansion code.
    block expansion1Hidden[2];
    uint8_t exp1Bits[2];
    expand(key, seed1Shares[0] ^ seed1Shares[1], expansion1Hidden, exp1Bits);

    // cout << "Seed 1:\nLeft Expansion: " << (expansion1Hidden[L]).b.to_string() << endl << "Right Expansion: " << (expansion1Hidden[R]).b.to_string() << endl;
    // Secret share the expansion's left side.
    // arc4random_buf(&expansion1Shares[0][L], sizeof(block));
    // expansion1Shares[1][L] = expansion1Shares[0][L] ^ expansion1Hidden[L];
    expansion1Shares[0][L] = 0;
    expansion1Shares[1][L] = expansion1Hidden[L];

    expansion1BitShares[0][L] = 0; //expansion1Shares[0][L].get_lsb();
    expansion1BitShares[1][L] = exp1Bits[L]; //expansion1Shares[1][L].get_lsb();
    // Secret share the expansion's right side.
    // arc4random_buf(&expansion1Shares[0][R], sizeof(block));
    // expansion1Shares[1][R] = expansion1Shares[0][R] ^ expansion1Hidden[R];
    expansion1Shares[0][R] = 0;
    expansion1Shares[1][R] = expansion1Hidden[R];

    expansion1BitShares[0][R] = 0; //expansion1Shares[0][R].get_lsb();
    expansion1BitShares[1][R] = exp1Bits[R]; // expansion1Shares[1][R].get_lsb();
    //***********************************************

    expansion1Shares[0][L] ^= multiplicationConstants[b1Shares[0]] & CW;
    expansion1Shares[1][L] ^= multiplicationConstants[b1Shares[1] ^ 1] & CW;
    expansion1Shares[0][R] ^= multiplicationConstants[b1Shares[0]] & CW;
    expansion1Shares[1][R] ^= multiplicationConstants[b1Shares[1] ^ 1] & CW;

    cout << "Left Side: " << (expansion0Shares[0][L] ^ expansion0Shares[1][L] ^ expansion1Shares[0][L] ^ expansion1Shares[1][L]).b.to_string() << endl
        << "Right Side: " << (expansion0Shares[0][R] ^ expansion0Shares[1][R] ^ expansion1Shares[0][R] ^ expansion1Shares[1][R]).b.to_string() << endl << endl;

    // Shares of the direction value
    uint8_t directionShares[2];
    directionShares[0] = expansion0BitShares[0][R] ^ expansion1BitShares[0][R];
    directionShares[1] = expansion0BitShares[1][R] ^ expansion1BitShares[1][R] ^ dpfkey[0].t[index][R];

    cout << "Direction: " << (directionShares[0] ^ directionShares[1]) << "\tShould be: " << (exp1Bits[R] ^ exp0Bits[R] ^ dpfkey[0].t[index][R]) << endl;
    // Use Du-Attalah multiplication to compute shares (1 - b) * L_0 + b * R_0 and to compute shares of b * L_0 + (1 - b) * R_0
    // TODO: Do this without the hardcoded 0 shares.
    block* bL0Shares = DuAttalahMultiplication(expansion0Shares[0][L], directionShares[0], expansion0Shares[1][L], directionShares[1]); // L_0 * direction
    block* bNotL0Shares = DuAttalahMultiplication(expansion0Shares[0][L], directionShares[0], expansion0Shares[1][L], directionShares[1] ^ 1); // L_0 * (direction - 1)
    block* bL1Shares = DuAttalahMultiplication(expansion1Shares[0][L], directionShares[0], expansion1Shares[1][L], directionShares[1]); // L_1 * direction
    block* bNotL1Shares = DuAttalahMultiplication(expansion1Shares[0][L], directionShares[0], expansion1Shares[1][L], directionShares[1] ^ 1); // L_1 * (direction - 1)
    block* bR0Shares = DuAttalahMultiplication(expansion0Shares[0][R], directionShares[0], expansion0Shares[1][R], directionShares[1]); // R_0 * direction
    block* bNotR0Shares = DuAttalahMultiplication(expansion0Shares[0][R], directionShares[0], expansion0Shares[1][R], directionShares[1] ^ 1); // R_0 * (direction - 1)
    block* bR1Shares = DuAttalahMultiplication(expansion1Shares[0][R], directionShares[0], expansion1Shares[1][R], directionShares[1]); // R_1 * direction
    block* bNotR1Shares = DuAttalahMultiplication(expansion1Shares[0][R], directionShares[0], expansion1Shares[1][R], directionShares[1] ^ 1); // R_1 * (direction - 1)

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
    //  bit0 = exp0Bits[direction] ^ (dpfkey[0].t[i][direction] & bit0);
    //  bit1 = exp1Bits[direction] ^ (dpfkey[1].t[i][direction] & bit1);
    b0Shares[0] = 0;
    b0Shares[1] = exp0Bits[directionShares[0] ^ directionShares[1]] ^ (dpfkey[0].t[index][directionShares[0] ^ directionShares[1]] & (b0Shares[0] ^ b0Shares[1]));
    b1Shares[0] = 0;
    b1Shares[1] = exp1Bits[directionShares[0] ^ directionShares[1]] ^ (dpfkey[1].t[index][directionShares[0] ^ directionShares[1]] & (b1Shares[0] ^ b1Shares[1]));
    // b0Shares[0] = seed0Shares[0].get_lsb();
    // b0Shares[1] = seed0Shares[1].get_lsb();
    // b1Shares[0] = seed1Shares[0].get_lsb();
    // b1Shares[1] = seed1Shares[1].get_lsb();

    printf("Expected: 0\nActual: %s\n", correctedLeftSide.b.to_string().c_str());
}

int main() {
    typedef __m256i __mX;
    LowMC key(1);


    // Define the constant all 1s and all 0s blocks used for multiplication by boolean values.
    block *allZeros = new block();
    block *allOnes = new block();
    allOnes->set();
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
    mpcFirstStage(key, dpfkey, dpfkey[0].root, dpfkey[0].root.get_lsb(), dpfkey[1].root, dpfkey[1].root.get_lsb(), seed0Shares, bit0Shares, seed1Shares, bit1Shares);
    for(ssize_t i = 1; i < depth; i++) {
        cout << endl
            << "Bit 0: " << (bit0Shares[0] ^ bit0Shares[1]) << endl
            << "Bit 1: " << (bit1Shares[0] ^ bit1Shares[1]) << endl
            << "Seed 0: " << (seed0Shares[0] ^ seed0Shares[1]).b.to_string() << endl
            << "Seed 1: " << (seed1Shares[0] ^ seed1Shares[1]).b.to_string() << endl
            << endl;
        mpcInnerStage(key, dpfkey, i, seed0Shares, bit0Shares, seed1Shares, bit1Shares);
    }
    
    cout << endl << endl << " --------------------------------------------------------------------------------\n\n";

    // Initialize the state information based on the key.
    block seed0 = dpfkey[0].root;
    block seed1 = dpfkey[1].root;
    uint8_t bit0 = dpfkey[0].root.get_lsb();
    uint8_t bit1 = dpfkey[1].root.get_lsb();
    block CW;
    uint8_t direction;
    
    for(ssize_t i = 0; i < 2; i++) {
        printf("Bit 0: %d\nBit 1: %d\n", bit0, bit1);
        cout
            << "Seed 0: " << seed0.b.to_string() << endl
            << "Seed 1: " << seed1.b.to_string() << endl
            << endl;
        CW = dpfkey[1].cw[i];

        // Expand seed 0.
        block expansion0[2];
        uint8_t exp0Bits[2];
        expand(key, seed0, expansion0, exp0Bits);
        // cout << "Seed0\nLeft expansion: " << expansion0[L].b.to_string() << endl << "Right expansion: " << expansion0[R].b.to_string() << endl;
        expansion0[L] ^= multiplicationConstants[bit0 ^ 1] & CW;
        expansion0[R] ^= multiplicationConstants[bit0 ^ 1] & CW;

        // Expand seed 1.
        block expansion1[2];
        uint8_t exp1Bits[2];
        expand(key, seed1, expansion1, exp1Bits);
        // cout << "Seed1\nLeft expansion: " << expansion1[L].b.to_string() << endl << "Right expansion: " << expansion1[R].b.to_string() << endl;
        expansion1[L] ^= multiplicationConstants[bit1 ^ 1] & CW;
        expansion1[R] ^= multiplicationConstants[bit1 ^ 1] & CW;

        direction = dpfkey[0].t[i][R] ^ exp0Bits[R] ^ exp1Bits[R]; // bit0 ^ bit1 ^ 1;
        // printf("Layer: %ld\nDirection: %d\nB0: %d\nB1: %d\n\n", i, direction, bit0, bit1);

        // Conditional on the bit, flip the right and left sides.
        block correctedLeft[2];
        correctedLeft[0] = (expansion0[L] & multiplicationConstants[direction]) 
            ^ (expansion0[R] & multiplicationConstants[direction ^ 1]);
        correctedLeft[1] = (expansion1[L] & multiplicationConstants[direction])
            ^ (expansion1[R] & multiplicationConstants[direction ^ 1]);
        block correctedRight[2];
        correctedRight[0] = (expansion0[R] & multiplicationConstants[direction])
            ^ (expansion0[L] & multiplicationConstants[direction ^ 1]);
        correctedRight[1] = (expansion1[R] & multiplicationConstants[direction])
            ^ (expansion1[L] & multiplicationConstants[direction ^ 1]);

        cout << "Expected: 0\nActual: " << (correctedLeft[0] ^ correctedLeft[1]).b.to_string() << endl
            << "Right Side: " << (correctedRight[0] ^ correctedRight[1]).b.to_string() << endl << endl;
        if((correctedLeft[0] ^ correctedLeft[1]).count() != 0) {
            cout << "Failed on layer " << i << endl;
            break;
        }

        if(i == depth - 1) {
            cout << endl << "Final Verification\n"
                << "Expected: 1\nActual: " << (correctedRight[0] ^ correctedRight[1]).b.to_string() << endl;
        }

        // Update the state based on the generated values and the key.
        bit0 = exp0Bits[direction] ^ (dpfkey[0].t[i][direction] & bit0);
        bit1 = exp1Bits[direction] ^ (dpfkey[1].t[i][direction] & bit1);
        seed0 = correctedRight[0];
        seed1 = correctedRight[1];
    }

    // cout << "Party 0:" << endl;
    // printTranscript(0);
    // cout << endl;
    // cout << "Party 1:" << endl;
    // printTranscript(1);
    // cout << endl;
    // cout << "Party 2:" << endl;
    // printTranscript(2);

    // delete[] result;
    // delete multiplicationConstants[0];
    // delete multiplicationConstants[1];
}


// Demo: Multiply 0x43 by the bit 1 using Du-Attalah Muliplication, and print out the complete set of transcripts.
// block x0  = 0xFD;
// block x1 = 0xBE;
// bool b0 = 1;
// bool b1 = 0;
// block expectedResult = ((x0 ^ x1) & multiplicationConstants[b0 ^ b1]);
// block *result = DuAttalahMultiplication(x0, b0, x1, b1);
// cout << "Expected: " << expectedResult.to_ulong() << endl
//     << "Actual: " << (result[0] ^ result[1]).to_ulong() << endl;

// x0  = 0xFD;
// x1 = 0xBE;
// b0 = 1;
// b1 = 1;
// expectedResult = ((x0 ^ x1) & multiplicationConstants[b0 ^ b1]);
// result = DuAttalahMultiplication(x0, b0, x1, b1);
// cout << "Expected: " << expectedResult.to_ulong() << endl
//     << "Actual: " << (result[0] ^ result[1]).to_ulong() << endl;
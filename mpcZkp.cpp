#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h> 
#include "dpf.h"


using namespace std;

const block multiplicationConstants[2] = {0x0000, 0xFFFF};

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

block* mpcFirstStage(LowMC key, block CW, bool b0, bool b1) {
    // P_0 computes L_0 || R_0 = G(seed_0) + b_0 * CW
    block expansion0[2];
    // Call expansion(...) here
    expansion0[L] = 0b11;
    expansion0[L] ^= multiplicationConstants[b0] & CW;
    expansion0[R] = 0b01;
    expansion0[R] ^= multiplicationConstants[b0] & CW;

    // P_1 computes L_1 || R_1 = G(seed_1) + b_1 * CW
    block expansion1[2];
    // Call expansion(...) here
    expansion1[L] = 0b11;
    expansion1[L] ^= multiplicationConstants[b1] & CW;
    expansion1[R] = 0b10;
    expansion1[R] ^= multiplicationConstants[b1] & CW;

    // Use Du-Attalah multiplication to compute shares (1 - b) * L + b * R and to compute shares of b * L + (1 - b) * R
    block* bLShares = DuAttalahMultiplication(expansion0[L], b0, expansion1[L], b1); // L * b
    block* bRShares = DuAttalahMultiplication(expansion0[R], b0, expansion1[R], b1); // R * b
    // Compute (1 - b) * L as L - b * L
    block bNotLShares[2];
    bNotLShares[0] = expansion0[L] ^ bLShares[0];
    bNotLShares[1] = expansion1[L] ^ bLShares[1];
    // Compute (1 - b) * R as R - b * R
    block bNotRShares[2];
    bNotRShares[0] = expansion0[R] ^ bRShares[0];
    bNotRShares[1] = expansion1[R] ^ bRShares[1];

    block correctedLeftSideShares[2];
    correctedLeftSideShares[0] = bNotLShares[0] ^ bRShares[0];
    correctedLeftSideShares[1] = bNotLShares[1] ^ bRShares[1];

    block* correctedRightSideShares = new block[2];
    correctedRightSideShares[0] = bNotRShares[0] ^ bLShares[0];
    correctedRightSideShares[1] = bNotRShares[1] ^ bLShares[1];

    // P0 and P1 reveal their shares of the corrected left side and calculate the result.
    sendData(0, 1, correctedLeftSideShares[0], "Corrected left side");
    sendData(1, 0, correctedLeftSideShares[1], "Corrected left side");

    block correctedLeftSide = correctedLeftSideShares[0] ^ correctedLeftSideShares[1];
    printf("Expected: 0\tActual: %ld\n", correctedLeftSide.to_ulong());

    // Cleanup
    delete[] bRShares;
    delete[] bLShares;
    
    return correctedRightSideShares; // Output the corrected right side shares to be used as shares of the seeds for the next layer.
}


int main() {
    typedef __m256i __mX;
    LowMC key(1);

    // The prover must generate the DPF keys and then the zero-knowledge proof of its correctness.
    const size_t nitems   =  1ULL << 15; // TODO: Change this value to reflect the actual number of items.
    dpf_key<__mX, nitems> dpfkey[2] = { 0 }; 
    gen(key, 2 , dpfkey);

    // Apply layer one of the proof generation MPC protocol
    block* result = mpcFirstStage(key, 0b11, 1, 0);
    cout << "Results: " << (result[0] ^ result[1]).to_ulong() << endl << endl;
    
    cout << "Party 0:" << endl;
    printTranscript(0);
    cout << endl;
    cout << "Party 1:" << endl;
    printTranscript(1);
    cout << endl;
    cout << "Party 2:" << endl;
    printTranscript(2);

    delete[] result;
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
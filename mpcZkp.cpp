#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>


using namespace std;

struct TranscriptEntry {
    bool generated; // True for data which is being randomly generated
    int sender; // The sender of the data (0 for generated data)
    int receiver; // The receiver of the data
    int data; // The data
};

// The transcript of the MPC.
TranscriptEntry transcripts[3][128]; /*Determine the correct number for this*/

TranscriptEntry** mpcFirstStage(int CW, bool b0, bool b1) {
    srand(time(NULL)); // Replace rand and srand with a secure PRG.

    int p0Index = 0;
    int p1Index = 0;
    int p2Index = 0;

    // P_0 computes L_0 || R_0 = G(seed_0) + b_0 * CW
    int tmp0 = /*PRG goes here +*/ b0 * CW;
    int L0 = (tmp0 >> (32/2)) & 0xFFFF;
    int R0 = tmp0  & 0xFFFF;

    // P_1 computes L_1 || R_1 = G(seed_1) + b_1 * CW
    int tmp1 = /*PRG goes here +*/ b1 * CW;
    int L1 = (tmp1 >> (32/2)) & 0xFFFF;
    int R1 = tmp1  & 0xFFFF;

    // (1 - b)
    int bNot0 = ~b0;
    int bNot1 = b1;

    // Use Du-Attalah multiplication to compute shares (1 - b) * L + b * R
    // P0 starts with L0 and bNot0
    // P1 starts with L1 and bNot1

    // P2 generates random values D0, d0, D1, d1 and alpha
    int D0 = rand(); // TODO: Replace rand
    transcripts[2][p2Index].generated = true;
    transcripts[2][p2Index].receiver = 2;
    transcripts[2][p2Index].data = D0;
    p2Index++;

    int D1 = rand(); // TODO: Replace rand
    transcripts[2][p2Index].generated = true;
    transcripts[2][p2Index].receiver = 2;
    transcripts[2][p2Index].data = D1;
    p2Index++;

    int d0 = rand(); // TODO: Replace rand
    transcripts[2][p2Index].generated = true;
    transcripts[2][p2Index].receiver = 2;
    transcripts[2][p2Index].data = d0;
    p2Index++;

    int d1 = rand(); // TODO: Replace rand
    transcripts[2][p2Index].generated = true;
    transcripts[2][p2Index].receiver = 2;
    transcripts[2][p2Index].data = d1;
    p2Index++;

    int alpha = rand();  // TODO: Replace rand
    transcripts[2][p2Index].generated = true;
    transcripts[2][p2Index].receiver = 2;
    transcripts[2][p2Index].data = alpha;
    p2Index++;

    int c0 = D0 * d1 + alpha;
    int c1 = D1*d0 - alpha;

    // P2 sends D0, d0, c0 to P0
    TranscriptEntry D0FromP2ToP0;
    D0FromP2ToP0.generated = false;
    D0FromP2ToP0.sender = 2;
    D0FromP2ToP0.receiver = 0;
    D0FromP2ToP0.data = D0;
    transcripts[0][p0Index] = D0FromP2ToP0;
    p0Index++;
    transcripts[2][p2Index] = D0FromP2ToP0;
    p2Index++;

    TranscriptEntry d0FromP2ToP0;
    d0FromP2ToP0.generated = false;
    d0FromP2ToP0.sender = 2;
    d0FromP2ToP0.receiver = 0;
    d0FromP2ToP0.data = d0;
    transcripts[0][p0Index] = d0FromP2ToP0;
    p0Index++;
    transcripts[2][p2Index] = d0FromP2ToP0;
    p2Index++;

    TranscriptEntry c0FromP2ToP0;
    c0FromP2ToP0.generated = false;
    c0FromP2ToP0.sender = 2;
    c0FromP2ToP0.receiver = 0;
    c0FromP2ToP0.data = c0;
    transcripts[0][p0Index] = c0FromP2ToP0;
    p0Index++;
    transcripts[2][p2Index] = c0FromP2ToP0;
    p2Index++;

    // P2 sends D1, d1, and c1 to P1
    TranscriptEntry D1FromP2ToP1;
    D1FromP2ToP1.generated = false;
    D1FromP2ToP1.sender = 2;
    D1FromP2ToP1.receiver = 1;
    D1FromP2ToP1.data = D1;
    transcripts[1][p0Index] = D1FromP2ToP1;
    p1Index++;
    transcripts[2][p2Index] = D1FromP2ToP1;
    p2Index++;

    TranscriptEntry d1FromP2ToP1;
    d1FromP2ToP1.generated = false;
    d1FromP2ToP1.sender = 2;
    d1FromP2ToP1.receiver = 1;
    d1FromP2ToP1.data = d1;
    transcripts[1][p0Index] = d1FromP2ToP1;
    p1Index++;
    transcripts[2][p2Index] = d1FromP2ToP1;
    p2Index++;

    TranscriptEntry c1FromP2ToP1;
    c1FromP2ToP1.generated = false;
    c1FromP2ToP1.sender = 2;
    c1FromP2ToP1.receiver = 1;
    c1FromP2ToP1.data = c1;
    transcripts[1][p0Index] = c1FromP2ToP1;
    p1Index++;
    transcripts[2][p2Index] = c1FromP2ToP1;
    p2Index++;

    // P0 Computes X0 = x0 + D0 and Y0 = b0 + d0 and sends them to P1
    int X0 = L0 + D0;
    int Y0 = bNot0 + d0;
    TranscriptEntry X0FromP0ToP1;
    X0FromP0ToP1.generated = false;
    X0FromP0ToP1.sender = 0;
    X0FromP0ToP1.receiver = 1;
    X0FromP0ToP1.data = X0;
    transcripts[0][p0Index] = X0FromP0ToP1;
    p0Index++;
    transcripts[1][p1Index] = X0FromP0ToP1;
    p1Index++;
    TranscriptEntry Y0FromP0ToP1;
    Y0FromP0ToP1.generated = false;
    Y0FromP0ToP1.sender = 0;
    Y0FromP0ToP1.receiver = 1;
    Y0FromP0ToP1.data = Y0;
    transcripts[0][p0Index] = Y0FromP0ToP1;
    p0Index++;
    transcripts[1][p1Index] = Y0FromP0ToP1;
    p1Index++;

    // P1 computes X1 = x1 + D1 and Y1 = b1 + d1 and sends them to P0
    int X1 = L1 + D1;
    int Y1 = bNot1 + d1;
    TranscriptEntry X1FromP1ToP0;
    X1FromP1ToP0.generated = false;
    X1FromP1ToP0.sender = 1;
    X1FromP1ToP0.receiver = 0;
    X1FromP1ToP0.data = X1;
    transcripts[1][p1Index] = X1FromP1ToP0;
    p1Index++;
    transcripts[0][p0Index] = X1FromP1ToP0;
    p0Index++;
    TranscriptEntry Y1FromP1ToP0;
    Y1FromP1ToP0.generated = false;
    Y1FromP1ToP0.sender = 1;
    Y1FromP1ToP0.receiver = 0;
    Y1FromP1ToP0.data = Y1;
    transcripts[1][p1Index] = Y1FromP1ToP0;
    p1Index++;
    transcripts[0][p0Index] = Y1FromP1ToP0;
    p0Index++;

    int z0 = L0 * (bNot0 + Y1) - d0 * X1 + c0;
    int z1 = L1 * (bNot1 + Y0) - d1 * X0 + c1;

    return (TranscriptEntry**) transcripts;
}
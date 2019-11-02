#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>


using namespace std;

struct TranscriptEntry {
    bool generated; // True for data which is being randomly generated
    int sender; // The sender of the data (same as receiver for generated data)
    int receiver; // The receiver of the data
    int data; // The data
    string label;
};

// The transcript of the MPC.
// TranscriptEntry transcripts[3][128]; /*Determine the correct number for this*/
vector<TranscriptEntry> transcripts[3];

void generateData(int party, int data, string label) {
	TranscriptEntry t;
    t.generated = true;
    t.sender = party;
    t.receiver = party;
    t.data = data;
    t.label = label;

    transcripts[party].push_back(t);
}

void sendData(int sender, int receiver, int data, string label) {
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
int* DuAttalahMultiplication(int x0, bool b0, int x1, bool b1)  {
    int* results = new int[2];

    // P2 generates random values D0, d0, D1, d1 and alpha
    int D0 = rand() && 0xFF; // TODO: Replace rand
    generateData(2, D0, "D_0");

    int D1 = rand() && 0xFF; // TODO: Replace rand
    generateData(2, D1, "D1");

    int d0 = rand() && 0xFF; // TODO: Replace rand
    generateData(2, d0, "d_0");

    int d1 = rand() && 0xFF; // TODO: Replace rand
    generateData(2, d1, "d_1");

    int alpha = rand() && 0xFF;  // TODO: Replace rand
    generateData(2, alpha, "alpha");

    // P2 calculates c0 and c1.
    int c0 = (D0 * d1) ^ alpha;
    int c1 = (D1*d0) ^ alpha;

    // P2 sends D0, d0, c0 to P0
    sendData(2, 0, D0, "D_0");
    sendData(2, 0, d0, "d_0");
    sendData(2, 0, c0, "c_0 = D_0 * d_1 + alpha");
    // P2 sends D1, d1, and c1 to P1
    sendData(2, 1, D1, "D_1");
    sendData(2, 1, d1, "d_1");
    sendData(2, 1, c1, "c_1 = D_1 * d_0 - alpha");

    // P0 Computes X0 = x0 + D0 and Y0 = b0 + d0 and sends them to P1
    int X0 = x0 ^ D0;
    int Y0 = b0 ^ d0;

    sendData(0, 1, X0, "X_0 = x_0 + D_0");
    sendData(0, 1, Y0, "Y_0 = b_0 + d_0");

    // P1 computes X1 = x1 + D1 and Y1 = b1 + d1 and sends them to P0
    int X1 = x1 ^ D1;
    int Y1 = b1 ^ d1;

    sendData(1, 0, X1, "X_1 = x_1 + D_1");
    sendData(1, 0, Y1, "Y_1 = b_1 + d_1");

    results[0] = (x0 * (b0 ^ Y1)) ^ (d0 * X1) ^ c0;
    results[1] = (x1 * (b1 ^ Y0)) ^ (d1 * X0) ^ c1;
    return results;
}

// vector<TranscriptEntry>* mpcFirstStage(int CW, bool b0, bool b1) {
//     srand(time(NULL)); // Replace rand and srand with a secure PRG.

//     int p0Index = 0;
//     int p1Index = 0;
//     int p2Index = 0;

//     // P_0 computes L_0 || R_0 = G(seed_0) + b_0 * CW
//     int tmp0 = /*PRG goes here +*/ b0 * CW;
//     int L0 = (tmp0 >> (32/2)) & 0xFFFF;
//     int R0 = tmp0  & 0xFFFF;

//     // P_1 computes L_1 || R_1 = G(seed_1) + b_1 * CW
//     int tmp1 = /*PRG goes here +*/ b1 * CW;
//     int L1 = (tmp1 >> (32/2)) & 0xFFFF;
//     int R1 = tmp1  & 0xFFFF;

//     // (1 - b)
//     int bNot0 = ~b0;
//     int bNot1 = b1;

//     // Use Du-Attalah multiplication to compute shares (1 - b) * L + b * R
//     // P0 starts with L0 and bNot0
//     // P1 starts with L1 and bNot1

    

//     return transcripts;
// }

int main() {
    int x0 = 0xFD;
    int x1 = 0xBE;
    bool b0 = 1;
    bool b1 = 0;
    int *result = DuAttalahMultiplication(x0, b0, x1, b1);
    cout << "Expected: " << ((x0 ^ x1) *(b0 ^ b1)) << endl
        << "Actual: " << (result[0] ^ result[1]) << endl;
    
    delete[] result;
}
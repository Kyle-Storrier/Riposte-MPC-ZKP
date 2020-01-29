#include<array>
#include <x86intrin.h>  // SSE and AVX intrinsics

std::array<std::array<uint64_t, 256>, 4> demoMatrix1;
 std::array<__m128i, 256> demoMatrix2;
std::array<__m128i, 256> resultMatrix;

inline void generateLUT(const int i, const std::array<__m128i, 256> & matrix, __m128i lut[16][16]) {
    for(int k = 0; k < 16; k++) {
        lut[k][0b0000] = _mm_setzero_si128();
        lut[k][0b0001] = matrix[64*i+4*k+0];
        lut[k][0b0010] = matrix[64*i+4*k+1];
        lut[k][0b0011] = _mm_xor_si128(lut[k][0b0010], lut[k][0b0001]);
        lut[k][0b0100] = matrix[64*i+4*k+2];
        lut[k][0b0101] = _mm_xor_si128(lut[k][0b0100], lut[k][0b0001]);
        lut[k][0b0110] = _mm_xor_si128(lut[k][0b0100], lut[k][0b0010]);
        lut[k][0b0111] = _mm_xor_si128(lut[k][0b0100], lut[k][0b0011]);
        lut[k][0b1000] = matrix[64*i+4*k+3];
        lut[k][0b1001] = _mm_xor_si128(lut[k][0b1000], lut[k][0b0001]);
        lut[k][0b1010] = _mm_xor_si128(lut[k][0b1000], lut[k][0b0010]);
        lut[k][0b1011] = _mm_xor_si128(lut[k][0b1000], lut[k][0b0011]);
        lut[k][0b1100] = _mm_xor_si128(lut[k][0b1000], lut[k][0b0100]);
        lut[k][0b1101] = _mm_xor_si128(lut[k][0b1000], lut[k][0b0101]);
        lut[k][0b1110] = _mm_xor_si128(lut[k][0b1000], lut[k][0b0110]);
        lut[k][0b1111] = _mm_xor_si128(lut[k][0b1000], lut[k][0b0111]);
    }
}

void matrixMultiply(
        const std::array<std::array<uint64_t, 256>, 4> matrix1,
        std::array<__m128i, 256> & matrix2,
        std::array<__m128i, 256> & answer) {
    __m128i lut[16][16];
    for(int i = 0; i < 4; i++) {
        generateLUT(i, matrix2, lut);
        for(int j = 0; j < 256; j++) {
            uint64_t tmp = matrix1[i][j];
            for(int k = 0; k < 16; k++, tmp >>= 4) {
                answer[j] = _mm_xor_si128(answer[j], lut[k][tmp & 0xf]);
            }
        }
    }
}

int main() {
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 256; j++) {
            demoMatrix1[i][j] = 0;
        }
    }
    demoMatrix1[0][0] = 1;

    for(int i = 0; i < 128; i++) {
        resultMatrix[i] = _mm_setzero_si128();
    }

    for(int j = 0; j < 256; j++) {
        demoMatrix2[j] = _mm_setzero_si128();
    }
    demoMatrix2[0] = _mm_set_epi32(0, 0, 0, 1 << 1);

    matrixMultiply(demoMatrix1, demoMatrix2, resultMatrix);

    printf("Matrix 1\n");
    for(int i = 0; i < 256; i++) {
        for(int j = 0; j < 4; j++) {
            for(int k = 0; k < 64; k++) {
                printf("%ld ", (demoMatrix1[j][i] >> k) & 0b1);
            }
        }
        printf("\n");
    }

    printf("\n\nMatrix 2:\n");

    for(int j = 0; j < 256; j++) {
        uint64_t data[2];
        data[0] = _mm_extract_epi64(demoMatrix2[j], 0);
        data[1] = _mm_extract_epi64(demoMatrix2[j], 1);
        for(int i = 0; i < 2; i++) {
            for(int k = 0; k < 64; k++) {
                printf("%ld ", (data[i] >> k) & 0b1);
            }
        }
        printf("\n");
    }

    printf("\n\nResults Matrix:\n");

    for(int j = 0; j < 256; j++) {
        uint64_t data[2];
        data[0] = _mm_extract_epi64(resultMatrix[j], 0);
        data[1] = _mm_extract_epi64(resultMatrix[j], 1);
        for(int i = 0; i < 2; i++) {
            for(int k = 0; k < 64; k++) {
                printf("%ld ", (data[i] >> k) & 0b1);
            }
        }
        printf("\n");
    }
}
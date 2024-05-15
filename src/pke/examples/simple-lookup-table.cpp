//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  Examples for using CKKS scheme to perform simple lookup table operations.
 */

#include "openfhe.h"

#include <complex>
#include <chrono>
#include <cassert>

using namespace lbcrypto;

// gets a vector of integer (0 to codedim-1), and codedim, and return exp(-2 pi i * x / codedim)
void GetExponent(std::vector<uint32_t> input, size_t logDim, std::vector<std::complex<double>> &result);
// inverse of GetExponent
std::vector<int32_t> GetLog(std::vector<std::complex<double>> intput, size_t logDim);
void FindLUTPoly(std::vector<uint32_t> table, size_t logDim, 
                 std::vector<std::complex<double>> &coeff); 
void EvalLUT(Ciphertext<DCRTPoly> &ct, size_t logDim, std::vector<std::complex<double>> coeffs, 
             std::vector<Ciphertext<DCRTPoly>> &dest);
void FindLUTPoly2D(std::vector<uint32_t> table, size_t logDim, 
                   std::vector<std::vector<std::complex<double>>> &coeff);
void EvalLUTs2D(Ciphertext<DCRTPoly> &ct0, Ciphertext<DCRTPoly> &ct1, 
                size_t logDim, std::vector<std::vector<std::vector<std::complex<double>>>> coeffs,
                std::vector<Ciphertext<DCRTPoly>> &dest);
void EvalLUTs2DAsymmetric(Ciphertext<DCRTPoly> &ct0, Ciphertext<DCRTPoly> &ct1, 
                size_t logDim0, size_t logDim1, std::vector<std::vector<std::vector<Plaintext>>> coeffs, 
                std::vector<Ciphertext<DCRTPoly>> &dest);
void EvalNoiseReductionInplace(Ciphertext<DCRTPoly> ct, size_t logDim, CryptoContext<DCRTPoly> cc);
void ifft(std::vector<std::complex<double>>& data);
void ifft2(std::vector<std::vector<std::complex<double>>>& data);

void testLUTNtoN();
void testLUTNtoBN();
void testLUTANtoBN();
void testLUTANtoBNAsymmetric();

void testFFT();
void testExponent();
void testComplexOperation();
void testLog();

int main() {
    // testLUTNtoN();
    // testLUTNtoBN();
    // testLUTANtoBN();
    testLUTANtoBNAsymmetric();
    // testLog();

    return 0;
}

// TODO: use ifft defined in src/core/include/math/dftransform.h
// Bit-reverse copy function
void BitReverse(std::vector<std::complex<double>>& a) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
}

// Iterative FFT function
void fft(std::vector<std::complex<double>>& a, bool invert) {
    int n = a.size();
    BitReverse(a);

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * M_PI / len * (invert ? 1 : -1);
        std::complex<double> wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (int j = 0; j < len / 2; j++) {
                std::complex<double> u = a[i + j];
                std::complex<double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (std::complex<double>& x : a) {
            x /= n;
        }
    }
}

// Wrapper function for IFFT
void ifft(std::vector<std::complex<double>>& a) {
    fft(a, true);
}

void ifft2(std::vector<std::vector<std::complex<double>>>& data) {
    // Apply IFFT on each row
    for (auto& row : data) {
        ifft(row);
    }

    // Apply IFFT on each column
    int rows = data.size();
    int cols = data[0].size();
    for (int col = 0; col < cols; ++col) {
        std::vector<std::complex<double>> columnData(rows);
        for (int row = 0; row < rows; ++row) {
            columnData[row] = data[row][col];
        }

        // Perform IFFT on this column
        ifft(columnData);

        // Place the result back into the data array
        for (int row = 0; row < rows; ++row) {
            data[row][col] = columnData[row];
        }
    }
}

// gets a vector of integer (0 to codedim-1), and codedim, and return exp(-2 pi i * x / codedim)
void GetExponent(std::vector<uint32_t> input, size_t logDim,
                 std::vector<std::complex<double>> &result){
    size_t codedim = 1 << logDim;
    for (size_t i = 0; i < input.size(); i++){
        result.push_back(std::exp(std::complex<double>(0, -2 * M_PI * input[i] / codedim)));
    }
}

std::vector<int32_t> GetLog(std::vector<std::complex<double>> intput, size_t logDim){
    size_t codedim = 1 << logDim;
    std::vector<int32_t> result;
    std::complex<double> scale = std::complex<double>(0, -2 * M_PI / codedim);
    for (size_t i = 0; i < intput.size(); i++){
        std::complex<double> res = std::log(intput[i]) / scale;
        int resInt = round(res.real());
        result.push_back((resInt + codedim) % codedim);
    }
    return result;
}

void FindLUTPoly(std::vector<uint32_t> table, size_t logDim, std::vector<std::complex<double>> &coeff){
    GetExponent(table, logDim, coeff);
    // Perform inverse FFT on coeff
    ifft(coeff);
}

void FindLUTPoly2D(std::vector<uint32_t> table, size_t logDim, 
                   std::vector<std::vector<std::complex<double>>> &coeff){
    size_t codedim = 1 << logDim;
    for (size_t i = 0; i < codedim; i++){
        std::vector<std::complex<double>> row;
        std::vector<uint32_t> chunk(table.begin() + i*codedim, table.begin() + (i+1)*codedim);
        GetExponent(chunk, logDim, row);
        coeff.push_back(row);
    }
    // Perform inverse FFT on coeff
    ifft2(coeff);
}

void FindLUTPoly2DAsymmetric(std::vector<uint32_t> table, size_t logDim0, size_t logDim1, size_t logDimOut, 
                   std::vector<std::vector<std::complex<double>>> &coeff){
    size_t codedim0 = 1 << logDim0;
    size_t codedim1 = 1 << logDim1;

    for (size_t i = 0; i < codedim1; i++){
        std::vector<std::complex<double>> row;
        std::vector<uint32_t> chunk(table.begin() + i*codedim0, table.begin() + (i+1)*codedim0);
        GetExponent(chunk, logDimOut, row);
        coeff.push_back(row);
    }
    // Perform inverse FFT on coeff
    ifft2(coeff);
}
// evaluate n-to-n LUT
void EvalLUT(Ciphertext<DCRTPoly> &ct, size_t logDim, 
             std::vector<Plaintext> coeffs, 
             Ciphertext<DCRTPoly> &dest){

    auto cc = ct->GetCryptoContext();
    auto algo = cc->GetScheme();
    
    size_t codedim = 1 << logDim;
    std::vector<Ciphertext<DCRTPoly>> b(codedim - 1);

    // find power basis, b[i] is for T_{i+1}
    b[0] = ct;
    for (size_t i = 2; i < codedim; i++){
        auto i0 = i/2;
        auto i1 = i - i0;
        b[i - 1] = cc->EvalMult(b[i0 - 1], b[i1 - 1]);
    }

    // match the level and depth
    for (size_t i = 0; i < codedim/2; i++){ 
        algo->AdjustLevelsAndDepthInPlace(b[i], b[codedim-2]);
    }

    // inner product with coefficients
    dest = cc->EvalMult(b[0], coeffs[1]);
    for (size_t i = 2; i < codedim; i++){
        auto tmp = cc->EvalMult(b[i - 1], coeffs[i]);
        dest = cc->EvalAdd(dest, tmp);
    }
    dest = cc->EvalAdd(dest, coeffs[0]);
}

// evaluate n-to-b*n LUT
void EvalLUTs(Ciphertext<DCRTPoly> &ct, size_t logDim, 
             std::vector<std::vector<Plaintext>> coeffs, 
             std::vector<Ciphertext<DCRTPoly>> &dest){

    auto cc = ct->GetCryptoContext();
    auto algo = cc->GetScheme();

    size_t codedim = 1 << logDim;
    std::vector<Ciphertext<DCRTPoly>> b(codedim - 1);

    // find power basis, b[i] is for T_{i+1}
    b[0] = ct;
    for (size_t i = 2; i < codedim; i++){
        auto i0 = i/2;
        auto i1 = i - i0;
        b[i - 1] = cc->EvalMult(b[i0 - 1], b[i1 - 1]);
    }

    // match the level and depth
    for (size_t i = 0; i < codedim/2; i++){ 
        algo->AdjustLevelsAndDepthInPlace(b[i], b[codedim-2]);
    }

    for(auto &table : coeffs){
        // inner product with coefficients
        dest.push_back(cc->EvalMult(b[0], table[1]));
        for (size_t i = 2; i < codedim; i++){
            auto tmp = cc->EvalMult(b[i - 1], table[i]);
            dest.back() = cc->EvalAdd(dest.back(), tmp);
        }
        dest.back() = cc->EvalAdd(dest.back(), table[0]);
    }
}

// evaluate 2n-to-b*n LUT
void EvalLUTs2D(Ciphertext<DCRTPoly> &ct0, Ciphertext<DCRTPoly> &ct1, 
                size_t logDim, std::vector<std::vector<std::vector<Plaintext>>> coeffs, 
                std::vector<Ciphertext<DCRTPoly>> &dest){
    auto cc = ct0->GetCryptoContext();
    auto algo = cc->GetScheme();

    size_t codedim = 1 << logDim;

    std::vector<Ciphertext<DCRTPoly>> b0(codedim - 1);
    std::vector<Ciphertext<DCRTPoly>> b1(codedim - 1);

    // find power basis
    b0[0] = ct0, b1[0] = ct1;
    for (size_t i = 2; i <= codedim/2; i++){
        auto i0 = i/2;
        auto i1 = i - i0;
        b0[i - 1] = cc->EvalMult(b0[i0 - 1], b0[i1 - 1]);
        b1[i - 1] = cc->EvalMult(b1[i0 - 1], b1[i1 - 1]);
    }

    // match the level and depth
    for (size_t i = 0; i < codedim/2 - 1; i++){ 
        algo->AdjustLevelsAndDepthInPlace(b0[i], b0[codedim/2 - 1]);
        algo->AdjustLevelsAndDepthInPlace(b1[i], b1[codedim/2 - 1]);
    }

    for (size_t i = codedim/2 + 1; i < codedim; i++){
        auto i0 = codedim - i;

        b0[i - 1] = cc->EvalConjugate(b0[i0 - 1]);
        b1[i - 1] = cc->EvalConjugate(b1[i0 - 1]);
    }

    // evaluate the polynomial
    for (auto &table : coeffs){
        // inner product with coefficients
        Ciphertext<DCRTPoly> tableEval;
        for (size_t i = 0; i < codedim; i++){
            // Evaluate i-th row
            auto rowEval = cc->EvalMult(b0[0], table[i][1]);
            for (size_t j = 2; j < codedim; j++){
                auto tmp = cc->EvalMult(b0[j - 1], table[i][j]);
                rowEval = cc->EvalAdd(rowEval, tmp);
            }
            rowEval = cc->EvalAdd(rowEval, table[i][0]);

            // Add to the result
            if (i == 0){
                tableEval = rowEval;
            } else {
                auto tmp = cc->EvalMult(rowEval, b1[i - 1]);
                tableEval = cc->EvalAdd(tableEval,tmp);
            }
        }
        dest.push_back(tableEval);
    }
}

// evaluate 2n-to-b*n LUT but two values has different bit length
void EvalLUTs2DAsymmetric(Ciphertext<DCRTPoly> &ct0, Ciphertext<DCRTPoly> &ct1, 
                size_t logDim0, size_t logDim1, std::vector<std::vector<std::vector<Plaintext>>> coeffs, 
                std::vector<Ciphertext<DCRTPoly>> &dest){
    auto cc = ct0->GetCryptoContext();
    auto algo = cc->GetScheme();

    size_t codedim0 = 1 << logDim0;
    size_t codedim1 = 1 << logDim1;

    std::vector<Ciphertext<DCRTPoly>> b0(codedim0 - 1);
    std::vector<Ciphertext<DCRTPoly>> b1(codedim1 - 1);

    // find power basis
    b0[0] = ct0;
    for (size_t i = 2; i <= codedim0/2; i++){
        auto i0 = i/2;
        auto i1 = i - i0;
        b0[i - 1] = cc->EvalMult(b0[i0 - 1], b0[i1 - 1]);
    }

    b1[0] = ct1;
    for (size_t i = 2; i <= codedim1/2; i++){
        auto i0 = i/2;
        auto i1 = i - i0;
        b1[i - 1] = cc->EvalMult(b1[i0 - 1], b1[i1 - 1]);
    }

    // match the level and depth
    for (size_t i = 0; i < codedim0/2 - 1; i++){ 
        algo->AdjustLevelsAndDepthInPlace(b0[i], b0[codedim0/2 - 1]);
    }
    for (size_t i = 0; i < codedim1/2 - 1; i++){ 
        algo->AdjustLevelsAndDepthInPlace(b1[i], b1[codedim1/2 - 1]);
    }

    for (size_t i = codedim0/2 + 1; i < codedim0; i++){
        auto i0 = codedim0 - i;
        b0[i - 1] = cc->EvalConjugate(b0[i0 - 1]);
    }
    for (size_t i = codedim1/2 + 1; i < codedim1; i++){
        auto i0 = codedim1 - i;
        b1[i - 1] = cc->EvalConjugate(b1[i0 - 1]);
    }

    // evaluate the polynomial
    for (auto &table : coeffs){
        // innerproduct with coefficients
        Ciphertext<DCRTPoly> tableEval;
        for (size_t i = 0; i < codedim1; i++){
            // Evaluate i-th row
            auto rowEval = cc->EvalMult(b0[0], table[i][1]);
            for (size_t j = 2; j < codedim0; j++){
                auto tmp = cc->EvalMult(b0[j - 1], table[i][j]);
                rowEval = cc->EvalAdd(rowEval, tmp);
            }
            rowEval = cc->EvalAdd(rowEval, table[i][0]);

            // Add to the result
            if (i == 0){
                tableEval = rowEval;
            } else {
                auto tmp = cc->EvalMult(rowEval, b1[i - 1]);
                tableEval = cc->EvalAdd(tableEval,tmp);
            }
        }
        dest.push_back(tableEval);
    }
}

// Evaluate -1/n x^{n+1} + (1 + 1/n) x
void EvalNoiseReductionInplace(CryptoContext<DCRTPoly> &cc
    , Ciphertext<DCRTPoly> &ct, size_t logDim){
    size_t codedim = 1 << logDim;

    auto low_term = cc->EvalMult(ct, 1. + 1./codedim);
    auto high_term = cc->EvalMult(ct, -1./codedim);
    
    for (size_t i = 1; i < codedim; i *= 2){
        ct = cc->EvalMult(ct, ct);
    }

    ct = cc->EvalMult(ct, high_term);
    ct = cc->EvalAdd(ct, low_term);
}

void testFFT(){
    int rows = 4;
    int cols = 4;
    std::vector<std::vector<std::complex<double>>> a(rows, std::vector<std::complex<double>>(cols));

    int value = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            a[i][j] = value++;
        }
    }

    // Print the 2D vector to verify
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << a[i][j] << " ";
        }
        std::cout << std::endl;
    }
    ifft2(a);
    // Print the 2D vector to verify
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << a[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void testLUTNtoN(){
    size_t logDim = 4;
    std::cout << "Test single LUT" << std::endl;
    // Note: put any logDim-to-logDim table you want to test here 
    std::vector<uint32_t> table = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::cout << "Input table: \n====================\n";
    for (size_t i = 0; i < table.size(); i++){
        std::cout << i << "\t" << table[i] << std::endl;
    }

    std::vector<std::complex<double>> coeffs;
    FindLUTPoly(table, logDim, coeffs);

    std::vector<uint32_t> mInt = {1,3,5,7,0,2,4,6,11,12,13,14};
    std::vector<std::complex<double>> mComplex;
    GetExponent(mInt, logDim, mComplex);
    
    // Create the crypto context
    uint32_t multDepth = 5;
    uint32_t scaleModSize = 50;
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    // Encode coefficients
    std::vector<Plaintext> coeffsPt(coeffs.size());
    for (size_t i = 0; i < coeffs.size(); i++){
        coeffsPt[i] = cc->MakeCKKSPackedPlaintext(std::vector<std::complex<double>>(cc->GetRingDimension()/2, coeffs[i]));
    }

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(mComplex);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keys.publicKey, ptxt1);

    // Evaluate the LUT and measure the time
    auto start = std::chrono::high_resolution_clock::now();

    EvalLUT(ctxt1, logDim, coeffsPt, ctxt1);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Time taken: " << duration.count() << "ms" << std::endl;

    Plaintext result;
    cc->Decrypt(keys.secretKey, ctxt1, &result);

    std::vector<std::complex<double>> complexResult = result->GetCKKSPackedValue();
    auto intResult = GetLog(complexResult, logDim);

    std::cout << "Output table: \n====================\n";
    std::cout << "Input\tOutput\n";
    for (size_t i = 0; i < mInt.size(); i++){
        std::cout << mInt[i] << "\t" << intResult[i] << std::endl;
    }
}

void testLUTNtoBN(){
    size_t logDim = 4;
    size_t outB = 2;
    std::cout << "Test n-to-bn LUT" << std::endl;
    // Note: put any logDim-to-2*logDim table you want to test here 
    // The output range is dim^2, e.g, logDim = 4, the output range is 16^2
    std::vector<uint32_t> table(1<<logDim);
    for (size_t i = 0; i < table.size(); i++){
        table[i] = (int)(128*sin(i))+128; // Any cool function 
    }
    std::cout << "Input table: \n====================\n";
    for (size_t i = 0; i < table.size(); i++){
        std::cout << i << "\t" << table[i] << std::endl;
    }

    std::vector<std::vector<std::complex<double>>> coeffs(outB);
    // divide tables by chunks of size logDim-bit
    for (size_t i = 0; i < outB; i++){
        std::vector<uint32_t> partialTable(table.size());
        uint32_t mask = (1<<logDim) - 1;
        mask <<= logDim * i;
        for (size_t j = 0; j < table.size(); j++){
            partialTable[j] = (table[j] & mask) >> (logDim * i);
        }
        FindLUTPoly(partialTable, logDim, coeffs[i]);
    }

    std::vector<uint32_t> mInt = {1,3,5,7,9,11,13,15};
    std::vector<std::complex<double>> mComplex;
    GetExponent(mInt, logDim, mComplex);
    
    // Create the crypto context
    uint32_t multDepth = 5;
    uint32_t scaleModSize = 50;
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    // Encode coeffs
    std::vector<std::vector<Plaintext>> coeffsPt(outB);
    for (size_t i = 0; i < outB; i++){
        coeffsPt[i].resize(coeffs[i].size());
        for (size_t j = 0; j < coeffs[i].size(); j++){
            coeffsPt[i][j] = cc->MakeCKKSPackedPlaintext(std::vector<std::complex<double>>(cc->GetRingDimension()/2, coeffs[i][j]));
        }
    }

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(mComplex);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keys.publicKey, ptxt1);
    std::vector<Ciphertext<DCRTPoly>> ctxtout;
    // Evaluate the LUT and measure the time
    auto start = std::chrono::high_resolution_clock::now();

    EvalLUTs(ctxt1, logDim, coeffsPt, ctxtout);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Time taken: " << duration.count() << "ms" << std::endl;

    std::vector<Plaintext> result(outB);
    for (size_t i = 0; i < outB; i++){
        cc->Decrypt(keys.secretKey, ctxtout[i], &result[i]);
    }

    std::vector<std::vector<int>> intResult(outB);
    for (size_t i = 0; i < outB; i++){
        std::vector<std::complex<double>> complexResult = result[i]->GetCKKSPackedValue();
        intResult[i] = GetLog(complexResult, logDim);
    }

    if (outB != 2) {
        std::cout << "Please modify the output range to 2 to see output printed\n";
        return;
    }

    std::cout << "Output table: \n====================\n";
    std::cout << "Input\tOut(R)\tOut(L)\tMerged\n";
    for (size_t i = 0; i < mInt.size(); i++){
        std::cout << mInt[i] << "\t" << intResult[0][i] << "\t" << intResult[1][i] << "\t" << intResult[0][i]+ intResult[1][i]*(1<<logDim)  << std::endl;
    }
}

inline uint32_t presetCryptoContext(CryptoContext<DCRTPoly>& cc){
    // Create the crypto context
    uint32_t multDepth = 5;
    uint32_t scaleModSize = 50;
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);

    cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    return multDepth;
}

inline uint32_t presetCryptoContextBootstrap(CryptoContext<DCRTPoly>& cc){
    CCParams<CryptoContextCKKSRNS> parameters;

    // Create the crypto context
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    parameters.SetSecurityLevel(HEStd_128_classic);
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint dcrtBits               = 50;
    usint firstMod               = 55;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    std::vector<uint32_t> levelBudget = {3, 3};
    std::vector<uint32_t> bsgsDim = {0, 0};

    uint32_t levelsAvailableAfterBootstrap = 5;
    usint depth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);

    cc = GenCryptoContext(parameters);

    if(cc->GetRingDimension() > (1<<16)){
        std::cout << "Ring dimension is too large for this example: " << std::dec << cc->GetRingDimension() << std::endl;
        exit(-1);
    }

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);
    cc->EvalBootstrapSetup(levelBudget, bsgsDim, 1<<15);

    return depth;
}

void testLUTANtoBN(){
    size_t logDim = 4;
    // size_t inA = 2, // inA > 2 is not implemented yet
    size_t outB = 2; 
    std::cout << "Test an-to-bn LUT" << std::endl;

    std::cout << "Input invsbox: \n====================\n";
    // Note: put any 2logDim-to-2logDim table you want to test here 
    std::vector<uint32_t> invsbox = {0x52,0x09,0x6a,0xd5,0x30,0x36,0xa5,0x38,0xbf,0x40,0xa3,0x9e,0x81,0xf3,0xd7,0xfb,
                                    0x7c,0xe3,0x39,0x82,0x9b,0x2f,0xff,0x87,0x34,0x8e,0x43,0x44,0xc4,0xde,0xe9,0xcb,
                                    0x54,0x7b,0x94,0x32,0xa6,0xc2,0x23,0x3d,0xee,0x4c,0x95,0x0b,0x42,0xfa,0xc3,0x4e,
                                    0x08,0x2e,0xa1,0x66,0x28,0xd9,0x24,0xb2,0x76,0x5b,0xa2,0x49,0x6d,0x8b,0xd1,0x25,
                                    0x72,0xf8,0xf6,0x64,0x86,0x68,0x98,0x16,0xd4,0xa4,0x5c,0xcc,0x5d,0x65,0xb6,0x92,
                                    0x6c,0x70,0x48,0x50,0xfd,0xed,0xb9,0xda,0x5e,0x15,0x46,0x57,0xa7,0x8d,0x9d,0x84,
                                    0x90,0xd8,0xab,0x00,0x8c,0xbc,0xd3,0x0a,0xf7,0xe4,0x58,0x05,0xb8,0xb3,0x45,0x06,
                                    0xd0,0x2c,0x1e,0x8f,0xca,0x3f,0x0f,0x02,0xc1,0xaf,0xbd,0x03,0x01,0x13,0x8a,0x6b,
                                    0x3a,0x91,0x11,0x41,0x4f,0x67,0xdc,0xea,0x97,0xf2,0xcf,0xce,0xf0,0xb4,0xe6,0x73,
                                    0x96,0xac,0x74,0x22,0xe7,0xad,0x35,0x85,0xe2,0xf9,0x37,0xe8,0x1c,0x75,0xdf,0x6e,
                                    0x47,0xf1,0x1a,0x71,0x1d,0x29,0xc5,0x89,0x6f,0xb7,0x62,0x0e,0xaa,0x18,0xbe,0x1b,
                                    0xfc,0x56,0x3e,0x4b,0xc6,0xd2,0x79,0x20,0x9a,0xdb,0xc0,0xfe,0x78,0xcd,0x5a,0xf4,
                                    0x1f,0xdd,0xa8,0x33,0x88,0x07,0xc7,0x31,0xb1,0x12,0x10,0x59,0x27,0x80,0xec,0x5f,
                                    0x60,0x51,0x7f,0xa9,0x19,0xb5,0x4a,0x0d,0x2d,0xe5,0x7a,0x9f,0x93,0xc9,0x9c,0xef,
                                    0xa0,0xe0,0x3b,0x4d,0xae,0x2a,0xf5,0xb0,0xc8,0xeb,0xbb,0x3c,0x83,0x53,0x99,0x61,
                                    0x17,0x2b,0x04,0x7e,0xba,0x77,0xd6,0x26,0xe1,0x69,0x14,0x63,0x55,0x21,0x0c,0x7d};
    std::cout << std::hex;
    for (size_t i = 0; i < 16; i++){
        std::cout << i << "\t" << invsbox[i] << std::endl;
    }
    std::cout << "..." << std::endl << std::endl;

    std::vector<std::vector<std::vector<std::complex<double>>>> coeffs2D(2);
    // divide tables by chunks of size logDim-bit
    for (size_t i = 0; i < outB; i++){
        std::vector<uint32_t> partialTable(invsbox.size());
        uint32_t mask = (1<<logDim) - 1;
        mask <<= logDim * i;
        for (size_t j = 0; j < invsbox.size(); j++){
            partialTable[j] = (invsbox[j] & mask) >> (logDim * i);
        }
        FindLUTPoly2D(partialTable, logDim, coeffs2D[i]);
    }

    // Create the crypto context
    CryptoContext<DCRTPoly> cc;
    uint32_t maxDepth = presetCryptoContextBootstrap(cc);

    auto keys = cc->KeyGen();
    uint32_t numSlots = 1 << 15;

    std::cout << "Generating keys..." << std::endl;
    cc->EvalMultKeyGen(keys.secretKey);
    cc->EvalConjugateKeyGen(keys.secretKey);
    cc->EvalBootstrapKeyGen(keys.secretKey, numSlots);
    std::cout << "Keys generated." << std::endl;

    // Encode coeffs
    std::vector<std::vector<std::vector<Plaintext>>> coeffsPt(outB);
    for (size_t i = 0; i < outB; i++){
        coeffsPt[i].resize(coeffs2D[i].size());
        for (size_t j = 0; j < coeffs2D[i].size(); j++){
            coeffsPt[i][j].resize(coeffs2D[i][j].size());
            for (size_t k = 0; k < coeffs2D[i][j].size(); k++){
                coeffsPt[i][j][k] = cc->MakeCKKSPackedPlaintext(std::vector<std::complex<double>>(numSlots, coeffs2D[i][j][k]), 1, maxDepth - 2);
            }
        }
    }

    std::vector<uint32_t> mInt0(numSlots);
    std::vector<uint32_t> mInt1(numSlots);
    for (size_t i = 0; i < numSlots; i++){
        uint32_t input = i % 256; // input range is 0 to 255
        mInt0[i] = input & 0xF; // right four bits
        mInt1[i] = input >> 4; // left four bits
    }
    
    std::vector<std::complex<double>> mComplex0, mComplex1;
    GetExponent(mInt0, logDim, mComplex0);
    GetExponent(mInt1, logDim, mComplex1);

    uint32_t initLevel = 1;

    Plaintext ptxt0 = cc->MakeCKKSPackedPlaintext(mComplex0, 1, maxDepth - initLevel);
    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(mComplex1, 1, maxDepth - initLevel);
    Ciphertext<DCRTPoly> ctxt0 = cc->Encrypt(keys.publicKey, ptxt0);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keys.publicKey, ptxt1);
    auto b_start = std::chrono::high_resolution_clock::now();

    ctxt0 = cc->EvalBootstrap(ctxt0);
    ctxt1 = cc->EvalBootstrap(ctxt1);

    auto b_end = std::chrono::high_resolution_clock::now();
    auto b_duration = std::chrono::duration_cast<std::chrono::milliseconds>(b_end - b_start);
    std::cout << "Time taken (boot 2 times): " << std::dec << b_duration.count() << " ms" << std::endl;
    
    std::vector<Ciphertext<DCRTPoly>> ctxtout;

    // Evaluate the LUT and measure the time
    auto start = std::chrono::high_resolution_clock::now();

    EvalLUTs2D(ctxt0, ctxt1, logDim, coeffsPt, ctxtout);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time taken (4x2-to-4x2 LUT): " << std::dec << duration.count() << " ms" << std::endl;

    std::vector<Plaintext> result(outB);
    for (size_t i = 0; i < outB; i++){
        cc->Decrypt(keys.secretKey, ctxtout[i], &result[i]);
    }

    std::vector<std::vector<int>> intResult(outB);
    for (size_t i = 0; i < outB; i++){
        std::vector<std::complex<double>> complexResult = result[i]->GetCKKSPackedValue();
        intResult[i] = GetLog(complexResult, logDim);
    }

    std::cout << "Output table: \n====================\n";
    std::cout << "In(L)\tIn(R)\tMerged\t|Out(L)\tOut(R)\tMerged\tExpected\n";
    std::cout << std::hex;
    for (size_t i = 0; i < 8; i++){
        std::cout << mInt1[i] << "\t" << mInt0[i] << "\t" << mInt1[i]*(1<<logDim) + mInt0[i] << 
                "\t|" << intResult[1][i] << "\t" << intResult[0][i] << "\t" << intResult[1][i]*(1<<logDim) + intResult[0][i] <<
                "\t" << invsbox[mInt1[i]*(1<<logDim) + mInt0[i]]<< std::endl;
    }

    std::vector<std::vector<uint32_t>> expected(outB);
    for (size_t i = 0; i < outB; i++){
        expected[i].resize(numSlots);
        uint32_t mask = (1<<logDim) - 1;
        mask <<= logDim * i;
        for (size_t j = 0; j < numSlots; j++){
            expected[i][j] = (invsbox[j % (invsbox.size())] & mask) >> (logDim * i);
        }
    }

    std::cout << "Int and complex value table\n================================" << std::endl;
    std::cout << "Enced(int)\tEnced\tExpect(int)\tExpect" << std::endl;
    for (size_t i = 0; i < 8; i++){
        std::cout << intResult[0][i] << "\t" << result[0]->GetCKKSPackedValue()[i] << "\t";
        // calculate and print exp(-2*pi* i / 16 *expected[0][i])
        std::cout << expected[0][i] << "\t" << std::exp(std::complex<double>(0, -2 * M_PI * expected[0][i] / 16)) << std::endl;
    }

    uint32_t count_failure = 0;
    std::cout << std::dec;
    for (size_t i = 0; i < outB; i++){
        for (size_t j = 0; j < numSlots; j++){
            if(expected[i][j] != intResult[i][j]){
                std::cout << "Mismatch at " << i << " " << j << " " << expected[i][j] << " " 
                << intResult[i][j] << "(" << result[i]->GetCKKSPackedValue()[j] << ")" << std::endl;
                count_failure++;
            }
        }  
    }
    std::cout << "Number of failures: " << count_failure << "/" << numSlots*outB << std::endl;

    std::cout << "Noise variance is: ";
    double variance = 0.;
    for (size_t i = 0; i < numSlots; i++){
        // variance of noise, noise: result[0]->GetCKKSPackedValue()[i] - std::exp(std::complex<double>(0, -2 * M_PI * expected[0][i] / 16))
        variance += std::norm(result[0]->GetCKKSPackedValue()[i] - std::exp(std::complex<double>(0, -2 * M_PI * (int)expected[0][i] / 16)));
        variance += std::norm(result[1]->GetCKKSPackedValue()[i] - std::exp(std::complex<double>(0, -2 * M_PI * (int)expected[1][i] / 16)));
    }
    std::cout << variance / (2*numSlots) << std::endl;
}

void testLUTANtoBNAsymmetric(){
    size_t outB = 2; 

    size_t logDim[] = {4, 5}; // size of outB, input and output has same size.

    size_t inputLen = 0;
    for (auto &n : logDim){
       inputLen += n; 
    }
    std::cout << "inputLen: " << inputLen << "\n";
    
    std::cout << "Test 2n-to-bn LUT, Asymmetric with " << logDim[0] << " and " << logDim[1] << " bits inputs" << std::endl;

    std::cout << "Input table: \n====================\n";
    // Note: put any 2logDim-to-2logDim table you want to test here 
    std::vector<uint32_t> input_table(1<<inputLen);

    for (size_t i = 0; i < 1<<inputLen; i++){ // put any cool table you want
        // mod 5
        input_table[i] = i % 5;
    }

    std::cout << std::hex;
    for (size_t i = 0; i < 16; i++){
        std::cout << i << "\t" << input_table[i] << std::endl;
    }
    std::cout << "..." << std::endl << std::endl;

    std::vector<std::vector<std::vector<std::complex<double>>>> coeffs2D(2);
    // divide tables by chunks of size logDim-bit
    size_t nShift = 0; 
    for (size_t i = 0; i < outB; i++){
        // std::cout << "Coeffs(int) for output " << i << std::endl;
        std::vector<uint32_t> partialTable(input_table.size());
        uint32_t mask = (1<<logDim[i]) - 1;
        mask <<= nShift;
        for (size_t j = 0; j < input_table.size(); j++){
            partialTable[j] = (input_table[j] & mask) >> (nShift);
            // std::cout << partialTable[j] << "(" << invsbox[j] << ",  " << mask << ", " <<  (nShift)<< ") ";
        }
        // std::cout << std::endl;
        FindLUTPoly2DAsymmetric(partialTable, logDim[0], logDim[1], logDim[i], coeffs2D[i]);
        nShift += logDim[i];
    }

    // print coeffs2D
    // for (size_t i = 0; i < outB; i++){
    //     std::cout << "Coeffs for output " << i << std::endl;
    //     for (size_t j = 0; j < coeffs2D[i].size(); j++){
    //         for (size_t k = 0; k < coeffs2D[i][j].size(); k++){
    //             std::cout << coeffs2D[i][j][k] << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    // Create the crypto context
    CryptoContext<DCRTPoly> cc;
    uint32_t maxDepth = presetCryptoContextBootstrap(cc);

    auto keys = cc->KeyGen();
    uint32_t numSlots = 1 << 15;

    std::cout << "Generating keys..." << std::endl;
    cc->EvalMultKeyGen(keys.secretKey);
    cc->EvalConjugateKeyGen(keys.secretKey);
    // cc->EvalBootstrapKeyGen(keys.secretKey, numSlots);
    std::cout << "Keys generated." << std::endl;

    // Encode coeffs
    std::vector<std::vector<std::vector<Plaintext>>> coeffsPt(outB);
    for (size_t i = 0; i < outB; i++){
        coeffsPt[i].resize(coeffs2D[i].size());
        for (size_t j = 0; j < coeffs2D[i].size(); j++){
            coeffsPt[i][j].resize(coeffs2D[i][j].size());
            for (size_t k = 0; k < coeffs2D[i][j].size(); k++){
                coeffsPt[i][j][k] = cc->MakeCKKSPackedPlaintext(std::vector<std::complex<double>>(numSlots, coeffs2D[i][j][k]), 1, maxDepth - 2);
            }
        }
    }

    std::vector<std::vector<uint32_t>> mInt(2);
    nShift = 0;
    for (size_t i = 0; i < outB; i++){
        mInt[i].resize(numSlots);
        uint32_t mask = (1<<logDim[i]) - 1;
        mask <<= nShift;
        for (size_t j = 0; j < numSlots; j++){
            uint32_t input = j % (1 << inputLen); // input range is 2^inputLen
            mInt[i][j] = (input & mask) >> (nShift);
        }
        nShift += logDim[i];
    }
    
    std::vector<std::vector<std::complex<double>>> mComplex(2);
    
    GetExponent(mInt[0], logDim[0], mComplex[0]);
    GetExponent(mInt[1], logDim[1], mComplex[1]);

    uint32_t initLevel = 6;

    Plaintext ptxt0 = cc->MakeCKKSPackedPlaintext(mComplex[0], 1, maxDepth - initLevel);
    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(mComplex[1], 1, maxDepth - initLevel);
    Ciphertext<DCRTPoly> ctxt0 = cc->Encrypt(keys.publicKey, ptxt0);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keys.publicKey, ptxt1);
    auto b_start = std::chrono::high_resolution_clock::now();

    // ctxt0 = cc->EvalBootstrap(ctxt0);
    // ctxt1 = cc->EvalBootstrap(ctxt1);

    auto b_end = std::chrono::high_resolution_clock::now();
    auto b_duration = std::chrono::duration_cast<std::chrono::milliseconds>(b_end - b_start);
    std::cout << "Time taken (boot 2 times): " << std::dec << b_duration.count() << " ms" << std::endl;
    
    std::vector<Ciphertext<DCRTPoly>> ctxtout;

    // Evaluate the LUT and measure the time
    auto start = std::chrono::high_resolution_clock::now();

    EvalLUTs2DAsymmetric(ctxt0, ctxt1, logDim[0], logDim[1], coeffsPt, ctxtout);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time taken (" << logDim[0] << "+" << logDim[1] << ")-to-(" << logDim[0] << "+" << logDim[1] << ") LUT): " << std::dec << duration.count() << " ms" << std::endl;

    std::vector<Plaintext> result(outB);
    for (size_t i = 0; i < outB; i++){
        cc->Decrypt(keys.secretKey, ctxtout[i], &(result[i]));
    }

    std::vector<std::vector<int>> intResult(outB);
    for (size_t i = 0; i < outB; i++){
        std::vector<std::complex<double>> complexResult = result[i]->GetCKKSPackedValue();
        intResult[i] = GetLog(complexResult, logDim[i]);
    }

    std::cout << "Output table: \n====================\n";
    std::cout << "In(L)\tIn(R)\tMerged\t|Out(L)\tOut(R)\tMerged\tExpected\n";
    std::cout << std::hex;
    for (size_t i = 0; i < 8; i++){
        std::cout << mInt[1][i] << "\t" << mInt[0][i] << "\t" << mInt[1][i]*(1<<logDim[0]) + mInt[0][i] << 
                "\t|" << intResult[1][i] << "\t" << intResult[0][i] << "\t" << intResult[1][i]*(1<<logDim[0]) + intResult[0][i] <<
                "\t" << input_table[mInt[1][i]*(1<<logDim[0]) + mInt[0][i]]<< std::endl;
    }

    std::vector<std::vector<uint32_t>> expected(outB);
    nShift = 0;
    for (size_t i = 0; i < outB; i++){
        expected[i].resize(numSlots);
        uint32_t mask = (1<<logDim[i]) - 1;
        mask <<= nShift;
        for (size_t j = 0; j < numSlots; j++){
            expected[i][j] = (input_table[j % (input_table.size())] & mask) >> (nShift);
        }
        nShift += logDim[i];
    }

    std::cout << "Int and complex value table\n================================" << std::endl;
    std::cout << "Enced(int)\tEnced\tExpect(int)\tExpect" << std::endl;
    for (size_t i = 0; i < 8; i++){
        std::cout << intResult[0][i] << "\t" << result[0]->GetCKKSPackedValue()[i] << "\t";
        // calculate and print exp(-2*pi* i / 16 *expected[0][i])
        std::cout << expected[0][i] << "\t" << std::exp(std::complex<double>(0, -2 * M_PI * expected[0][i] / (1 << logDim[0]))) << std::endl;
    }

    uint32_t count_failure = 0;
    std::cout << std::dec;
    for (size_t i = 0; i < outB; i++){
        for (size_t j = 0; j < numSlots; j++){
            if(expected[i][j] != intResult[i][j]){
                std::cout << "Mismatch at " << i << " " << j << " " << expected[i][j] << " " 
                << intResult[i][j] << "(" << result[i]->GetCKKSPackedValue()[j] << ")" << std::endl;
                count_failure++;
            }
        }  
    }
    std::cout << "Number of failures: " << count_failure << "/" << numSlots*outB << std::endl;

    std::cout << "Noise variance is: ";
    double variance = 0.;
    for (size_t i = 0; i < numSlots; i++){
        // variance of noise, noise: result[0]->GetCKKSPackedValue()[i] - std::exp(std::complex<double>(0, -2 * M_PI * expected[0][i] / 16))
        variance += std::norm(result[0]->GetCKKSPackedValue()[i] - std::exp(std::complex<double>(0, -2 * M_PI * (int)expected[0][i] / (1 << logDim[0]))));
        variance += std::norm(result[1]->GetCKKSPackedValue()[i] - std::exp(std::complex<double>(0, -2 * M_PI * (int)expected[1][i] / (1 << logDim[1]))));
    }
    std::cout << variance / (2*numSlots) << std::endl;
}

void testExponent(){
    std::vector<uint32_t> mInt = {0,1,2,3,4,5,6,7};
    std::vector<std::complex<double>> mComplex;
    GetExponent(mInt, 3, mComplex);
    auto mInt2 = GetLog(mComplex, 3);
    for (size_t i = 0; i < mInt.size(); i++){
        std::cout << mInt[i] << " " << mComplex[i] << " "  << mInt2[i] << std::endl;
    }
}

void testLog(){
    // make some vector of complex numbner (size 4) and test getLog
    std::vector<std::complex<double>> mComplex = {std::complex<double>(0.000110165,0.999442), std::complex<double>(0, 1), std::complex<double>(-1, 0), std::complex<double>(0, -1)};
    auto mInt = GetLog(mComplex, 2);
    for (size_t i = 0; i < mComplex.size(); i++){
        std::cout << mComplex[i] << " " << mInt[i] << std::endl;
    }
}
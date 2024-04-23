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
#include <fftw3.h>
#include <chrono>


using namespace lbcrypto;

// gets a vector of integer (0 to codedim-1), and codedim, and return exp(-2 pi i * x / codedim)
void GetExponent(std::vector<uint32_t> input, size_t logDim, std::vector<std::complex<double>> &result);
// inverse of GetExponent
std::vector<int32_t> GetLog(std::vector<std::complex<double>> intput, size_t logDim);
void FindLUTPoly(std::vector<uint32_t> table, size_t logDim, 
                 std::vector<std::complex<double>> &coeff); 
void EvalLUT(Ciphertext<DCRTPoly> &ct, size_t logDim, std::vector<std::complex<double>> coeffs, 
             std::vector<Ciphertext<DCRTPoly>> &dest);
void findLUTPoly2D(std::vector<uint32_t> table, size_t logDim, 
                   std::vector<std::vector<std::complex<double>>> &coeff);

void EvalNoiseReductionInplace(Ciphertext<DCRTPoly> ct, size_t logDim, CryptoContext<DCRTPoly> cc);
void ifft(std::vector<std::complex<double>>& data);
void ifft2(std::vector<std::vector<std::complex<double>>>& data);

void testLUTNtoN();
void testLUTNtoBN();
void testLUTANtoBN();

void testFFT();
void testExponent();
void testComplexOperation();

int main() {
    // testLUTNtoN();
    // testLUTNtoBN();
    testLUTANtoBN();
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
    for (size_t i = 0; i < intput.size(); i++){
        double angle = std::arg(intput[i]);
        // if (angle < 0){
        //     angle += 2 * M_PI;
        // }
        result.push_back((codedim - 1) & (int32_t)std::round(angle * codedim / (-2 * M_PI)));
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
        GetExponent(std::vector<uint32_t>(table.begin() + i*codedim, table.begin() + (i+1)*codedim), logDim,
                    coeff[i]);
    }
    // Perform inverse FFT on coeff
    ifft2(coeff);
}

inline Ciphertext<DCRTPoly> EvalMultComplex(Ciphertext<DCRTPoly> &ct, 
                       std::complex<double> val){
    auto cc = ct->GetCryptoContext();
    return cc->EvalMult(ct, cc->MakeCKKSPackedPlaintext(std::vector<std::complex<double>>(cc->GetRingDimension()/2, val)));
}

inline Ciphertext<DCRTPoly> EvalAddComplex(Ciphertext<DCRTPoly> &ct, 
                       std::complex<double> val){
    auto cc = ct->GetCryptoContext();
    return cc->EvalAdd(ct, cc->MakeCKKSPackedPlaintext(std::vector<std::complex<double>>(cc->GetRingDimension()/2, val)));
}

// evaluate n-to-n LUT
void EvalLUT(Ciphertext<DCRTPoly> &ct, size_t logDim, 
             std::vector<std::complex<double>> coeffs, 
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
    dest = EvalMultComplex(b[0], coeffs[1]);
    for (size_t i = 2; i < codedim; i++){
        auto tmp = EvalMultComplex(b[i - 1], coeffs[i]);
        dest = cc->EvalAdd(dest, tmp);
    }
    dest = EvalAddComplex(dest, coeffs[0]);
}

// evaluate n-to-b*n LUT
void EvalLUTs(Ciphertext<DCRTPoly> &ct, size_t logDim, 
             std::vector<std::vector<std::complex<double>>> coeffs, 
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

    for(auto &table : coeffs){
        // inner product with coefficients
        dest.push_back(EvalMultComplex(b[0], table[1]));
        for (size_t i = 2; i < codedim; i++){
            auto tmp = EvalMultComplex(b[i - 1], table[i]);
            dest.back() = cc->EvalAdd(dest.back(), tmp);
        }
        dest.back() = EvalAddComplex(dest.back(), table[0]);
        std::cout << "size: " << dest.size() << std::endl;
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

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(mComplex);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keys.publicKey, ptxt1);

    // Evaluate the LUT and measure the time
    auto start = std::chrono::high_resolution_clock::now();

    EvalLUT(ctxt1, logDim, coeffs, ctxt1);

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
    std::cout << "Test single LUT" << std::endl;
    // Note: put any logDim-to-logDim table you want to test here 
    // The output range is dim^2, e.g, logDim = 4, the output range is 16^2
    std::vector<uint32_t> table(1<<logDim);
    for (size_t i = 0; i < table.size(); i++){
        table[i] = (int)(128*sin(i))+128; // Any cool function 
    }
    std::cout << "Input table: \n====================\n";
    for (size_t i = 0; i < table.size(); i++){
        std::cout << i << "\t" << table[i] << std::endl;
    }

    std::vector<std::vector<std::complex<double>>> coeffs2D(2);
    // divide tables by chunks of size logDim-bit
    for (size_t i = 0; i < 2; i++){
        std::vector<uint32_t> partialTable(table.size());
        uint32_t mask = (1<<logDim) - 1;
        mask <<= logDim * i;
        for (size_t j = 0; j < table.size(); j++){
            partialTable[j] = (table[j] & mask) >> (logDim * i);
        }
        FindLUTPoly(partialTable, logDim, coeffs2D[i]);
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

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(mComplex);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keys.publicKey, ptxt1);
    std::vector<Ciphertext<DCRTPoly>> ctxtout;
    // Evaluate the LUT and measure the time
    auto start = std::chrono::high_resolution_clock::now();

    EvalLUTs(ctxt1, logDim, coeffs2D, ctxtout);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Time taken: " << duration.count() << "ms" << std::endl;

    std::vector<Plaintext> result(2);
    for (size_t i = 0; i < 2; i++){
        cc->Decrypt(keys.secretKey, ctxtout[i], &result[i]);
    }

    std::vector<std::vector<int>> intResult(2);
    for (size_t i = 0; i < 2; i++){
        std::vector<std::complex<double>> complexResult = result[i]->GetCKKSPackedValue();
        intResult[i] = GetLog(complexResult, logDim);
    }

    std::cout << "Output table: \n====================\n";
    std::cout << "Input\tOut(R)\tOut(L)\tMerged\n";
    for (size_t i = 0; i < mInt.size(); i++){
        std::cout << mInt[i] << "\t" << intResult[0][i] << "\t" << intResult[1][i] << "\t" << intResult[0][i]+ intResult[1][i]*(1<<logDim)  << std::endl;
    }
}


void testLUTANtoBN(){
    size_t logDim = 4;
    std::cout << "Test single LUT" << std::endl;
    // Note: put any logDim-to-logDim table you want to test here 
    // The output range is dim^2, e.g, logDim = 4, the output range is 16^2
    std::vector<uint32_t> table(1<<logDim);
    for (size_t i = 0; i < table.size(); i++){
        table[i] = (int)(128*sin(i))+128; // Any cool function 
    }
    std::cout << "Input table: \n====================\n";
    for (size_t i = 0; i < table.size(); i++){
        std::cout << i << "\t" << table[i] << std::endl;
    }

    std::vector<std::vector<std::complex<double>>> coeffs2D(2);
    // divide tables by chunks of size logDim-bit
    for (size_t i = 0; i < 2; i++){
        std::vector<uint32_t> partialTable(table.size());
        uint32_t mask = (1<<logDim) - 1;
        mask <<= logDim * i;
        for (size_t j = 0; j < table.size(); j++){
            partialTable[j] = (table[j] & mask) >> (logDim * i);
        }
        FindLUTPoly(partialTable, logDim, coeffs2D[i]);
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

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(mComplex);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keys.publicKey, ptxt1);
    std::vector<Ciphertext<DCRTPoly>> ctxtout;
    // Evaluate the LUT and measure the time
    auto start = std::chrono::high_resolution_clock::now();

    EvalLUTs(ctxt1, logDim, coeffs2D, ctxtout);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Time taken: " << duration.count() << "ms" << std::endl;

    std::vector<Plaintext> result(2);
    for (size_t i = 0; i < 2; i++){
        cc->Decrypt(keys.secretKey, ctxtout[i], &result[i]);
    }

    std::vector<std::vector<int>> intResult(2);
    for (size_t i = 0; i < 2; i++){
        std::vector<std::complex<double>> complexResult = result[i]->GetCKKSPackedValue();
        intResult[i] = GetLog(complexResult, logDim);
    }

    std::cout << "Output table: \n====================\n";
    std::cout << "Input\tOut(R)\tOut(L)\tMerged\n";
    for (size_t i = 0; i < mInt.size(); i++){
        std::cout << mInt[i] << "\t" << intResult[0][i] << "\t" << intResult[1][i] << "\t" << intResult[0][i]+ intResult[1][i]*(1<<logDim)  << std::endl;
    }
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

void testComplexOperation(){
    std::vector<uint32_t> mInt = {0,1,2,3,4,5,6,7};
    std::vector<std::complex<double>> mComplex;
    GetExponent(mInt, 3, mComplex);
    std::complex<double> c(1,2);

    // Create the crypto context
    uint32_t multDepth = 1;
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

    // encrypt mInt
    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(mComplex);
    Ciphertext<DCRTPoly> ctxt1 = cc->Encrypt(keys.publicKey, ptxt1);

    auto ctxt2 = EvalAddComplex(ctxt1, c);

    // Decrypt ctxt2 and check the result
    Plaintext result;
    cc->Decrypt(keys.secretKey, ctxt2, &result);

    std::vector<std::complex<double>> complexResult = result->GetCKKSPackedValue();
    for (size_t i = 0; i < mInt.size(); i++){
        std::cout << mComplex[i] << " " << complexResult[i] << std::endl;
    }
}
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
 * This file benchmarks FHEW-GINX gate evaluation operations
 */

#include "benchmark/benchmark.h"
#include "binfhecontext.h"

using namespace lbcrypto;

/*
 * Context setup utility methods
 */

BinFHEContext GenerateFHEWContext(usint numAutoKeys) {
    BinFHEContextParams lmkcdey_param = 
        {   28,     2048,      447, 2048,   16384, 3.19,  32,      1024,  64,       numAutoKeys,        GAUSSIAN };
    auto cc = BinFHEContext();
    cc.GenerateBinFHEContext(lmkcdey_param, LMKCDEY);
    return cc;
}

/*
 * FHEW benchmarks
 */

void FHEW_BTKEYGEN(benchmark::State& state, usint numAutoKeys) {
    BinFHEContext cc = GenerateFHEWContext(numAutoKeys);

    for (auto _ : state) {
        LWEPrivateKey sk = cc.KeyGen();
        cc.BTKeyGen(sk);
    }
}

// function name, title to show, arguments
BENCHMARK_CAPTURE(FHEW_BTKEYGEN, STD128_LMKCDEY, 1u)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BTKEYGEN, STD128_LMKCDEY, 2u)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BTKEYGEN, STD128_LMKCDEY, 4u)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BTKEYGEN, STD128_LMKCDEY, 5u)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BTKEYGEN, STD128_LMKCDEY, 6u)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BTKEYGEN, STD128_LMKCDEY, 7u)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BTKEYGEN, STD128_LMKCDEY, 8u)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BTKEYGEN, STD128_LMKCDEY, 9u)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BTKEYGEN, STD128_LMKCDEY, 10u)->Unit(benchmark::kMicrosecond)->Iterations(50);

// benchmark for binary gates, such as AND, OR, NAND, NOR
template <class BinGate>
void FHEW_BINGATE(benchmark::State& state, usint numAutoKeys, BinGate bin_gate) {
    BINGATE gate(bin_gate);

    BinFHEContext cc = GenerateFHEWContext(numAutoKeys);

    LWEPrivateKey sk = cc.KeyGen();

    cc.BTKeyGen(sk);

    LWECiphertext ct1 = cc.Encrypt(sk, 1);
    LWECiphertext ct2 = cc.Encrypt(sk, 1);

    for (auto _ : state) {
        LWECiphertext ct11 = cc.EvalBinGate(gate, ct1, ct2);
    }
}

BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 1u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 2u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 3u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 4u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 5u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 6u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 7u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 8u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 9u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);
BENCHMARK_CAPTURE(FHEW_BINGATE, STD128_LMKCDEY_NAND, 10u, NAND)->Unit(benchmark::kMicrosecond)->Iterations(50);


BENCHMARK_MAIN();

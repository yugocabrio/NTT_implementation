# Rust test

## Introduction

Implement and benchmark in Rust a routine that performs the forward and backward Number Theoretic Transform (NTT) over $Z_{Q}[X]/X^N+1$ where $N$ is a power of two and $Q$ is a prime satisfying $Q\equiv 1\mod 2N$.

We will judge the submission by i) performance ii) code readability and iii) documentation.

### Base

Provide a generic implementation for the NTT for any prime fitting in an `u64` word.
Since do not expect the code to provide an API to find suitable primitive roots for a given prime, we provide the prime `0x1fffffffffe00001` and it's 2N-th root of unity `0x15eb043c7aa2b01f` for $N=2^{16}$ (hard coded in the template, feel free to change it).

### Bonus 1

Improve the generic implementation by extending it with platform specific optimize code (e.g. for an Intel CPUs with AVX2 or AVX512 extensions). It is ok if it requires bounding the size of the prime (e.g. 61 or 51 bits depending on the optimization).

### Bonus 2

Provide an implementation optimized for the Goldilocks prime `0xffffffff00000001` = $2^{64} - 2^{32} + 1$ with the 2N-th root of unity `0xabd0a6e8aa3d8a0e` for $N=2^{16}$.
For additional information check this link [](). 

### Resources for NTT
- [Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography](https://eprint.iacr.org/2016/504)
- [Number Theoretic Transform and Its Applications in Lattice-based Cryptosystems: A Survey](https://arxiv.org/pdf/2211.13546)





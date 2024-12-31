# Rust test

## Introduction

Implement and benchmark a rust routine that performs the forward and backward Number Theoretic Transform (NTT) over $Z_{Q}[X]/X^N+1$ where $N$ is a power of two and $Q$ is an NTT friendly prime satisfying $Q\equiv 1\mod 2N$.

We will judge the submission by i) performance ii) code readability and iii) documentation. We highly suggest adding unit test to ensure the solution performs correctly.

### Base task

Implement NTT for prime $q \leq 2^{61}$.

As a reference, we've provided a 61 bit prime `0x1fffffffffe00001` and its 2N-th root of unity `0x15eb043c7aa2b01f`, required for NTT, for $N=2^{16}$. Note that it's not necessary to stick with reference values.

To calculate a 2N-th root of unity for an NTT friendly prime, refer to the following [link](https://crypto.stackexchange.com/a/63616).

### Bonus 1

Provide a second implementation that leverages CPU specific instructions (e.g. Intel CPUs with AVX2 or AVX512 extensions). It's acceptable to limit bit-width of primes (for ex, to $\lt 61$ or $\lt 51$) if needed.

### Bonus 2

Provide a third implementation optimized for the Goldilocks prime `0xffffffff00000001` = $2^{64} - 2^{32} + 1$ (2N-th root of unity `0xabd0a6e8aa3d8a0e`). Check this [link](https://cp4space.hatsya.com/2021/09/01/an-efficient-prime-for-number-theoretic-transforms/) for additional information.

### Resources for NTT

-   [Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography](https://eprint.iacr.org/2016/504)
-   [Number Theoretic Transform and Its Applications in Lattice-based Cryptosystems: A Survey](https://arxiv.org/pdf/2211.13546)

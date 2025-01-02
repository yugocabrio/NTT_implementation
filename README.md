# Rust test

## Introduction

This is the coding interview for the Rust Engineer role at Phantom Zone (see [Linkedin](https://www.linkedin.com/jobs/view/4110239195/)).

Implement and benchmark a rust routine that performs the forward and backward Number Theoretic Transform (NTT) over $Z_{Q}[X]/X^N+1$ where $N$ is a power of two and $Q$ is an NTT friendly prime satisfying $Q\equiv 1\mod 2N$.

We will judge the submission by i) performance ii) code readability and iii) documentation. We highly suggest adding unit test to ensure the solution performs correctly.

### Base task

Implement NTT for prime $q \leq 2^{61}$.

As a reference, we've provided a 61 bit prime `0x1fffffffffe00001` and its 2^17-th root of unity `0x15eb043c7aa2b01f`, required for NTT, for $N=2^{16}$. Note that it's not necessary to stick with reference values.

To calculate a 2N-th root of unity for an NTT friendly prime, refer to the following [link](https://crypto.stackexchange.com/a/63616).

### Bonus 1

Provide a second implementation that leverages CPU specific instructions (e.g. Intel CPUs with AVX2 or AVX512 extensions). It's acceptable to limit bit-width of primes (for ex, to $\lt 61$ or $\lt 51$ or even $\lt 32$) if needed.

This implementation is expected to be faster than the base case.

### Bonus 2

Provide a third implementation optimized for the Goldilocks prime `0xffffffff00000001` = $2^{64} - 2^{32} + 1$ (2^17-th root of unity `0xabd0a6e8aa3d8a0e`). Check this [link](https://cp4space.hatsya.com/2021/09/01/an-efficient-prime-for-number-theoretic-transforms/) for additional information.

This implementation is expected to be faster than one accepting a generic prime (either w.r.t. the base case or the vectorized implementation).

## Validity Conditions

A solution will **only** be valid if:
- It performs correctly (the code correctly evaluates the NTT and INTT).
- It does so efficiently, even for the base case (e.g. by using specialized finite field arithmetic).
- It was written solely by the candidate and can be fully explained by the later.
- It is licensed under Apache 2.0.

Bonuses 1 and 2 are not needed for a submission to be valid, but are a necessary condition to be eligible for the compensation (see job post).

### Resources for NTT

-   [Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography](https://eprint.iacr.org/2016/504)
-   [Number Theoretic Transform and Its Applications in Lattice-based Cryptosystems: A Survey](https://arxiv.org/pdf/2211.13546)

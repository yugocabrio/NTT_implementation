# Rust NTT Library

This repository provides several implementations of the Number Theoretic Transform (NTT) in Rust.  
In particular, it implements the negacyclic version of the NTT over the ring $\mathbb{Z}_{Q}[X]/\bigl(X^{N} + 1\bigr)$, where $N$ is a power of two and $Q$ is an NTT-friendly prime satisfying $Q \equiv 1 \pmod{2N}$.

Note: This code was created for educational purposes by a learner.

---

## Type of NTTs

1. `shoup_ntt.rs`
   - A general-purpose NTT for up to $64$-bit primes, using Shoup multiplication.  
   - On macOS, its performance is roughly comparable to concrete-ntt.

2. `mont_ntt.rs` 
   - A general-purpose NTT for up to $64$-bit primes, using Montgomery multiplication.  
   - Currently has some overhead due to converting to/from the Montgomery form.

3. `goldilocks_ntt.rs` 
   - Specialized for the so-called Goldilocks prime $2^{64} - 2^{32} + 1$ amd used PZT22 reduction.  
   - Its performance is roughly on par with the plonky2 implementation.

4. `barrett_scalar_ntt.rs` 
   - A general-purpose NTT for up to $32$-bit primes, using Barrett multiplication.

5. `barrett_vector_ntt.rs`
   - A general-purpose NTT for up to $32$-bit primes, partially vectorized with NEON instructions (on AArch64) for addition and subtraction.  
   - Currently, its performance is similar to the scalar version because `barrett_mul` is still done in scalar form. There is a plan to vectorize `barrett_mul` as well.

---

## Benchmark

### Run
```
cargo bench
```
Note that `barrett_vector.rs` code runs only on target_arch = "aarch64".

Here is a speed comparison of `forward_inplace`.

### A general-purpose NTT for up to 64-bit prime  
| log_n | `mont_ntt.rs`| `shoup_ntt.rs`| `concrete_ntt` |
|-------|----------:|----------:|-------------:|
| 11    | 50.197µs  | 10.585µs  | 10.623µs     |
| 12    | 112.19µs  | 23.880µs  | 22.586µs     |
| 13    | 225.49µs  | 49.489µs  | 48.661µs     |
| 14    | 475.96µs  | 106.66µs  | 106.25µs     |
| 15    | 996.50µs  | 240.23µs  | 230.52µs     |
| 16    | 2.0449ms  | 477.64µs  | 477.72µs     |

(from `benches/ntt.rs`)

### A NTT for goldilocks prime  
| log_n | `goldilocks_ntt.rs` | `plonky2` |
|-------|---------:|----------:|
| 11    | 23.192µs | 31.813µs  |
| 12    | 50.804µs | 63.902µs  |
| 13    | 108.68µs | 128.01µs  |
| 14    | 251.59µs | 281.26µs  |
| 15    | 510.63µs | 591.31µs  |
| 16    | 1.1101ms | 1.2720ms  |

(from `benches/ntt.rs`)

### A NTT for 32 bit prime(Neon/SIMD experiment) 

| log_n | `barrett_scalar.rs` (un-optimized add/sub) | `barrett_scalar.rs` (optimized add/sub) | `barrett_vector.rs` (vectorized add/sub) | `concrete_ntt` |
|-------|-------------------------:|--------------------------:|----------------------:|--------------------:|
| 11    | 17.583µs                | 10.027µs                  | 10.403µs             | 9.5152µs           |
| 12    | 37.852µs                | 24.04µs                   | 23.507µs             | 20.331µs           |
| 13    | 82.653µs                | 54.333µs                  | 49.696µs             | 39.564µs           |
| 14    | 177.18µs                | 96.052µs                  | 102.73µs             | 82.942µs           |
| 15    | 432.07µs                | 208.71µs                  | 215.82µs             | 166.73µs           |
| 16    | 878.98µs                | 427.84µs                  | 457.99µs             | 347.83µs           |

(from `benches/vectorized.rs`)

---

## Test
```
cargo test
```
This command runs unit tests for various finite field arithmetic, as well as the NTT round-trip and polynomial multiplication tests.

---

## TODO

- [ ] Address any warnings from `cargo clippy`.
- [ ] Apply parallelization to the butterfly algorithm.
use crate::dft::barrett_field_32bit::{
    add as add_32, exp as exp_32, inv as inv_32, mul as mul_32, sub as sub_32,
};
use crate::dft::field::{
    add as add_64, exp as exp_64, inv as inv_64, mul as mul_64, sub as sub_64,
};
use rand::Rng;

/// Finds a primitive 2n-th root of unity modulo p for NWC based NTT by random sampling.
/// Checks (p - 1) is divisible by 2n, then construct g = x^(p-1/2n) (mod p) with ramdom x.
/// Checks if g^n is -1 (mod p) and g^(2n) is 1 (mod p).
/// It returns (g, g_inv)
#[inline(always)]
pub fn find_primitive_2nth_root_of_unity_64(p: u64, n: usize) -> Option<(u64, u64)> {
    let two_n = 2 * (n as u64);
    if (p - 1) % two_n != 0 {
        return None;
    }
    let exponent = (p - 1) / two_n;
    let mut rng = rand::thread_rng();

    for _ in 0..100 {
        let x_rand = rng.gen_range(1..p);
        let g = exp_64(x_rand, exponent, p);
        // g^n = p-1 => -1 mod p
        let g_n = exp_64(g, n as u64, p);
        if g_n == p.wrapping_sub(1) {
            // g^(2n) = 1 mod p
            let g_2n = exp_64(g, 2 * (n as u64), p);
            if g_2n == 1 {
                if let Some(g_inv) = inv_64(g, p) {
                    return Some((g, g_inv));
                }
            }
        }
    }
    None
}

#[inline(always)]
pub fn find_primitive_2nth_root_of_unity_32(p: u32, n: usize) -> Option<(u32, u32)> {
    let two_n = 2 * (n as u32);
    if (p - 1) % two_n != 0 {
        return None;
    }
    let exponent = (p - 1) / two_n;
    let mut rng = rand::thread_rng();

    for _ in 0..100 {
        let x = rng.gen_range(1..p);
        let g = exp_32(x, exponent, p);
        let g_n = exp_32(g, n as u32, p);
        if g_n == p - 1 {
            let g_2n = exp_32(g, 2 * (n as u32), p);
            if g_2n == 1 {
                if let Some(g_inv) = inv_32(g, p) {
                    return Some((g, g_inv));
                }
            }
        }
    }
    None
}

/// Builds bit-reversed twiddle-factor tables.
/// fwd[i] = ψ^i mod q (with bit-reversed index mapping)
/// inv[i] = ψ_inv^i mod q
#[inline(always)]
pub fn build_bitrev_tables_u64(n: usize, q: u64, psi: u64, psi_inv: u64) -> (Vec<u64>, Vec<u64>) {
    let mut fwd = vec![0u64; n];
    let mut inv = vec![0u64; n];
    let log_n = n.trailing_zeros();

    let mut cur_f = 1u64;
    let mut cur_i = 1u64;

    for i in 0..n {
        let r = (i as u32).reverse_bits() >> (32 - log_n);
        let r = r as usize;

        fwd[r] = cur_f;
        inv[r] = cur_i;

        // (psi^(i+1)), (psi_inv^(i+1))
        cur_f = mul_64(cur_f, psi, q);
        cur_i = mul_64(cur_i, psi_inv, q);
    }
    (fwd, inv)
}

#[inline(always)]
pub fn build_bitrev_tables_u32(n: usize, q: u32, psi: u32, psi_inv: u32) -> (Vec<u32>, Vec<u32>) {
    let mut fwd = vec![0u32; n];
    let mut inv = vec![0u32; n];
    let log_n = n.trailing_zeros();

    let mut cur_f = 1u32;
    let mut cur_i = 1u32;

    for i in 0..n {
        let r = (i as u32).reverse_bits() >> (32 - log_n);
        let r = r as usize;

        fwd[r] = cur_f;
        inv[r] = cur_i;

        // (psi^(i+1)), (psi_inv^(i+1))
        cur_f = mul_32(cur_f, psi, q);
        cur_i = mul_32(cur_i, psi_inv, q);
    }
    (fwd, inv)
}

/// Naive negacyclic polynomial multiplication in 64-bit arithmetic.
/// Interprets x^n = -1 and subtracts terms when the index exceeds n.
/// This is used for the NTT-based polynomial multiplication test.
#[inline]
pub fn naive_negacyclic_u64(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
    let n = a.len();
    let mut c = vec![0u64; n];

    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            let prod = mul_64(ai, bj, p);
            let idx = i + j;

            // If idx < n, just add to c[idx].
            // If idx >= n, we interpret x^n = -1, so c[idx - n] -= product.
            if idx < n {
                c[idx] = add_64(c[idx], prod, p);
            } else {
                c[idx - n] = sub_64(c[idx - n], prod, p);
            }
        }
    }
    c
}

#[inline]
pub fn naive_negacyclic_u32(a: &[u32], b: &[u32], p: u32) -> Vec<u32> {
    let n = a.len();
    let mut c = vec![0u32; n];

    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            let prod = mul_32(ai, bj, p);
            let idx = i + j;

            if idx < n {
                c[idx] = add_32(c[idx], prod, p);
            } else {
                c[idx - n] = sub_32(c[idx - n], prod, p);
            }
        }
    }
    c
}

/// Performs pointwise multiplication on two vectors under mod p.
/// This is used for the NTT-based polynomial multiplication test.
#[inline]
pub fn pointwise_u64(a: &mut [u64], b: &[u64], p: u64) {
    for i in 0..a.len() {
        a[i] = mul_64(a[i], b[i], p);
    }
}

#[inline]
pub fn pointwise_u32(a: &mut [u32], b: &[u32], p: u32) {
    for i in 0..a.len() {
        let prod = (a[i] as u64) * (b[i] as u64) % (p as u64);
        a[i] = prod as u32;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_naive_negacyclic_64_simple() {
        let p = 19u64;
        let a = [1u64, 2, 3, 4];
        let b = [2u64, 2, 2, 2];
        let c = naive_negacyclic_u64(&a, &b, p);
        assert_eq!(c, [3, 11, 4, 1]);
    }

    #[test]
    fn test_pointwise_64_simple() {
        let p = 19u64;
        let mut a = [3u64, 5];
        let b = [7u64, 2];
        pointwise_u64(&mut a, &b, p);
        assert_eq!(a, [2, 10]);
    }

    #[test]
    fn test_find_primitive_2n_root_of_unity() {
        let q = 7681u64;
        let n = 16;

        let got = find_primitive_2nth_root_of_unity_64(q, n);
        assert!(got.is_some(), "cannot find");
        let (g, g_inv) = got.unwrap();

        let inv_check = mul_64(g, g_inv, q);
        assert_eq!(inv_check, 1, "g*g_inv != 1");

        let g_2n = exp_64(g, (2 * n) as u64, q);
        assert_eq!(g_2n, 1, "g^(2n) != 1 mod q");

        let g_n = exp_64(g, n as u64, q);
        assert_ne!(g_n, 1, "g^n == 1 mod q");
    }
}

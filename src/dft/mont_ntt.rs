use crate::dft::field::inv;
use crate::dft::mont_field::{mont_add, mont_inv, mont_mul, mont_sub, MontgomeryContext};
use crate::dft::util::{build_bitrev_tables_u64, find_primitive_2nth_root_of_unity_64};
use crate::dft::DFT;

/// A structure for performing NTT using Montgomery multiplication.
pub struct MontTable {
    /// NTT friendly prime
    q: u64,
    /// n, the power of 2
    n: usize,
    /// 2n-th root of unity
    psi: u64,
    /// inverse of psi
    psi_inv: u64,

    /// Bit-reversed twiddle factors
    fwd_twid: Vec<u64>,
    /// same but inverse
    inv_twid: Vec<u64>,

    /// n^-1 mod q in mont form (used for the backward transform)
    inv_n: u64,

    /// Montgomery context for operations mod q.
    mont: MontgomeryContext,
}

impl MontTable {
    /// Creates a default MontTable with q=0x1fffffffffe00001, n=2^16, and a known 2n-th root of unity.
    #[inline(always)]
    pub fn new() -> Self {
        let q = 0x1fffffffffe00001u64;
        let n = 1 << 16;
        // 2^17-th root of unity
        let psi = 0x15eb043c7aa2b01fu64;
        let psi_inv = inv(psi, q).expect("cannot calc invere of psi");

        // Build bitrev twiddle factors but not in Montgomery form yet
        let (mut fwd_twid, mut inv_twid) = build_bitrev_tables_u64(n, q, psi, psi_inv);

        let k: u32 = 62;
        let mont = MontgomeryContext::new(q, k);

        // Convert twiddle factors into Montgomery form
        for x in fwd_twid.iter_mut() {
            *x = mont.to_mont(*x);
        }
        for x in inv_twid.iter_mut() {
            *x = mont.to_mont(*x);
        }

        // Precompute n^-1 in Montgomery form
        let n_mont = mont.to_mont(n as u64);
        let inv_n = mont_inv(n_mont, &mont).expect("cannot calc inverse of n");

        Self {
            q,
            n,
            psi,
            psi_inv,
            fwd_twid,
            inv_twid,
            inv_n,
            mont,
        }
    }

    /// Builds a MontTable from user-specified prime q and size n,
    /// automatically finding a suitable 2n-th root of unity.
    #[inline(always)]
    pub fn with_params(q: u64, n: usize) -> Option<Self> {
        if !n.is_power_of_two() {
            return None;
        }
        if (q - 1) % (2 * n as u64) != 0 {
            return None;
        }

        let (psi, psi_inv) = find_primitive_2nth_root_of_unity_64(q, n)?;

        let (mut fwd_twid, mut inv_twid) = build_bitrev_tables_u64(n, q, psi, psi_inv);

        let k: u32 = 62;
        if q >= (1u64 << k) {
            return None;
        }
        let mont = MontgomeryContext::new(q, k);

        // Convert twiddle factors into Montgomery form
        for x in fwd_twid.iter_mut() {
            *x = mont.to_mont(*x);
        }
        for x in inv_twid.iter_mut() {
            *x = mont.to_mont(*x);
        }

        let n_mont = mont.to_mont(n as u64);
        let inv_n = mont_inv(n_mont, &mont)?;

        Some(Self {
            q,
            n,
            psi,
            psi_inv,
            fwd_twid,
            inv_twid,
            inv_n,
            mont,
        })
    }

    #[inline(always)]
    pub fn q(&self) -> u64 {
        self.q
    }
    #[inline(always)]
    pub fn size(&self) -> usize {
        self.n
    }

    #[inline(always)]
    pub fn forward_inplace(&self, a: &mut [u64]) {
        // Convert all a[i] into Montgomery form
        for x in a.iter_mut() {
            *x = self.mont.to_mont(*x);
        }

        let mut t = self.n;
        let mut m = 1;
        while m < self.n {
            t >>= 1;
            for i in 0..m {
                let j1 = 2 * i * t;
                let j2 = j1 + t - 1;
                let s_mont = self.fwd_twid[m + i];
                for j in j1..=j2 {
                    let u = unsafe { *a.get_unchecked(j) };
                    let tv = unsafe { *a.get_unchecked(j + t) };
                    let v = mont_mul(tv, s_mont, &self.mont);

                    let sum_ = mont_add(u, v, &self.mont);
                    let diff_ = mont_sub(u, v, &self.mont);

                    unsafe {
                        *a.get_unchecked_mut(j) = sum_;
                        *a.get_unchecked_mut(j + t) = diff_;
                    }
                }
            }
            m <<= 1;
        }
        // Convert back from Montgomery form
        for x in a.iter_mut() {
            *x = self.mont.from_mont(*x);
        }
    }

    #[inline(always)]
    pub fn backward_inplace(&self, a: &mut [u64]) {
        // Convert all a[i] into Montgomery form
        for x in a.iter_mut() {
            *x = self.mont.to_mont(*x);
        }

        let mut t = 1;
        let mut m = self.n;
        while m > 1 {
            let h = m >> 1;
            for i in 0..h {
                let j1 = 2 * i * t;
                let j2 = j1 + t - 1;
                let s_mont = self.inv_twid[h + i];
                for j in j1..=j2 {
                    let u = unsafe { *a.get_unchecked(j) };
                    let v = unsafe { *a.get_unchecked(j + t) };

                    let sum_ = mont_add(u, v, &self.mont);
                    let diff_ = mont_sub(u, v, &self.mont);
                    let diff_m = mont_mul(diff_, s_mont, &self.mont);

                    unsafe {
                        *a.get_unchecked_mut(j) = sum_;
                        *a.get_unchecked_mut(j + t) = diff_m;
                    }
                }
            }
            t <<= 1;
            m = h;
        }
        // Multiply each entry by inv_n (Montgomery form)
        for x in a.iter_mut() {
            *x = mont_mul(*x, self.inv_n, &self.mont);
        }
        // Convert back from Montgomery form
        for x in a.iter_mut() {
            *x = self.mont.from_mont(*x);
        }
    }
}

impl DFT<u64> for MontTable {
    #[inline(always)]
    fn forward_inplace(&self, a: &mut [u64]) {
        self.forward_inplace(a);
    }
    #[inline(always)]
    fn backward_inplace(&self, a: &mut [u64]) {
        self.backward_inplace(a);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dft::util::{naive_negacyclic_u64, pointwise_u64};
    use rand::thread_rng;
    use rand::Rng;

    #[test]
    fn forward_backward_small_prime() {
        let q = 7681u64;
        let n = 16usize;
        let table = MontTable::with_params(q, n).unwrap();

        let mut rng = thread_rng();
        let mut data = vec![0u64; n];
        for x in data.iter_mut() {
            *x = rng.gen_range(0..q);
        }
        let orig = data.clone();

        table.forward_inplace(&mut data);
        table.backward_inplace(&mut data);

        assert_eq!(data, orig);
    }

    #[test]
    fn forward_backward_default() {
        let table = MontTable::new();
        let q = table.q();
        let n = table.size();
        let mut rng = thread_rng();
        let mut data = vec![0u64; n];
        for x in data.iter_mut() {
            *x = rng.gen_range(0..q);
        }
        let orig = data.clone();
        table.forward_inplace(&mut data);
        table.backward_inplace(&mut data);
        assert_eq!(data, orig);
    }

    #[test]
    fn test_ntt_polymul_small_prime() {
        let q = 7681u64;
        let n = 8usize;

        let table = MontTable::with_params(q, n).unwrap();

        let mut rng = rand::thread_rng();
        let mut a = vec![0u64; n];
        let mut b = vec![0u64; n];
        for i in 0..n {
            a[i] = rng.gen_range(0..q);
            b[i] = rng.gen_range(0..q);
        }

        let result_naive = naive_negacyclic_u64(&a, &b, q);

        table.forward_inplace(&mut a);
        table.forward_inplace(&mut b);
        pointwise_u64(&mut a, &mut b, q);
        table.backward_inplace(&mut a);

        assert_eq!(a, result_naive);
    }

    #[test]
    #[ignore]
    fn test_ntt_polymul_default() {
        let table = MontTable::new();
        let q = table.q();
        let n = table.size();

        let mut rng = rand::thread_rng();
        let mut a = vec![0u64; n];
        let mut b = vec![0u64; n];
        for i in 0..n {
            a[i] = rng.gen_range(0..q);
            b[i] = rng.gen_range(0..q);
        }

        let result_naive = naive_negacyclic_u64(&a, &b, q);

        table.forward_inplace(&mut a);
        table.forward_inplace(&mut b);
        pointwise_u64(&mut a, &mut b, q);
        table.backward_inplace(&mut a);

        assert_eq!(a, result_naive);
    }
}

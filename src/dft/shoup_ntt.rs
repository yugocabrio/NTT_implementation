use crate::dft::field::{add, inv, mul, sub};
use crate::dft::shoup_field::{shoup_mul, shoup_precompute};
use crate::dft::util::find_primitive_2nth_root_of_unity_64;
use crate::dft::DFT;

type ShoupTwiddleTables = (Vec<(u64, u64)>, Vec<(u64, u64)>);

/// A structure for performing NTT using Shoup multiplication.
pub struct ShoupTable {
    /// NTT-friendly prime
    q: u64,
    /// n = power-of-two
    n: usize,

    /// Bit-reversed twiddle factors (forward). Each element is (w, w_shoup),
    /// so we can do (a * w) mod q via Shoup multiplication.
    fwd_twid: Vec<(u64, u64)>,
    /// Bit-reversed twiddle factors (inverse)
    inv_twid: Vec<(u64, u64)>,

    /// n^-1 in normal form + its shoup counterpart.
    /// Used for final scaling in inverse NTT.
    inv_n: (u64, u64),
}

impl Default for ShoupTable {
    fn default() -> Self {
        Self::new()
    }
}

impl ShoupTable {
    /// Creates a default ShoupTable with q = 0x1fffffffffe00001, n = 2^16, and a known 2n-th root of unity.
    /// It also sets up the (w, w_shoup) pairs in a bit-reversed order.
    #[inline]
    pub fn new() -> Self {
        let q = 0x1fffffffffe00001u64;
        let n = 1 << 16;
        let psi = 0x15eb043c7aa2b01fu64;
        let psi_inv = inv(psi, q).expect("cannot invert psi");

        // inv(n) in normal form, then precompute its shoup pair
        let inv_n_val = inv(n as u64, q).expect("cannot invert n");
        let inv_n_pair = shoup_precompute(inv_n_val, q);

        // Build bit-reversed tables of twiddle factors, each stored as (w, w_shoup).
        let (fwd_twid, inv_twid) = shoup_build_bitrev_tables(q, n, psi, psi_inv);

        Self {
            q,
            n,
            fwd_twid,
            inv_twid,
            inv_n: inv_n_pair,
        }
    }

    /// Builds a ShoupTable from user-specified prime q and size n,
    /// automatically finding a suitable 2n-th root of unity.
    /// Precomputes the (w, w_shoup) pairs for bit-reversed order.
    #[inline]
    pub fn with_params(q: u64, n: usize) -> Option<Self> {
        if !n.is_power_of_two() {
            return None;
        }
        if (q - 1) % (2 * n as u64) != 0 {
            return None;
        }

        // find psi
        let (psi, psi_inv) = find_primitive_2nth_root_of_unity_64(q, n)?;

        let inv_n_val = inv(n as u64, q)?;
        let inv_n_pair = shoup_precompute(inv_n_val, q);

        let (fwd_twid, inv_twid) = shoup_build_bitrev_tables(q, n, psi, psi_inv);

        Some(Self {
            q,
            n,
            fwd_twid,
            inv_twid,
            inv_n: inv_n_pair,
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
        let q = self.q;
        let mut half = self.n;
        let mut step = 1;

        while step < self.n {
            half >>= 1;
            for i in 0..step {
                let (w, w_shoup) = self.fwd_twid[step + i];
                let base = 2 * i * half;
                let end = base + half;

                for j in base..end {
                    let u = unsafe { *a.get_unchecked(j) };
                    let tv = unsafe { *a.get_unchecked(j + half) };

                    let v = shoup_mul(tv, (w, w_shoup), q);

                    let sum_ = add(u, v, q);
                    let diff_ = sub(u, v, q);

                    unsafe {
                        *a.get_unchecked_mut(j) = sum_;
                        *a.get_unchecked_mut(j + half) = diff_;
                    }
                }
            }
            step <<= 1;
        }
    }

    #[inline(always)]
    pub fn backward_inplace(&self, a: &mut [u64]) {
        let q = self.q;
        let (inv_n_val, inv_n_shoup) = self.inv_n;

        let mut step = self.n;
        let mut half = 1;
        while step > 1 {
            let halfstep = step >> 1;
            for i in 0..halfstep {
                let (w, w_shoup) = self.inv_twid[halfstep + i];
                let base = 2 * i * half;
                let end = base + half;

                for j in base..end {
                    let u = unsafe { *a.get_unchecked(j) };
                    let v = unsafe { *a.get_unchecked(j + half) };

                    let sum_ = add(u, v, q);
                    let diff_ = sub(u, v, q);
                    let diff_m = shoup_mul(diff_, (w, w_shoup), q);

                    unsafe {
                        *a.get_unchecked_mut(j) = sum_;
                        *a.get_unchecked_mut(j + half) = diff_m;
                    }
                }
            }
            half <<= 1;
            step = halfstep;
        }

        // multiply everything by inv_n via Shoup
        for x in a.iter_mut() {
            *x = shoup_mul(*x, (inv_n_val, inv_n_shoup), q);
        }
    }
}

impl DFT<u64> for ShoupTable {
    #[inline(always)]
    fn forward_inplace(&self, a: &mut [u64]) {
        self.forward_inplace(a);
    }

    #[inline(always)]
    fn backward_inplace(&self, a: &mut [u64]) {
        self.backward_inplace(a);
    }
}

/// Builds bit-reversed tables of (w, w_shoup) pairs for Shoup multiplication
#[inline]
fn shoup_build_bitrev_tables(q: u64, n: usize, psi: u64, psi_inv: u64) -> ShoupTwiddleTables {
    let log_n = n.trailing_zeros();
    let mut fwd = vec![(0u64, 0u64); n];
    let mut inv = vec![(0u64, 0u64); n];

    let mut power_psi = 1u64;
    let mut power_psi_inv = 1u64;

    for i in 0..n {
        let r = (i as u32).reverse_bits() >> (32 - log_n);
        let ridx = r as usize;

        fwd[ridx] = shoup_precompute(power_psi, q);
        inv[ridx] = shoup_precompute(power_psi_inv, q);

        power_psi = mul(power_psi, psi, q);
        power_psi_inv = mul(power_psi_inv, psi_inv, q);
    }
    (fwd, inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dft::util::{naive_negacyclic_u64, pointwise_u64};
    use rand::thread_rng;
    use rand::Rng;

    #[test]
    fn shoup_forward_backward_small_prime() {
        let q = 7681u64;
        let n = 16usize;
        let table = ShoupTable::with_params(q, n).unwrap();

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
    fn shoup_forward_backward_default() {
        let table = ShoupTable::new();
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
    fn shoup_test_ntt_polymul_small_prime() {
        let q = 7681u64;
        let n = 8usize;

        let table = ShoupTable::with_params(q, n).unwrap();

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
    fn shoup_test_ntt_polymul_default() {
        let table = ShoupTable::new();
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

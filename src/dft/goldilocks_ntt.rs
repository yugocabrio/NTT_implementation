use crate::dft::DFT;
use crate::dft::goldilocks_field::{GOLDILOCKS_P, mul_mod, add_mod, sub_mod, mod_inv, pow_mod};
use rand::Rng;

pub struct GoldilocksNttTable {
    q: u64,
    n: usize,
    psi: u64,
    psi_inv: u64,
    fwd_twid: Vec<u64>,
    inv_twid: Vec<u64>,
    inv_n: u64,
}

impl GoldilocksNttTable {
    /// デフォルト生成: n=2^16, ψ=0xabd0a6e8aa3d8a0e で固定
    pub fn new() -> Self {
        let q = GOLDILOCKS_P;
        let n = 1 << 16;
        let psi = 0xabd0a6e8aa3d8a0e_u64;
        let psi_inv = mod_inv(psi);
        let inv_n = mod_inv(n as u64);

        let (fwd_twid, inv_twid) = build_bitrev_tables(n, psi, psi_inv);
        Self { q, n, psi, psi_inv, fwd_twid, inv_twid, inv_n }
    }

    /// 任意の n=2^k に対して 2n-th root を探索してテーブル生成
    pub fn with_params(n: usize) -> Option<Self> {
        if !n.is_power_of_two() {
            return None;
        }
        let (psi, psi_inv) = find_primitive_2n_root_of_unity(n)?;
        let (fwd_twid, inv_twid) = build_bitrev_tables(n, psi, psi_inv);
        let inv_n = mod_inv(n as u64);
        Some(Self { q: GOLDILOCKS_P, n, psi, psi_inv, fwd_twid, inv_twid, inv_n })
    }

    #[inline(always)]
    pub fn q(&self) -> u64 { self.q }
    #[inline(always)]
    pub fn size(&self) -> usize { self.n }

    #[inline(always)]
    pub fn forward_inplace(&self, a: &mut [u64]) {
        let n = self.n;
        let mut half = n;
        let mut step = 1;
        while step < n {
            half >>= 1;
            for i in 0..step {
                let w = self.fwd_twid[step + i];
                let base = 2 * i * half;
                unsafe {
                    let a_ptr = a.as_mut_ptr();
                    for j in 0..half {
                        let idx = base + j;
                        let u = *a_ptr.add(idx);
                        let tv = *a_ptr.add(idx + half);
                        let v = mul_mod(tv, w);
                        *a_ptr.add(idx) = add_mod(u, v);
                        *a_ptr.add(idx + half) = sub_mod(u, v);
                    }
                }
            }
            step <<= 1;
        }
    }

    #[inline(always)]
    pub fn backward_inplace(&self, a: &mut [u64]) {
        let n = self.n;
        let mut step = n;
        let mut half = 1;
        while step > 1 {
            let halfstep = step / 2;
            for i in 0..halfstep {
                let w = self.inv_twid[halfstep + i];
                let base = 2 * i * half;
                unsafe {
                    let a_ptr = a.as_mut_ptr();
                    for j in 0..half {
                        let idx = base + j;
                        let u = *a_ptr.add(idx);
                        let v = *a_ptr.add(idx + half);
                        let sum_ = add_mod(u, v);
                        let diff_ = sub_mod(u, v);
                        let diffm = mul_mod(diff_, w);
                        *a_ptr.add(idx) = sum_;
                        *a_ptr.add(idx + half) = diffm;
                    }
                }
            }
            half *= 2;
            step = halfstep;
        }
        for x in a.iter_mut() {
            *x = mul_mod(*x, self.inv_n);
        }
    }
}

impl DFT<u64> for GoldilocksNttTable {
    #[inline(always)]
    fn forward_inplace(&self, a: &mut [u64]) {
        self.forward_inplace(a)
    }
    #[inline(always)]
    fn backward_inplace(&self, a: &mut [u64]) {
        self.backward_inplace(a)
    }
}

pub fn find_primitive_2n_root_of_unity(n: usize) -> Option<(u64, u64)> {
    let two_n = 2 * (n as u64);
    if (GOLDILOCKS_P - 1) % two_n != 0 {
        return None;
    }
    let exponent = (GOLDILOCKS_P - 1) / two_n;
    let mut rng = rand::thread_rng();
    const TRIES_LIMIT: usize = 500_000;

    for _ in 0..TRIES_LIMIT {
        let x = rng.gen_range(1..GOLDILOCKS_P);
        let g = pow_mod(x, exponent);
        // g^n = p-1  (≡ -1 mod p) かどうか
        if pow_mod(g, n as u64) == GOLDILOCKS_P.wrapping_sub(1) {
            // かつ g^(2n) = 1
            if pow_mod(g, 2 * (n as u64)) == 1 {
                return Some((g, mod_inv(g)));
            }
        }
    }
    None
}

fn build_bitrev_tables(n: usize, psi: u64, psi_inv: u64) -> (Vec<u64>, Vec<u64>) {
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
        cur_f = mul_mod(cur_f, psi);
        cur_i = mul_mod(cur_i, psi_inv);
    }
    (fwd, inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dft::DFT;
    use crate::dft::goldilocks_field::{mul_mod, add_mod, sub_mod, GOLDILOCKS_P};
    use rand::{thread_rng, Rng};

    fn poly_negacyclic_naive(a: &[u64], b: &[u64]) -> Vec<u64> {
        let n = a.len();
        let mut c = vec![0u64; n];
        for i in 0..n {
            for j in 0..n {
                let prod = mul_mod(a[i], b[j]);
                let idx = i + j;
                if idx < n {
                    c[idx] = add_mod(c[idx], prod);
                } else {
                    // x^n = -1
                    c[idx - n] = sub_mod(c[idx - n], prod);
                }
            }
        }
        c
    }

    #[test]
    fn goldi_test_forward_backward_n8() {
        let n = 8;
        let table = GoldilocksNttTable::with_params(n)
            .expect("cannot find it");
        let mut rng = thread_rng();

        let mut data = vec![0u64; n];
        for x in data.iter_mut() {
            *x = rng.gen_range(0..GOLDILOCKS_P);
        }
        let orig = data.clone();

        table.forward_inplace(&mut data);
        table.backward_inplace(&mut data);
        assert_eq!(data, orig);
    }

    #[test]
    fn goldi_test_poly_mul_n8() {
        let n = 8;
        let table = GoldilocksNttTable::with_params(n).unwrap();
        let mut rng = thread_rng();

        let mut a = vec![0u64; n];
        let mut b = vec![0u64; n];
        for i in 0..n {
            a[i] = rng.gen_range(0..GOLDILOCKS_P);
            b[i] = rng.gen_range(0..GOLDILOCKS_P);
        }

        let c_naive = poly_negacyclic_naive(&a, &b);

        let mut fa = a.clone();
        let mut fb = b.clone();
        table.forward_inplace(&mut fa);
        table.forward_inplace(&mut fb);
        for i in 0..n {
            fa[i] = mul_mod(fa[i], fb[i]);
        }
        table.backward_inplace(&mut fa);

        assert_eq!(fa, c_naive);
    }

    #[test]
    fn goldi_test_poly_mul_n4() {
        let table = GoldilocksNttTable::with_params(4).unwrap();
        let a = [1u64, 10, 20, 3];
        let b = [5u64, 6, 7, 11];

        let c_naive = poly_negacyclic_naive(&a, &b);

        let mut fa = a.to_vec();
        let mut fb = b.to_vec();
        table.forward_inplace(&mut fa);
        table.forward_inplace(&mut fb);
        for i in 0..4 {
            fa[i] = mul_mod(fa[i], fb[i]);
        }
        table.backward_inplace(&mut fa);

        assert_eq!(fa, c_naive);
    }

    #[test]
    fn goldi_test_forward_backward_defaul() {
        let table = GoldilocksNttTable::new();
        let n = table.size();
        let mut rng = thread_rng();
        let mut data = vec![0u64; n];
        for x in data.iter_mut() {
            *x = rng.gen_range(0..GOLDILOCKS_P);
        }
        let orig = data.clone();

        table.forward_inplace(&mut data);
        table.backward_inplace(&mut data);
        assert_eq!(data, orig);
    }

    #[test]
    #[ignore]
    fn goldi_test_default_polymul_2pow16() {
        let table = GoldilocksNttTable::new();
        let n = table.size();
        let mut rng = thread_rng();

        let mut a = vec![0u64; n];
        let mut b = vec![0u64; n];
        for i in 0..n {
            a[i] = rng.gen_range(0..GOLDILOCKS_P);
            b[i] = rng.gen_range(0..GOLDILOCKS_P);
        }

        let c_naive = {
            let mut c = vec![0u64; n];
            for i in 0..n {
                for j in 0..n {
                    let prod = mul_mod(a[i], b[j]);
                    let idx = i + j;
                    if idx < n {
                        c[idx] = add_mod(c[idx], prod);
                    } else {
                        c[idx - n] = sub_mod(c[idx - n], prod);
                    }
                }
            }
            c
        };

        let mut fa = a.clone();
        let mut fb = b.clone();
        table.forward_inplace(&mut fa);
        table.forward_inplace(&mut fb);
        for i in 0..n {
            fa[i] = mul_mod(fa[i], fb[i]);
        }
        table.backward_inplace(&mut fa);

        assert_eq!(fa, c_naive);
    }
}

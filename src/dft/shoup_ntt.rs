use rand::Rng;
use crate::dft::DFT;
use crate::dft::field::{mul as field_mul, inv as field_inv, exp as field_exp, add as field_add, sub as field_sub};
use crate::dft::shoup_field::{shoup_mul, shoup_precompute};

/// Shoup 方式での NTT テーブル
pub struct ShoupTable {
    /// NTT-friendly prime
    q: u64,
    /// n = power-of-two
    n: usize,
    /// 2n-th root of unity
    psi: u64,
    /// inverse of psi
    psi_inv: u64,

    /// Bit-reversed twiddle factors, powers of psi (or psi^2)
    fwd_twid: Vec<(u64,u64)>,
    /// same but inverse
    inv_twid: Vec<(u64,u64)>,

    inv_n: (u64, u64),
}

impl ShoupTable {
    pub fn new() -> Self {
        let q = 0x1fffffffffe00001u64;
        let n = 1 << 16;
        let psi = 0x15eb043c7aa2b01fu64;
        let psi_inv = field_inv(psi, q).expect("cannot invert psi");

        let inv_n_val = field_inv(n as u64, q).expect("cannot invert n");
        let inv_n_pair = shoup_precompute(inv_n_val, q);

        let (mut fwd_twid, mut inv_twid) = build_bitrev_tables(q, n, psi, psi_inv);

        Self {
            q,
            n,
            psi,
            psi_inv,
            fwd_twid,
            inv_twid,
            inv_n: inv_n_pair,
        }
    }

    /// dynamic params
    pub fn with_params(q: u64, n: usize) -> Option<Self> {
        if !n.is_power_of_two() {
            return None;
        }
        if (q - 1) % (2*n as u64) != 0 {
            return None;
        }

        // find psi
        let (psi, psi_inv) = find_primitive_2nth_root_of_unity(q, n)?;

        let inv_n_val = field_inv(n as u64, q)?;
        let inv_n_pair = shoup_precompute(inv_n_val, q);

        let (mut fwd_twid, mut inv_twid) = build_bitrev_tables(q, n, psi, psi_inv);

        Some(Self {
            q,
            n,
            psi,
            psi_inv,
            fwd_twid,
            inv_twid,
            inv_n: inv_n_pair,
        })
    }

    /// query the prime field
    pub fn q(&self) -> u64 {
        self.q
    }

    /// query the n
    pub fn size(&self) -> usize {
        self.n
    }

    pub fn forward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {
        let q = self.q;
        let mut t = self.n;
        let mut m = 1;
        while m < self.n {
            t >>= 1;
            for i in 0..m {
                let (w, w_shoup) = self.fwd_twid[m + i];
                let j1 = 2*i*t;
                let j2 = j1 + t - 1;
                for j in j1..=j2 {
                    let u = a[j];
                    let v = shoup_mul(a[j + t], (w, w_shoup), q);
                    let sum_ = field_add(u, v, q);
                    let diff_ = field_sub(u, v, q);
                    a[j] = sum_;
                    a[j + t] = diff_;
                }
            }
            m <<= 1;
        }
    }

    pub fn backward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {
        let q = self.q;
        let mut t = 1;
        let mut m = self.n;
        while m > 1 {
            let h = m >> 1;
            for i in 0..h {
                let (w, w_shoup) = self.inv_twid[h + i];
                let j1 = 2*i*t;
                let j2 = j1 + t - 1;
                for j in j1..=j2 {
                    let u = a[j];
                    let v = a[j + t];
                    let sum_ = field_add(u, v, q);
                    let diff_ = field_sub(u, v, q);
                    let diff_m = shoup_mul(diff_, (w, w_shoup), q);
                    a[j] = sum_;
                    a[j + t] = diff_m;
                }
            }
            t <<= 1;
            m = h;
        }

        // multiply by inv_n
        let (inv_n_val, inv_n_shoup) = self.inv_n;
        for x in a {
            *x = shoup_mul(*x, (inv_n_val, inv_n_shoup), q);
        }
    }
}

impl DFT<u64> for ShoupTable {
    /// NTT forward routine
    ///
    /// - `a`: vector with each element in range `[0, q)`
    fn forward_inplace(&self, a: &mut [u64]) {
        self.forward_inplace_core::<false>(a)
    }

    /// NTT forward lazy routine
    ///
    /// - `a`: vector with each element in range `[0, 2q)`
    fn forward_inplace_lazy(&self, a: &mut [u64]) {
        self.forward_inplace_core::<true>(a)
    }

    /// NTT backward routine
    ///
    /// - `a`: vector with each element in range `[0, q)`
    fn backward_inplace(&self, a: &mut [u64]) {
        self.backward_inplace_core::<false>(a)
    }

    /// NTT backward lazy routine
    ///
    /// - `a`: vector with each element in range `[0, 2q)`
    fn backward_inplace_lazy(&self, a: &mut [u64]) {
        self.backward_inplace_core::<true>(a)
    }
}

// twidの定義
fn build_bitrev_tables(q: u64, n: usize, psi: u64, psi_inv: u64) -> (Vec<(u64,u64)>, Vec<(u64,u64)>) {
    let log_n = n.trailing_zeros();
    let mut fwd = vec![(0u64,0u64); n];
    let mut inv = vec![(0u64,0u64); n];

    let mut power_psi = 1u64;
    let mut power_psi_inv = 1u64;

    for i in 0..n {
        let r = (i as u32).reverse_bits() >> (32 - log_n);
        let ridx = r as usize;

        fwd[ridx] = shoup_precompute(power_psi, q);
        inv[ridx] = shoup_precompute(power_psi_inv, q);

        power_psi = field_mul(power_psi, psi, q);
        power_psi_inv = field_mul(power_psi_inv, psi_inv, q);
    }
    (fwd, inv)
}

// dynamic paramの探索
fn find_primitive_2nth_root_of_unity(q: u64, n: usize) -> Option<(u64,u64)> {
    let mut rng = rand::thread_rng();

    // 2n = 2 * n
    let two_n = 2*(n as u64);
    if (q - 1) % two_n != 0 {
        return None;
    }
    let exponent = (q - 1) / two_n;

    loop {
        let x_random = rng.gen_range(1..q);

        let g = field_exp(x_random, exponent, q);

        // g^n == q-1" (≡ -1 mod q)
        let g_n = field_exp(g, n as u64, q);
        if g_n == q.wrapping_sub(1) {
            // "g^(2n) == 1"
            let g_2n = field_exp(g, two_n, q);
            if g_2n == 1 {
                if let Some(g_inv) = field_inv(g, q) {
                    return Some((g, g_inv));
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dft::DFT;
    use crate::dft::field::{
        add as field_add, sub as field_sub, mul as field_mul, exp as field_exp
    };
    use rand::thread_rng;
    use rand::Rng;

    #[test]
    fn test_find_primitive_2n_root_of_unity() {
        let q = 7681u64;
        let n = 16; 
        // 2n=32
        // (q-1)=7680 は 2n=32を割り切る => 7680/32=240
    
        let got = find_primitive_2nth_root_of_unity(q, n);
        assert!(got.is_some(), "cannot find");
        let (g, g_inv) = got.unwrap();
    
        let inv_check = field_mul(g, g_inv, q);
        assert_eq!(inv_check, 1, "g*g_inv != 1");

        let g_2n = field_exp(g, (2 * n) as u64, q);
        assert_eq!(g_2n, 1, "g^(2n) != 1 mod q");

        let g_n = field_exp(g, n as u64, q);
        assert_ne!(g_n, 1, "g^n == 1 mod q");       
    }

    #[test]
    fn shoup_forward_backward_small_prime() {
        let q = 7681u64;
        let n = 16usize;
        let table = ShoupTable::with_params(q, n)
            .expect("cannot build");
    
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
    fn shoup_forward_backward_default(){
        let table=ShoupTable::new();
        let q=table.q();
        let n=table.size();
        let mut rng=thread_rng();
        let mut data=vec![0u64;n];
        for x in data.iter_mut(){
            *x=rng.gen_range(0..q);
        }
        let orig=data.clone();
        table.forward_inplace(&mut data);
        table.backward_inplace(&mut data);
        assert_eq!(data,orig);
    }

    /// c(x) = a(x)*b(x) (mod x^n + 1).
    /// (i + j >= n) のとき、c[(i+j) - n] にマイナスを加える
    pub fn negacyclic_mult(a: &[u64], b: &[u64], q: u64) -> Vec<u64> {
        let n = a.len();
        let mut c = vec![0u64; n];

        for i in 0..n {
            for j in 0..n {
                // index = i + j
                let idx = i + j;
                let product = field_mul(a[i], b[j], q);

                if idx < n {
                    // c[idx] += a[i]*b[j]
                    c[idx] = field_add(c[idx], product, q);
                } else {
                    // c[idx - n] -= a[i]*b[j]
                    let wrapped_idx = idx - n;
                    c[wrapped_idx] = field_sub(c[wrapped_idx], product, q);
                }
            }
        }

        c
    }

    fn point_multiply(a:&mut[u64],b:&mut[u64],q:u64){
        for i in 0..a.len(){
            a[i]=field_mul(a[i],b[i],q);
        }
    }

    #[test]
    fn shoup_test_ntt_polymul_small_prime() {
        let q = 7681u64;
        let n = 8usize;

        let table = ShoupTable::with_params(q, n)
            .expect("error");

        let mut rng = rand::thread_rng();
        let mut a = vec![0u64; n];
        let mut b = vec![0u64; n];
        for i in 0..n {
            a[i] = rng.gen_range(0..q);
            b[i] = rng.gen_range(0..q);
        }

        let result_naive = negacyclic_mult(&a, &b, q);

        table.forward_inplace(&mut a);
        table.forward_inplace(&mut b);
        point_multiply(&mut a, &mut b, q);
        table.backward_inplace(&mut a);

        assert_eq!(a, result_naive);
    }

    #[test]
    #[ignore]
    fn shoup_test_ntt_polymul_default() {
        let table = ShoupTable::new();
        let q=table.q();
        let n=table.size();

        let mut rng = rand::thread_rng();
        let mut a = vec![0u64; n];
        let mut b = vec![0u64; n];
        for i in 0..n {
            a[i] = rng.gen_range(0..q);
            b[i] = rng.gen_range(0..q);
        }

        let result_naive = negacyclic_mult(&a, &b, q);

        table.forward_inplace(&mut a);
        table.forward_inplace(&mut b);
        point_multiply(&mut a, &mut b, q);
        table.backward_inplace(&mut a);

        assert_eq!(a, result_naive);
    }

    /// a(x) = 1 + 10x + 20x^2 + 3x^3
    /// b(x) = 5 +  6x +  7x^2 + 11x^3
    ///
    /// x^4 + 1 で割るnegacyclic乗算
    /// c = [28, 9, 37, 22] (mod 97)
    #[test]
    fn shoup_test_ntt_polymul_small_prime_evaluation() {
        let q = 97u64;
        let n = 4usize;

        let table = ShoupTable::with_params(q, n)
            .expect("cannot build");

        let mut a = vec![1u64, 10, 20, 3];
        let mut b = vec![5u64,  6,  7, 11];
        let expected = vec![28u64, 9, 37, 22];

        table.forward_inplace(&mut a);
        table.forward_inplace(&mut b);
        for i in 0..n {
            a[i] = field_mul(a[i], b[i], q);
        }
        table.backward_inplace(&mut a);

        assert_eq!(a, expected);
    }

}
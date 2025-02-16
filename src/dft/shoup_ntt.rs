use rand::Rng;
use crate::dft::DFT;
use crate::dft::field::{mul, inv, exp, add, sub};
use crate::dft::shoup_field::{shoup_mul, shoup_precompute};
use crate::dft::util::{find_primitive_2nth_root_of_unity_64};

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
    #[inline]
    pub fn new() -> Self {
        let q = 0x1fffffffffe00001u64;
        let n = 1 << 16;
        let psi = 0x15eb043c7aa2b01fu64;
        let psi_inv = inv(psi, q).expect("cannot invert psi");

        let inv_n_val = inv(n as u64, q).expect("cannot invert n");
        let inv_n_pair = shoup_precompute(inv_n_val, q);

        let (mut fwd_twid, mut inv_twid) = shoup_build_bitrev_tables(q, n, psi, psi_inv);

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
    #[inline]
    pub fn with_params(q: u64, n: usize) -> Option<Self> {
        if !n.is_power_of_two() {
            return None;
        }
        if (q - 1) % (2*n as u64) != 0 {
            return None;
        }

        // find psi
        let (psi, psi_inv) = find_primitive_2nth_root_of_unity_64(q, n)?;

        let inv_n_val = inv(n as u64, q)?;
        let inv_n_pair = shoup_precompute(inv_n_val, q);

        let (mut fwd_twid, mut inv_twid) = shoup_build_bitrev_tables(q, n, psi, psi_inv);

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
    #[inline(always)]
    pub fn q(&self) -> u64 {
        self.q
    }

    /// query the n
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
                let base = 2 * i * half; // j1
                let end = base + half;   // j2+1
  
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

        // multiply by inv_n
        for x in a.iter_mut() {
            *x = shoup_mul(*x, (inv_n_val, inv_n_shoup), q);
        }
    }
}

impl DFT<u64> for ShoupTable {
    /// NTT forward routine
    ///
    /// - `a`: vector with each element in range `[0, q)`
    #[inline(always)]
    fn forward_inplace(&self, a: &mut [u64]) {
        self.forward_inplace(a);
    }

    /// NTT backward routine
    ///
    /// - `a`: vector with each element in range `[0, q)`
    #[inline(always)]
    fn backward_inplace(&self, a: &mut [u64]) {
        self.backward_inplace(a);
    }
}

// twidの定義
#[inline]
fn shoup_build_bitrev_tables(q: u64, n: usize, psi: u64, psi_inv: u64) -> (Vec<(u64,u64)>, Vec<(u64,u64)>) {
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

        power_psi = mul(power_psi, psi, q);
        power_psi_inv = mul(power_psi_inv, psi_inv, q);
    }
    (fwd, inv)
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

}
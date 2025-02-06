use crate::dft::DFT;
use crate::dft::mont_field::{
    MontgomeryContext, mont_add, mont_sub, mont_mul, mont_inv
};
use crate::dft::field::{
    mul as field_mul, exp as field_exp, inv as field_inv
};
use rand::Rng;

pub struct Table{
    /// NTT friendly prime
    q: u64,
    /// n, the power of 2
    n : usize,
    /// 2n-th root of unity
    psi: u64,
    /// inverse of psi
    psi_inv: u64,
    
    /// Bit-reversed twiddle factors, powers of psi (or psi^2)
    pub fwd_twid: Vec<u64>,
    /// same but inverse
    pub inv_twid: Vec<u64>,

    /// n^-1 mod q in mont form for INTT
    inv_n: u64,

    mont: MontgomeryContext,
}

impl Table {
    pub fn new() -> Self {
        let q = 0x1fffffffffe00001u64;
        let n =  1<<16;
        let psi =  0x15eb043c7aa2b01fu64; //2^17th root of unity
        let psi_inv = field_inv(psi, q).expect("cannot calc invere of psi");

        let (fwd_twid, inv_twid) = build_bitrev_tables(q, n, psi, psi_inv);

        let k: u32 = 62;
        let mont = MontgomeryContext::new(q, k);
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

    /// dynamic params
    pub fn with_params(q: u64, n: usize) -> Option<Self> {
        if !n.is_power_of_two() {
            return None;
        }
        if (q - 1) % (2 * n as u64) != 0 {
            return None;
        }

        // find psi
        let (psi, psi_inv) = find_primitive_2nth_root_of_unity(q, n)?;

        let (fwd_twid, inv_twid) = build_bitrev_tables(q, n, psi, psi_inv);

        let k: u32 = 62;
        if q >= (1u64 << k) {
            return None;
        }
        let mont = MontgomeryContext::new(q, k);
        let n_mont = mont.to_mont(n as u64);
        let inv_n = mont_inv(n_mont, &mont)?;

        Some(Self {
            q, n,
            psi,
            psi_inv,
            fwd_twid,
            inv_twid,
            inv_n,
            mont,
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

    pub fn forward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {}
    pub fn backward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {}

}


impl DFT<u64> for Table {
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


impl Default for Table {
    fn default() -> Self {
        Self::new()
    }
}

// twidの定義
fn build_bitrev_tables(_q: u64, _n: usize, _psi: u64, _psi_inv: u64) -> (Vec<u64>, Vec<u64>) {
    (vec![], vec![])
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

        if n > 1 {
            let check_half = field_exp(g, (n/2) as u64, q);
            if check_half != 1 {
                let g_inv = field_inv(g, q)?;
                return Some((g, g_inv));
            }
        }
        else {
            break;
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
        assert!(got.is_some(), "should be able to find 2n-th root for (7681,16)");
        let (g, g_inv) = got.unwrap();
    
        let inv_check = field_mul(g, g_inv, q);
        assert_eq!(inv_check, 1, "g*g_inv != 1");
    }
    
}
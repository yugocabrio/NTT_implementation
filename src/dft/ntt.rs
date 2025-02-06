use crate::dft::DFT;
use crate::dft::mont_field::{
    MontgomeryContext, mont_add, mont_sub, mont_mul, mont_inv
};
use crate::dft::field::{
    mul as field_mul, exp as field_exp, inv as field_inv
};

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
        Self {
            q: 0x1fffffffffe00001u64,
            n: 0,
            psi: 0x15eb043c7aa2b01fu64, //2^17th root of unity
            psi_inv: 0,
            fwd_twid: vec![],
            inv_twid: vec![],
            inv_n: 0,
            mont: MontgomeryContext::new(1, 1),
        }
    }

    /// dynamic params
    pub fn with_params(_q: u64, _n: usize) -> Option<Self> {
        None
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
fn find_2n_root_of_unity(_q: u64, _n: usize) -> Option<(u64,u64)> {
    None
}

fn gcd(_a:u64, _b:u64)->u64 {
    1
}
fn is_generator_candidate(_g:u64,_q:u64)->bool {
    false
}
fn factorize_small(_x: &mut u64)->Vec<u64>{
    vec![]
}
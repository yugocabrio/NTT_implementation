/// MongomeryContext
#[derive(Debug)]
pub struct MongomeryContext {
    pub m: u64,
    pub r: u64,
    pub n_prime: u64,
    pub k: u32,
}

impl MongomeryContext {
    pub fn new(m: u64, k: u32) -> Self {
        // R = 2^k
        let r = 1u64 << k;
        // m < 2^k and m is odd
        assert!(m < r, "m must be < R; 2^k");
        assert_eq!(m& 1, 1, "m must be odd number, gcd(m, 2) = 1");

        // n' = -m^-1 mod R
        let m_inv_mod_r = inv_mod_u64(m, r).expect("must invert under 2^k");
        let n_prime = r.wrapping_sub(m_inv_mod_r);

        MongomeryContext {m, r, n_prime, k }
    }
}

fn inv_mod_u64(_a: u64, _m: u64) -> Option<u64> {
    None
}
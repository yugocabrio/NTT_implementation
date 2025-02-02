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

/// モンテゴメリのmの逆元を求めるため
fn inv_mod_u64(a: u64, m: u64) -> Option<u64> {
    // gcd(a, m)とx,yを拡張ユークリッドで求める
    // aが逆元を求めたい
    // mがその時のmod
    // gcdが1なら
    // xがa^-1 mod bに該当する
    let (gcd, x, _) = extended_gcd(a as i64, m as i64);

    if gcd != 1 {
        return None;
    }

    let mut x_mod_m = x % (m as i64);
    if x_mod_m < 0 {
        x_mod_m += m as i64;
    }
    Some(x_mod_m as u64)
}

/// gcd(a, b) = a*s0 + b*t0, 最終的に (gcd, si, ti) を返す
fn extended_gcd(mut r0: i64, mut r1: i64) -> (i64, i64, i64) {
    let (mut s0, mut s1) = (1, 0);
    let (mut t0, mut t1) = (0, 1);

    while r1 != 0 {
        let q = r0 / r1;

        let next_r = r0 - q * r1;
        r0 = r1;
        r1 = next_r;

        let next_s = s0 - q * s1;
        s0 = s1;
        s1 = next_s;

        let next_t = t0 - q * t1;
        t0 = t1;
        t1 = next_t;
    }

    (r0, s0, t0)
}
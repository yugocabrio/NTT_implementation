/// MongomeryContext
#[derive(Debug)]
pub struct MontgomeryContext {
    pub m: u64,
    pub r: u64,
    pub n_prime: u64,
    pub k: u32,
}

impl MontgomeryContext {
    pub fn new(m: u64, k: u32) -> Self {
        // R = 2^k
        let r = 1u64 << k;
        // m < 2^k and m is odd
        assert!(m < r, "m must be < R; 2^k");
        assert_eq!(m& 1, 1, "m must be odd number, gcd(m, 2) = 1");

        // n' = -m^-1 mod R
        let m_inv_mod_r = inv_mod_u64(m, r).expect("must invert under 2^k");
        let n_prime = r.wrapping_sub(m_inv_mod_r);

        MontgomeryContext {m, r, n_prime, k }
    }

    pub fn mont_reduce(&self, t: u128) -> u64 {
        // R=2^k, mask = R -1 = 2^k - 1は、下位kビットが全部1
        let mask = self.r.wrapping_sub(1);
        // T mod R
        let t_mod_r = (t as u64) & mask;

        let u = (t_mod_r as u128).wrapping_mul(self.n_prime as u128) & (mask as u128);
        // tmpがrの倍数
        let tmp = t.wrapping_add(u.wrapping_mul(self.m as u128));

        // (tmp / 2^k)
        // Rの1回の割り算に相当
        let big = tmp >> self.k;
        let mut res = big as u64;
        if res >= self.m {
            res -= self.m;
        }
        res
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_montgomery_context_new_basic() {
        // m=17, k=5 => r=32
        let ctx = MontgomeryContext::new(17, 5);
        assert_eq!(ctx.m, 17);
        assert_eq!(ctx.r, 32);

        let expected_m_inv = inv_mod_u64(17, 32).unwrap();
        let expected_n_prime = 32u64.wrapping_sub(expected_m_inv);
        assert_eq!(ctx.n_prime, expected_n_prime);
    }

    #[test]
    fn test_extended_gcd() {
        // gcd = gcd(a,b)
        // g = a*x + b*y
        let (gcd, x, y) = extended_gcd(30, 18);
        assert_eq!(gcd, 6);
        assert_eq!(30*x + 18*y, 6);
    }

    #[test]
    fn test_inv_mod_u64() {
        // a=5, m=17 => 5^-1 mod 17=7 (5*7=35≡1 mod 17)
        let inv = inv_mod_u64(5, 17).unwrap();
        assert_eq!((5 * inv) % 17, 1);

        // gcd(6,12)=6=>逆元なし
        assert_eq!(inv_mod_u64(6, 12), None);
    }
    
}
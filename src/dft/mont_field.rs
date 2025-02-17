/// Holds parameters for Montgomery reduction.
/// m: the modulus
/// r: R = 2^k
/// n_prime: precomputed factor for montgomery reduction
/// k: the exponent for R = 2^k
/// mask: bitmask for mod R operations (2^k - 1)
#[derive(Debug)]
pub struct MontgomeryContext {
    pub m: u64,
    pub r: u64,
    pub n_prime: u64,
    pub k: u32,
    pub mask: u64,
}

impl MontgomeryContext {
    /// Constructs a MontgomeryContext for a prime m < 2^k (odd) and the given k.
    /// This precomputes n_prime which is -m^-1 mod R, used for montgomery reduction.
    #[inline]
    pub fn new(m: u64, k: u32) -> Self {
        // R = 2^k
        let r = 1u64 << k;
        // m < 2^k and m is odd
        assert!(m < r, "m must be < R; 2^k");
        assert_eq!(m & 1, 1, "m must be odd number, gcd(m, 2) = 1");

        // n' = -m^-1 mod R
        let m_inv_mod_r = inv_mod_u64(m, r).expect("must invert under 2^k");
        let n_prime = r.wrapping_sub(m_inv_mod_r);
        // mask = R - 1
        let mask = r.wrapping_sub(1);

        MontgomeryContext {
            m,
            r,
            n_prime,
            k,
            mask,
        }
    }

    /// Performs the Montgomery reduction of a 128-bit value t into a 64-bit number mod m.
    /// Internally uses the mask to emulate mod R=2^k and then subtracts m if necessary.
    #[inline(always)]
    pub fn mont_reduce(&self, t: u128) -> u64 {
        // R=2^k, mask = R -1 = 2^k - 1, the lower k bits are all 1
        let t_mod_r = (t as u64) & self.mask;

        let u = (t_mod_r as u128).wrapping_mul(self.n_prime as u128) & (self.mask as u128);
        // tmp is multiple of r
        let tmp = t.wrapping_add(u.wrapping_mul(self.m as u128));

        // tmp divide by 2^k
        let big = tmp >> self.k;
        let mut res = big as u64;
        if res >= self.m {
            res = res.wrapping_sub(self.m);
        }
        res
    }

    /// Converts a normal integer x to its Montgomery form (x * R) mod m.
    #[inline(always)]
    pub fn to_mont(&self, x: u64) -> u64 {
        let xr = (x as u128).wrapping_mul(self.r as u128);
        (xr % (self.m as u128)) as u64
    }

    /// Converts a number in Montgomery form x back to normal form (x / R) mod m.
    #[inline(always)]
    pub fn from_mont(&self, x: u64) -> u64 {
        self.mont_reduce(x as u128)
    }
}

/// Multiplies two numbers in Montgomery form and reduces them
/// x_mont * y_mont)/R mod m)
#[inline(always)]
pub fn mont_mul(x_mont: u64, y_mont: u64, ctx: &MontgomeryContext) -> u64 {
    let t = (x_mont as u128) * (y_mont as u128);
    ctx.mont_reduce(t)
}

/// Adds two numbers in Montgomery form (mod m).
#[inline(always)]
pub fn mont_add(x_mont: u64, y_mont: u64, ctx: &MontgomeryContext) -> u64 {
    let s = x_mont.wrapping_add(y_mont);
    if s >= ctx.m {
        s.wrapping_sub(ctx.m)
    } else {
        s
    }
}

/// Subtracts two numbers in Montgomery form (mod m).
#[inline(always)]
pub fn mont_sub(x_mont: u64, y_mont: u64, ctx: &MontgomeryContext) -> u64 {
    let d = x_mont.wrapping_sub(y_mont);
    if x_mont < y_mont {
        d.wrapping_add(ctx.m)
    } else {
        d
    }
}

/// Exponentiates a Montgomery-form base base_mont by exp (mod m).
#[inline(always)]
pub fn mont_exp(base_mont: u64, exp: u64, ctx: &MontgomeryContext) -> u64 {
    let mut result = ctx.to_mont(1);
    let mut cur = base_mont;
    let mut e = exp;
    while e > 0 {
        if (e & 1) == 1 {
            result = mont_mul(result, cur, ctx);
        }
        cur = mont_mul(cur, cur, ctx);
        e >>= 1;
    }
    result
}

/// Computes the inverse of x_mont (in Montgomery form) using Fermat's little theorem.
#[inline]
pub fn mont_inv(x_mont: u64, ctx: &MontgomeryContext) -> Option<u64> {
    if x_mont == 0 {
        return None;
    }
    let exp = ctx.m.wrapping_sub(2);
    let result_mont = mont_exp(x_mont, exp, ctx);
    Some(result_mont)
}

/// Finds the inverse of a modulo m (here m is 2^k) by extended Euclid. 
/// This is used for Montgomery setup.
#[inline]
fn inv_mod_u64(a: u64, m: u64) -> Option<u64> {
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

/// Extended GCD on (r0, r1); returns (gcd, s0, t0) with gcd = r0*s0 + r1*t0.
/// Used for computing modular inverses.
#[inline]
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
    fn test_to_mont_from_mont() {
        let ctx = MontgomeryContext::new(17, 5);
        // 10*32 mod 17 = 14
        let x_mont = ctx.to_mont(10); 
        assert_eq!(x_mont, 14);

        // 10 mod 17
        let back = ctx.from_mont(x_mont);
        assert_eq!(back, 10);
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

        // gcd(6,12)=6=>None
        assert_eq!(inv_mod_u64(6, 12), None);
    }
    
    #[test]
    fn test_mont_mul() {
        // m = 17, k = 5, R = 32
        let ctx = MontgomeryContext::new(17, 5);
        let x = 3;
        let y = 5;
        let x_mont = ctx.to_mont(x); 
        let y_mont = ctx.to_mont(y);

        // (3*5)%17=15
        let z_mont = mont_mul(x_mont, y_mont, &ctx);
        let z = ctx.from_mont(z_mont);

        assert_eq!(z, (x * y) % 17);
    }

    #[test]
    fn test_mont_add() {
        let ctx = MontgomeryContext::new(17, 5);

        let x = 12;
        let y = 4;
        let x_mont = ctx.to_mont(x);
        let y_mont = ctx.to_mont(y);

        let z_mont = mont_add(x_mont, y_mont, &ctx);
        let z = ctx.from_mont(z_mont);

        assert_eq!(z, (x + y) % 17);
    }

    #[test]
    fn test_mont_sub() {
        let ctx = MontgomeryContext::new(17, 5);

        let x = 12;
        let y = 4;
        let x_mont = ctx.to_mont(x);
        let y_mont = ctx.to_mont(y);

        let z_mont = mont_sub(x_mont, y_mont, &ctx);
        let z = ctx.from_mont(z_mont);

        // (12 - 4) mod 17 = 8
        assert_eq!(z, (x + 17 - y) % 17);
    }

    #[test]
    fn test_mont_exp() {
        let ctx = MontgomeryContext::new(17, 5);

        let base = 3;
        let exp = 4;

        // 3^4=81 = 81%17=13
        let base_mont = ctx.to_mont(base);
        let result_mont = mont_exp(base_mont, exp, &ctx);
        let result = ctx.from_mont(result_mont);

        assert_eq!(result, 13);
    }

    #[test]
    fn test_mont_inv() {
        let ctx = MontgomeryContext::new(17, 5);

        // inverse of 5 is 7
        let x_mont = ctx.to_mont(5);
        let result_mont = mont_inv(x_mont, &ctx).expect("must have inverse");
        let result = ctx.from_mont(result_mont);

        assert_eq!((5 * result) % 17, 1);

        // inverse of 0
        let x_mont = ctx.to_mont(0);
        assert_eq!(mont_inv(x_mont, &ctx), None);
    }
}
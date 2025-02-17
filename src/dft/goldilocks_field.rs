/// Goldilocks prime: 2^64 - 2^32 + 1
pub const GOLDILOCKS_P: u64 = 0xFFFF_FFFF_0000_0001;

/// Performs a special 128→64 bit reduction using the PZT22 method.
/// For Goldilocks Prime, a single subtraction is required to obtain the correct remainder.
/// Interpret x = x3*(2^32)^3 + x2*(2^32)^2 + x1*(2^32) + x0
/// Then apply some manipulations (t0, t1), adding or subtracting (2^32 -1).
pub fn reduce128(x: u128) -> u64 {
    const EPSILON: u64 = (1 << 32) - 1;

    let x_lo = x as u64; // lower 64 bits (x0..x1)
    let x_hi = (x >> 64) as u64; // upper 64 bits (x2..x3)

    // break x_hi into x3 (upper 32 bits) and x2 (lower 32 bits)
    let x3 = x_hi >> 32;
    let x2 = x_hi & EPSILON;

    // t0 = x0 - x3, with possible borrow => subtract (2^32 - 1)
    let (mut t0, borrow) = x_lo.overflowing_sub(x3);
    if borrow {
        t0 = t0.wrapping_sub(EPSILON);
    }

    // t1 = x2 * (2^32 - 1)
    let t1 = x2.wrapping_mul(EPSILON);

    // add them with possible carry => also add (2^32 -1)
    let (res_wrapped, carry) = t0.overflowing_add(t1);
    let mut r = res_wrapped.wrapping_add(EPSILON * (carry as u64));

    // single subtraction by GOLDILOCKS_P if needed
    if r >= GOLDILOCKS_P {
        r -= GOLDILOCKS_P;
    }

    r
}

/// (a*b mod GOLDILOCKS_P) using reduce128
#[inline(always)]
pub fn mul(a: u64, b: u64) -> u64 {
    reduce128((a as u128) * (b as u128))
}

/// (a+b mod GOLDILOCKS_P)
#[inline(always)]
pub fn add(a: u64, b: u64) -> u64 {
    let (res, carry) = a.overflowing_add(b);
    let mut sum = res;
    if carry || sum >= GOLDILOCKS_P {
        sum = sum.wrapping_sub(GOLDILOCKS_P);
    }
    sum
}

/// (a-b mod GOLDILOCKS_P)
pub fn sub(a: u64, b: u64) -> u64 {
    let (diff, borrow) = a.overflowing_sub(b);
    if borrow {
        diff.wrapping_add(GOLDILOCKS_P)
    } else {
        diff
    }
}

/// (base^exp mod GOLDILOCKS_P)
#[inline(always)]
pub fn exp(mut base: u64, mut exp: u64) -> u64 {
    let mut result = 1u64;
    while exp > 0 {
        if (exp & 1) != 0 {
            result = mul(result, base);
        }
        base = mul(base, base);
        exp >>= 1;
    }
    result
}

/// (a^(p-2) mod GOLDILOCKS_P)
#[inline(always)]
pub fn inv(a: u64) -> u64 {
    exp(a, GOLDILOCKS_P - 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce128() {
        // (p-1)*(p-11) = p^2 -2p +1
        let a = (GOLDILOCKS_P - 1) as u128;
        let prod = a * a;
        let r = reduce128(prod);
        assert_eq!(r, 1);
    }

    #[test]
    fn test_add_mod() {
        let a = GOLDILOCKS_P - 13;
        let b = 13u64;
        let s = add(a, b);
        assert_eq!(s, 0);
    }

    #[test]
    fn test_sub_mod() {
        // 0 - 1 = p-1
        let a = 0u64;
        let b = 111u64;
        let d = sub(a, b);
        assert_eq!(d, GOLDILOCKS_P - 111);
    }

    #[test]
    fn test_mul_mod() {
        // (p-1)*(p-1) = p^2 -2p +1 ≡ 1 (mod p)
        let a = GOLDILOCKS_P - 1;
        let m = mul(a, a);
        assert_eq!(m, 1);
    }

    #[test]
    fn test_mod_inv() {
        // a*(a^-1) ≡ 1
        let a = 12345;
        let inv_a = inv(a);
        let chk = mul(a, inv_a);
        assert_eq!(chk, 1);
    }
}

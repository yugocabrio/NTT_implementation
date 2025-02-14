pub const GOLDILOCKS_P: u64 = 0xFFFF_FFFF_0000_0001;

/// PZT22方式
/// x = x3 * (2^32)^3 + x2 * (2^32)^2 + x1 * 2^32 + x0 とみなし、
///   t0 = x0 - x3  (borrowが発生した場合、 t0 -= (2^32 - 1))
///   t1 = x2 * (2^32 - 1)
///   r = t0 + t1 (+ carry分 * (2^32 - 1)) を p で一度だけ減算
#[inline(always)]
pub fn reduce128(x: u128) -> u64 {
    let lo = x as u64;
    let hi = (x >> 64) as u64;

    let a = hi >> 32;
    let b = hi & 0xFFFF_FFFF;
    let c = lo;

    // 2^96 = -1
    let mut t = (c as i128) - (a as i128);
    if t < 0 {
        t += GOLDILOCKS_P as i128;
    }
    let r = t as u64;

    // 2^64 = 2^32 - 1
    let add_part = ((b as u128) << 32).wrapping_sub(b as u128);
    let big = (r as u128).wrapping_add(add_part);

    let mut r2 = big as u64;
    let carry = big >> 64;
    if carry != 0 {
        r2 = r2.wrapping_sub(GOLDILOCKS_P);
    }
    if r2 >= GOLDILOCKS_P {
        r2 -= GOLDILOCKS_P;
    }
    r2
}

/// (a*b mod p)
#[inline(always)]
pub fn mul_mod(a: u64, b: u64) -> u64 {
    reduce128((a as u128) * (b as u128))
}

/// (a+b mod p)
#[inline(always)]
pub fn add_mod(a: u64, b: u64) -> u64 {
    let (res, carry) = a.overflowing_add(b);
    let mut sum = res;
    if carry || sum >= GOLDILOCKS_P {
        sum = sum.wrapping_sub(GOLDILOCKS_P);
    }
    sum
}

/// (a-b mod p)
#[inline(always)]
pub fn sub_mod(a: u64, b: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a.wrapping_sub(b).wrapping_add(GOLDILOCKS_P)
    }
}

/// (base^exp mod p)
pub fn pow_mod(mut base: u64, mut exp: u64) -> u64 {
    let mut result = 1u64;
    while exp > 0 {
        if (exp & 1) != 0 {
            result = mul_mod(result, base);
        }
        base = mul_mod(base, base);
        exp >>= 1;
    }
    result
}

/// (a^(p-2) mod p)
pub fn mod_inv(a: u64) -> u64 {
    pow_mod(a, GOLDILOCKS_P - 2)
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
        let s = add_mod(a, b);
        assert_eq!(s, 0);
    }

    #[test]
    fn test_sub_mod() {
        // 0 - 1 = p-1
        let a = 0u64;
        let b = 111u64;
        let d = sub_mod(a, b);
        assert_eq!(d, GOLDILOCKS_P - 111);
    }

    #[test]
    fn test_mul_mod() {
        // (p-1)*(p-1) = p^2 -2p +1 ≡ 1 (mod p)
        let a = GOLDILOCKS_P - 1;
        let m = mul_mod(a, a);
        assert_eq!(m, 1);
    }

    #[test]
    fn test_mod_inv() {
        // a*(a^-1) ≡ 1
        let a = 12345;
        let inv_a = mod_inv(a);
        let chk = mul_mod(a, inv_a);
        assert_eq!(chk, 1);
    }
}

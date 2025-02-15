/// (a + b) mod p
#[inline(always)]
pub fn add(a: u64, b: u64, p: u64) -> u64 {
    let (res, carry) = a.overflowing_add(b);
    let mut s = res;
    if carry || s >= p {
        s = s.wrapping_sub(p);
    }
    s
}

/// (a - b) mod p
#[inline(always)]
pub fn sub(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a.wrapping_sub(b).wrapping_add(p)
    }
}

/// (a * b) mod p
#[inline(always)]
pub fn mul(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128) * (b as u128) % (p as u128)) as u64
}

/// a^e mod p
#[inline(always)]
pub fn exp(mut base: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
    base = base % p;
    while e > 0 {
        if (e & 1) != 0 {
            r = mul(r, base, p);
        }
        base = mul(base, base, p);
        e >>= 1;
    }
    r
}

/// a^(p-2) mod p, Fermat
#[inline(always)]
pub fn inv(a: u64, p: u64) -> Option<u64> {
    if a == 0 {
        None
    } else {
        Some(exp(a, p - 2, p))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let q = 5u64;
        assert_eq!(add(4, 1, q), 0);
    }

    #[test]
    fn test_sub() {
        let q = 5u64;
        assert_eq!(sub(2, 3, q), 4);
    }

    #[test]
    fn test_mul() {
        let q = 5u64;
        assert_eq!(mul(2, 3, q), 1);
    }

    #[test]
    fn test_exp() {
        let q = 5u64;
        assert_eq!(exp(2, 4, q), 1);
    }

    #[test]
    fn test_inv() {
        let q = 5u64;
        assert_eq!(inv(4, q), Some(4));
        assert_eq!(inv(0, q), None);
    }
}

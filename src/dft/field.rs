/// a+b mod q
#[inline(always)]
pub fn add(a: u64, b: u64, p: u64) -> u64 {
    let c = a.wrapping_add(b);
    if c >= p { c - p } else { c }
}

/// a-b mod q
#[inline(always)]
pub fn sub(a: u64, b: u64, p: u64) -> u64 {
    // a + q - b
    let c = a.wrapping_add(p).wrapping_sub(b);
    if c >= p { c - p } else { c }
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
        let p= 5u64;
        assert_eq!(add(4, 1, p), 0);
        assert_eq!(add(p-1,1,p), 0);
        assert_eq!(add(0,3,p), 3);
    }

    #[test]
    fn test_sub() {
        let p= 5u64;
        assert_eq!(sub(2, 3, p), 4);
        assert_eq!(sub(p-1,p-1,p),0);
        assert_eq!(sub(0,0,p),0);
    }

    #[test]
    fn test_mul() {
        let p= 5u64;
        assert_eq!(mul(2,3,p), 1);
        assert_eq!(mul(p-1,p-1,p), 1);
        assert_eq!(mul(0,3,p),0);
    }

    #[test]
    fn test_exp() {
        let p= 5u64;
        assert_eq!(exp(2,4,p), 1);
        assert_eq!(exp(3,3,p), 2);
    }

    #[test]
    fn test_inv() {
        let p= 5u64;
        assert_eq!(inv(4,p), Some(4));
        assert_eq!(inv(0,p), None);
        assert_eq!(inv(1,p), Some(1));
    }
}
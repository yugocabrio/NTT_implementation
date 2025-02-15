/// p < 2^31 前提で "floor(2^64 / p)" を計算
#[inline(always)]
pub fn barrett_precompute(p: u32) -> u64 {
    let one = 1u128 << 64;
    (one / (p as u128)) as u64
}

/// Barrett mul(32bit版)
///   z = a*b
///   q = floor( (z * p_bar) >> 64 )
///   r = z - q*p
///   if r>=p => r-p
#[inline(always)]
pub fn barrett_mul(a: u32, (b,_): (u32,u32), p: u32, p_bar: u64) -> u32 {
    // (a*b) in 64bit
    let z = (a as u64).wrapping_mul(b as u64);

    // q = floor( (z * p_bar) >> 64 )
    let q = ((z as u128).wrapping_mul(p_bar as u128) >> 64) as u64;

    // r = z - q*p
    let r = z.wrapping_sub(q.wrapping_mul(p as u64)) as u32;

    // if r>=p => r-p
    if r >= p { r - p } else { r }
}

/// (a + b) mod p
#[inline(always)]
pub fn add(a: u32, b: u32, p: u32) -> u32 {
    let (res, carry) = a.overflowing_add(b);
    let mut s = res;
    if carry || s >= p {
        s = s.wrapping_sub(p);
    }
    s
}

/// (a - b) mod p
#[inline(always)]
pub fn sub(a: u32, b: u32, p: u32) -> u32 {
    if a >= b { a - b }
    else { a.wrapping_sub(b).wrapping_add(p) }
}

/// a*b mod p
#[inline(always)]
pub fn mul(a:u32,b:u32,p:u32)->u32 {
    ((a as u64)*(b as u64) % (p as u64)) as u32
}

/// a^e mod p
#[inline(always)]
pub fn exp(mut base: u32, mut e: u32, p: u32) -> u32 {
    let mut r = 1u32;
    base %= p; 
    while e > 0 {
        if (e & 1) != 0 {
            r = mul(r, base, p);
        }
        base = mul(base, base, p);
        e >>= 1;
    }
    r
}

/// a^(p-2) mod p
#[inline(always)]
pub fn inv(a: u32, p: u32) -> Option<u32> {
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
        let p = 5u32;
        assert_eq!(add(4, 1, p), 0);
        assert_eq!(add(p-1, 1, p), 0);
        assert_eq!(add(0,3,p), 3);
    }

    #[test]
    fn test_sub() {
        let p = 5u32;
        assert_eq!(sub(2, 3, p), 4);
        assert_eq!(sub(p-1, p-1, p), 0);
        assert_eq!(sub(0,0,p), 0);
    }

    #[test]
    fn test_mul() {
        let p = 5u32;
        assert_eq!(mul(2,3,p), 1);
        assert_eq!(mul(p-1,p-1,p), 1);
        assert_eq!(mul(0,3,p),0);
    }

    #[test]
    fn test_exp() {
        let p = 5u32;
        assert_eq!(exp(2,4,p), 1);
        assert_eq!(exp(3,3,p), 2);
    }

    #[test]
    fn test_inv() {
        let p= 5u32;
        assert_eq!(inv(4, p), Some(4));
        assert_eq!(inv(0, p), None);
        assert_eq!(inv(1, p), Some(1));
    }

    #[test]
    fn test_barrett_mul() {
        let p=13u32;
        let p_bar= barrett_precompute(p);
        let b_tuple= (4u32,0);
        let r0= barrett_mul(5, b_tuple, p, p_bar);
        assert_eq!(r0 as u64, 7);

        let r1 = barrett_mul(p-1, (p-1,0), p, p_bar);
        assert_eq!(r1 as u64, 1);

        let r2= barrett_mul(0, (4,0), p, p_bar);
        assert_eq!(r2, 0);
    }
}
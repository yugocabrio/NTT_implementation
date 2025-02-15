/// b_shoup = floor((w << 64) / p) とし、(w, b_shoup)
#[inline(always)]
pub fn shoup_precompute(w: u64, p: u64) -> (u64, u64) {
    let w_shoup = ((w as u128) << 64) / (p as u128);
    (w, w_shoup as u64)
}

/// Shoup 乗算: a * b mod p.
/// b は (b, b_shoup)
#[inline(always)]
pub fn shoup_mul(a: u64, (b, b_shoup): (u64, u64), p: u64) -> u64 {
    // hi = floor((a * b_shoup) / 2^64)
    let hi = (((a as u128) * (b_shoup as u128)) >> 64) as u64;

    // (a*b) - hi*p
    let prod_128 = (a as u128) * (b as u128);
    let sub_128  = (hi as u128) * (p as u128);
    let tmp = prod_128.wrapping_sub(sub_128);

    let ret = tmp as u64;
    if ret >= p { ret - p } else { ret }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shoup_mul() {
        let p= 13u64;
        let b_tuple= shoup_precompute(5, p);
        let shoup_res= shoup_mul(4, b_tuple, p);
        assert_eq!(shoup_res, 7);

    
        let (bw, bw_sh)= shoup_precompute(p-1, p);
        let r= shoup_mul(p-1, (bw,bw_sh), p);
        assert_eq!(r, 1);
        
        let (bw,bw_sh)= shoup_precompute(3, p);
        let r0= shoup_mul(0, (bw,bw_sh), p);
        assert_eq!(r0, 0);

        let (b0,b0_sh)= shoup_precompute(0, p);
        let r1= shoup_mul(7, (b0,b0_sh), p);
        assert_eq!(r1, 0);
    }
}
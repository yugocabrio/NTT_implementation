
/// b_shoup = floor((w << 64) / q) とし、(w, b_shoup)
#[inline(always)]
pub fn shoup_precompute(w: u64, q: u64) -> (u64, u64) {
    let w_shoup = ((w as u128) << 64) / (q as u128);
    (w, w_shoup as u64)
}

/// Shoup 乗算: a * b mod q.
/// b は (b, b_shoup)
#[inline(always)]
pub fn shoup_mul(a: u64, (b, b_shoup): (u64, u64), q: u64) -> u64 {
    // hi = floor((a * b_shoup) / 2^64)
    let hi = (((a as u128) * (b_shoup as u128)) >> 64) as u64;

    // (a*b) - hi*q
    let prod_128 = (a as u128) * (b as u128);
    let sub_128  = (hi as u128) * (q as u128);
    let tmp = prod_128.wrapping_sub(sub_128);

    let ret = tmp as u64;
    if ret >= q { ret - q } else { ret }
}
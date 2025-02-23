#[cfg(target_arch = "aarch64")]
use core::arch::aarch64::{
    uint32x4_t, vaddq_u32, vandq_u32, vcgeq_u32, vdupq_n_u32, vld1q_u32, vmvnq_u32, vst1q_u32,
    vsubq_u32, vget_low_u32, vget_high_u32, vgetq_lane_u64, vdup_n_u32, vmull_u32
};

#[cfg(target_arch = "aarch64")]
use core::mem::transmute;

/// Precomputes p_bar = floor((1 << 64) / p) for Barrett reduction.
/// This helps accelerate (a*b) mod p.
#[inline(always)]
pub fn barrett_precompute(p: u32) -> u64 {
    let one = 1u128 << 64;
    (one / (p as u128)) as u64
}

/// Barrett multiplication.
/// z = a*b (as 64-bit)
/// q = floor((z * p_bar) >> 64)
/// r = z - q*p
/// if r >= p => r -= p
#[inline(always)]
pub fn barrett_mul(a: u32, b: u32, p: u32, p_bar: u64) -> u32 {
    let z = (a as u64).wrapping_mul(b as u64);
    let q = ((z as u128).wrapping_mul(p_bar as u128) >> 64) as u64;
    let r = z.wrapping_sub(q.wrapping_mul(p as u64)) as u32;
    if r >= p {
        r - p
    } else {
        r
    }
}

/// (a + b) mod p
#[inline(always)]
pub fn add(a: u32, b: u32, p: u32) -> u32 {
    let s = a.wrapping_add(b);
    if s >= p {
        s - p
    } else {
        s
    }
}

/// (a - b) mod p
#[inline(always)]
pub fn sub(a: u32, b: u32, p: u32) -> u32 {
    let s = a.wrapping_add(p).wrapping_sub(b);
    if s >= p {
        s - p
    } else {
        s
    }
}

/// (a * b) mod p
#[inline(always)]
pub fn mul(a: u32, b: u32, p: u32) -> u32 {
    ((a as u64) * (b as u64) % (p as u64)) as u32
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
        return None;
    }
    Some(exp(a, p.wrapping_sub(2), p))
}

/// Vectorized addition using NEON (lane-wise).
/// If (u+v) >= p in a given lane, subtract p from that lane.
#[cfg(target_arch = "aarch64")]
#[inline(always)]
pub unsafe fn vec_add(u: uint32x4_t, v: uint32x4_t, p_vec: uint32x4_t) -> uint32x4_t {
    // sum = u + v (lane-wise)
    let sum = vaddq_u32(u, v);
    // cmp = (sum >= p_vec)? => 0xFFFF_FFFF if true, else 0
    let cmp = vcgeq_u32(sum, p_vec);
    // sub_p = cmp & p_vec => pick 'p' for lanes where sum >= p
    let sub_p = vandq_u32(cmp, p_vec);
    // final result: sum - p on lanes that overflowed, otherwise sum
    vsubq_u32(sum, sub_p)
}

/// Vectorized subtraction using NEON (lane-wise).
/// If (u < v) in a given lane, add p to that lane for wrap-around.
#[cfg(target_arch = "aarch64")]
#[inline(always)]
pub unsafe fn vec_sub(u: uint32x4_t, v: uint32x4_t, p_vec: uint32x4_t) -> uint32x4_t {
    // diff = u - v (lane-wise)
    let diff = vsubq_u32(u, v);
    // cmp = (u >= v)? => 0xFFFF_FFFF if true, else 0
    let cmp = vcgeq_u32(u, v);
    // vmvnq_u32(cmp) flips bits, lanes where (u < v) become 0xFFFF_FFFF
    let underflow_mask = vmvnq_u32(cmp);
    // add_p = underflow_mask & p_vec, pick p for lanes that underflowed
    let add_p = vandq_u32(underflow_mask, p_vec);
    // final result: diff + p on lanes that underflowed, otherwise diff
    vaddq_u32(diff, add_p)
}

// reduce for  max 64 bit z
#[cfg(target_arch="aarch64")]
#[inline(always)]
fn reduce_barrett_64(z: u64, p: u32, p_bar: u64) -> u32 {
    let q = ((z as u128).wrapping_mul(p_bar as u128) >> 64) as u64;
    let r = z.wrapping_sub(q.wrapping_mul(p as u64)) as u32;
    if r >= p { r - p } else { r }
}

/// (u[i] * w) mod p
#[cfg(target_arch="aarch64")]
#[inline(always)]
pub unsafe fn vec_barrett_mul_scalar(u: uint32x4_t, w: u32, p: u32, p_bar: u64) -> uint32x4_t {

    let w_vec = vdup_n_u32(w);

    let u_low  = vget_low_u32(u);
    let u_high = vget_high_u32(u);

    let lo_prod = vmull_u32(u_low, w_vec);
    let hi_prod = vmull_u32(u_high, w_vec);

    let mut out_arr = [0u32;4];
    out_arr[0] = reduce_barrett_64(vgetq_lane_u64::<0>(lo_prod), p, p_bar);
    out_arr[1] = reduce_barrett_64(vgetq_lane_u64::<1>(lo_prod), p, p_bar);
    out_arr[2] = reduce_barrett_64(vgetq_lane_u64::<0>(hi_prod), p, p_bar);
    out_arr[3] = reduce_barrett_64(vgetq_lane_u64::<1>(hi_prod), p, p_bar);

    transmute(out_arr)
}

// slow
/*
#[inline(always)]
pub fn add(a: u32, b: u32, p: u32) -> u32 {
    let (res, carry) = a.overflowing_add(b);
    let mut s = res;
    if carry || s >= p {
        s = s.wrapping_sub(p);
    }
    s
}

#[inline(always)]
pub fn sub(a: u32, b: u32, p: u32) -> u32 {
    let (res, borrow) = a.overflowing_sub(b);
    let mut s = res;
    if borrow {
        s = s.wrapping_add(p);
    }
    s
}
*/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let p = 5u32;
        assert_eq!(add(4, 1, p), 0);
        assert_eq!(add(p - 1, 1, p), 0);
        assert_eq!(add(0, 3, p), 3);
    }

    #[test]
    fn test_sub() {
        let p = 5u32;
        assert_eq!(sub(2, 3, p), 4);
        assert_eq!(sub(p - 1, p - 1, p), 0);
        assert_eq!(sub(0, 0, p), 0);
    }

    #[test]
    fn test_mul() {
        let p = 5u32;
        assert_eq!(mul(2, 3, p), 1);
        assert_eq!(mul(p - 1, p - 1, p), 1);
        assert_eq!(mul(0, 3, p), 0);
    }

    #[test]
    fn test_exp() {
        let p = 5u32;
        assert_eq!(exp(2, 4, p), 1);
        assert_eq!(exp(3, 3, p), 2);
    }

    #[test]
    fn test_inv() {
        let p = 5u32;
        assert_eq!(inv(4, p), Some(4));
        assert_eq!(inv(0, p), None);
        assert_eq!(inv(1, p), Some(1));
    }

    #[test]
    fn test_barrett_mul_32() {
        let p = 13u32;
        let p_bar = barrett_precompute(p);

        let r = barrett_mul(7, 2, p, p_bar);
        assert_eq!(r, 1);

        let r2 = barrett_mul(12, 12, p, p_bar);
        assert_eq!(r2, 1);
    }

    #[test]
    fn test_vadd_mod() {
        #[cfg(target_arch = "aarch64")]
        unsafe {
            let p = 13u32;
            let p_vec = vdupq_n_u32(p);

            let a_arr = [12u32, 0u32, 5u32, 13u32];
            let b_arr = [1u32, 13u32, 8u32, 1u32];

            let a_vec = vld1q_u32(a_arr.as_ptr());
            let b_vec = vld1q_u32(b_arr.as_ptr());

            let sum_vec = vec_add(a_vec, b_vec, p_vec);
            let sum_arr: [u32; 4] = transmute(sum_vec);

            assert_eq!(sum_arr, [0, 0, 0, 1]);
        }
    }

    #[test]
    fn test_vsub_mod() {
        #[cfg(target_arch = "aarch64")]
        unsafe {
            let p = 13u32;
            let p_vec = vdupq_n_u32(p);

            let a_arr = [0u32, 12u32, 5u32, 13u32];
            let b_arr = [1u32, 0u32, 8u32, 1u32];

            let a_vec = vld1q_u32(a_arr.as_ptr());
            let b_vec = vld1q_u32(b_arr.as_ptr());

            let diff_vec = vec_sub(a_vec, b_vec, p_vec);
            let diff_arr: [u32; 4] = transmute(diff_vec);

            assert_eq!(diff_arr, [12, 12, 10, 12]);
        }
    }

    #[test]
    fn test_vec_barrett_mul() {
        #[cfg(target_arch="aarch64")]
        unsafe {
            let p = 17u32;
            let p_bar = barrett_precompute(p);

            let a_arr = [1u32, 8u32, 16u32, 10u32];
            let w = 9u32; 
            let a_vec = transmute::<[u32;4], uint32x4_t>(a_arr);

            let w_vec = vec_barrett_mul_scalar(a_vec, w, p, p_bar);
            let w_arr: [u32;4] = transmute(w_vec);

            //  (1*9)%17=9, (8*9)=72%17=72-68=4, (16*9)=144%17=8, (10*9)=90%17=5
            assert_eq!(w_arr, [9, 4, 8, 5]);
        }
    }
}

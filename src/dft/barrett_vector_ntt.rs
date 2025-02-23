use crate::dft::barrett_field_32bit::{
    add, barrett_mul, barrett_precompute, inv, sub,
    vec_add, vec_sub, vec_barrett_mul_scalar
};
use crate::dft::util::{build_bitrev_tables_u32, find_primitive_2nth_root_of_unity_32};
use crate::dft::DFT;

#[cfg(target_arch = "aarch64")]
use core::arch::aarch64::{
    vdupq_n_u32, vld1q_u32, vst1q_u32,
};

/// NTT implementation for p to a 32-bit prime using Barrett reduction and NEON vectorization.
pub struct BarrettVectorNtt {
    q: u32,
    n: usize,
    psi: u32,
    psi_inv: u32,
    fwd_twid: Vec<u32>,
    inv_twid: Vec<u32>,
    inv_n: u32,
    barrett: u64,
}

impl BarrettVectorNtt {
    pub fn with_params(q: u32, n: usize) -> Option<Self> {
        if !n.is_power_of_two() {
            return None;
        }
        let two_n = 2 * (n as u32);
        if (q - 1) % two_n != 0 {
            return None;
        }

        let (psi, psi_inv) = find_primitive_2nth_root_of_unity_32(q, n)?;
        let inv_n_val = inv(n as u32, q)?;
        let bar = barrett_precompute(q);

        let (fwd_twid, inv_twid) = build_bitrev_tables_u32(n, q, psi, psi_inv);

        Some(Self {
            q,
            n,
            psi,
            psi_inv,
            fwd_twid,
            inv_twid,
            inv_n: inv_n_val,
            barrett: bar,
        })
    }

    #[inline(always)]
    pub fn q(&self) -> u32 {
        self.q
    }

    #[inline(always)]
    pub fn size(&self) -> usize {
        self.n
    }

    #[inline(always)]
    pub fn forward_inplace(&self, data: &mut [u32]) {
        #[cfg(target_arch = "aarch64")]
        unsafe {
            self.forward_inplace_neon(data);
        }
        #[cfg(not(target_arch = "aarch64"))]
        panic!("no neon support!");
    }

    #[inline(always)]
    pub fn backward_inplace(&self, data: &mut [u32]) {
        #[cfg(target_arch = "aarch64")]
        unsafe {
            self.backward_inplace_neon(data);
        }
        #[cfg(not(target_arch = "aarch64"))]
        panic!("no neon support!");
    }

    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "neon")]
    unsafe fn forward_inplace_neon(&self, a: &mut [u32]) {
        let q = self.q;
        let p_bar = self.barrett;
        let n = self.n;
        // broadcast q into a NEON register
        let p_vec = vdupq_n_u32(q);

        let mut half = n;
        let mut step = 1;
        while step < n {
            half >>= 1;
            for i in 0..step {
                let w = self.fwd_twid[step + i];
                let base = 2 * i * half;
                let end = base + half;

                let mut j = base;
                // Vectorized butterfly, 4 elements at a time
                while j + 3 < end {
                    // load 4 items from a[j], a[j+1], a[j+2], a[j+3] into a NEON register
                    let u_vec = vld1q_u32(a.as_ptr().add(j));
                    // load from a[j+half..]
                    let tv_vec = vld1q_u32(a.as_ptr().add(j + half));

                    let v_vec = vec_barrett_mul_scalar(tv_vec, w, q, p_bar);

                    // lane-wise add/sub (mod p) using vec_add / vec_sub
                    let sum_vec = vec_add(u_vec, v_vec, p_vec);
                    let diff_vec = vec_sub(u_vec, v_vec, p_vec);

                    // store back
                    vst1q_u32(a.as_mut_ptr().add(j), sum_vec);
                    vst1q_u32(a.as_mut_ptr().add(j + half), diff_vec);

                    j += 4;
                }
                // handle remainder with scalar loop
                while j < end {
                    let u = unsafe { *a.get_unchecked(j) };
                    let tv = unsafe { *a.get_unchecked(j + half) };

                    let v = barrett_mul(tv, w, q, p_bar);

                    let sum_ = add(u, v, q);
                    let diff_ = sub(u, v, q);

                    unsafe {
                        *a.get_unchecked_mut(j) = sum_;
                        *a.get_unchecked_mut(j + half) = diff_;
                    }
                    j += 1;
                }
            }
            step <<= 1;
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "neon")]
    unsafe fn backward_inplace_neon(&self, a: &mut [u32]) {
        let q = self.q;
        let p_bar = self.barrett;
        let n = self.n;
        let invn_val = self.inv_n;
        let p_vec = vdupq_n_u32(q);

        let mut step = n;
        let mut half = 1;
        while step > 1 {
            let halfstep = step >> 1;
            for i in 0..halfstep {
                let w = self.inv_twid[halfstep + i];
                let base = 2 * i * half;
                let end = base + half;

                let mut j = base;
                while j + 3 < end {
                    let u_vec = vld1q_u32(a.as_ptr().add(j));
                    let tv_vec = vld1q_u32(a.as_ptr().add(j + half));

                    let sum_vec = vec_add(u_vec, tv_vec, p_vec);
                    let tmp_diff_vec = vec_sub(u_vec, tv_vec, p_vec);

                    let diffm_vec = vec_barrett_mul_scalar(tmp_diff_vec, w, q, p_bar);

                    vst1q_u32(a.as_mut_ptr().add(j), sum_vec);
                    vst1q_u32(a.as_mut_ptr().add(j + half), diffm_vec);

                    j += 4;
                }
                while j < end {
                    let u = unsafe { *a.get_unchecked(j) };
                    let v = unsafe { *a.get_unchecked(j + half) };

                    let sum_ = add(u, v, q);
                    let diff_ = sub(u, v, q);
                    let diffm = barrett_mul(diff_, w, q, p_bar);

                    unsafe {
                        *a.get_unchecked_mut(j) = sum_;
                        *a.get_unchecked_mut(j + half) = diffm;
                    }
                    j += 1;
                }
            }
            half <<= 1;
            step = halfstep;
        }
        for idx in 0..n {
            let val = unsafe { *a.get_unchecked(idx) };
            let out = barrett_mul(val, invn_val, q, p_bar);
            unsafe {
                *a.get_unchecked_mut(idx) = out;
            }
        }
    }
}

impl DFT<u32> for BarrettVectorNtt {
    fn forward_inplace(&self, x: &mut [u32]) {
        self.forward_inplace(x);
    }
    fn backward_inplace(&self, x: &mut [u32]) {
        self.backward_inplace(x);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dft::util::{naive_negacyclic_u32, pointwise_u32};
    use rand::thread_rng;
    use rand::Rng;

    #[test]
    fn test_barret_vector_forward_backward() {
        #[cfg(target_arch = "aarch64")]
        {
            let q = 2013265921u32;
            let n = 8;
            let table = BarrettVectorNtt::with_params(q, n).unwrap();

            let mut rng = thread_rng();
            let mut data = vec![0u32; n];
            for x in data.iter_mut() {
                *x = rng.gen_range(0..q);
            }
            let orig = data.clone();

            table.forward_inplace(&mut data);
            table.backward_inplace(&mut data);
            assert_eq!(data, orig);
        }
    }

    #[test]
    fn test_barret_vector_polymul() {
        #[cfg(target_arch = "aarch64")]
        {
            let q = 1062862849u32;
            let n = 8;
            let table = BarrettVectorNtt::with_params(q, n).unwrap();

            let mut rng = thread_rng();
            let mut a = vec![0u32; n];
            let mut b = vec![0u32; n];
            for i in 0..n {
                a[i] = rng.gen_range(0..q);
                b[i] = rng.gen_range(0..q);
            }
            let c_naive = naive_negacyclic_u32(&a, &b, q);

            let mut fa = a.clone();
            let mut fb = b.clone();
            table.forward_inplace(&mut fa);
            table.forward_inplace(&mut fb);

            pointwise_u32(&mut fa, &fb, q);

            table.backward_inplace(&mut fa);

            assert_eq!(fa, c_naive);
        }
    }
}

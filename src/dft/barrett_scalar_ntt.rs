use rand::Rng;
use crate::dft::DFT;
use crate::dft::barrett_field_32bit::{
    barrett_precompute, barrett_mul, add, sub, mul, exp, inv,
};
use crate::dft::util::{
    find_primitive_2nth_root_of_unity_32,
    build_bitrev_tables_u32,
};

pub struct BarrettScalarNtt {
    q: u32,
    n: usize,
    psi: u32,
    psi_inv: u32,
    fwd_twid: Vec<u32>,
    inv_twid: Vec<u32>,
    inv_n: u32,
    barrett: u64,
}

impl BarrettScalarNtt {
    pub fn with_params(q: u32, n: usize) -> Option<Self> {
        if !n.is_power_of_two() {
            return None;
        }
        let two_n = 2*(n as u32);
        if (q-1) % two_n != 0 {
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
    pub fn q(&self) -> u32 { self.q }

    #[inline(always)]
    pub fn size(&self) -> usize { self.n }

    #[inline(always)]
    pub fn forward_inplace(&self, a: &mut [u32]) {
        let p_bar = self.barrett;
        let q = self.q;

        let mut half = self.n;
        let mut step = 1;
        while step < self.n {
            half >>= 1;
            for i in 0..step {
                let w = self.fwd_twid[step + i];
                let base = 2*i*half;
                let end  = base + half;
                for j in base..end {
                    let u  = unsafe { *a.get_unchecked(j) };
                    let tv = unsafe { *a.get_unchecked(j+half) };

                    let v = barrett_mul(tv, w, q, p_bar);

                    let sum_  = add(u, v, q);
                    let diff_ = sub(u, v, q);

                    unsafe {
                        *a.get_unchecked_mut(j)        = sum_;
                        *a.get_unchecked_mut(j + half) = diff_;
                    }
                }
            }
            step <<= 1;
        }
    }

    #[inline(always)]
    pub fn backward_inplace(&self, a: &mut [u32]) {
        let p_bar = self.barrett;
        let q = self.q;
        let invn_val = self.inv_n;

        let mut step = self.n;
        let mut half = 1;
        while step > 1 {
            let halfstep = step >> 1;
            for i in 0..halfstep {
                let w = self.inv_twid[halfstep + i];
                let base = 2*i*half;
                let end  = base + half;
                for j in base..end {
                    let u = unsafe { *a.get_unchecked(j) };
                    let v = unsafe { *a.get_unchecked(j+half) };

                    let sum_  = add(u, v, q);
                    let diff_ = sub(u, v, q);
                    let diffm = barrett_mul(diff_, w, q, p_bar);

                    unsafe {
                        *a.get_unchecked_mut(j)        = sum_;
                        *a.get_unchecked_mut(j + half) = diffm;
                    }
                }
            }
            half <<= 1;
            step = halfstep;
        }
        for j in 0..self.n {
            let val = unsafe { *a.get_unchecked(j) };
            let out = barrett_mul(val, invn_val, q, p_bar);
            unsafe {
                *a.get_unchecked_mut(j) = out;
            }
        }
    }
}

impl DFT<u32> for BarrettScalarNtt {
    fn forward_inplace(&self, x: &mut [u32]) { self.forward_inplace(x); }
    fn backward_inplace(&self, x: &mut [u32]) { self.backward_inplace(x); }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dft::util::{ naive_negacyclic_u32, pointwise_u32 };

    use rand::thread_rng;
    use rand::Rng;

    #[test]
    fn test_barret_scalar_forward_backward() {
        let q= 2013265921u32;
        let n= 8;
        let table= BarrettScalarNtt::with_params(q, n).expect("cannot build");

        let mut rng= thread_rng();
        let mut data= vec![0u32; n];
        for x in data.iter_mut() {
            *x= rng.gen_range(0..q);
        }
        let orig= data.clone();

        table.forward_inplace(&mut data);
        table.backward_inplace(&mut data);
        assert_eq!(data, orig);
    }

    #[test]
    fn test_barret_scalar_polymul() {
        let q= 1062862849u32;
        let n= 8;
        let table= BarrettScalarNtt::with_params(q, n).expect("cannot build");

        let mut rng= thread_rng();
        let mut a= vec![0u32; n];
        let mut b= vec![0u32; n];
        for i in 0..n {
            a[i]= rng.gen_range(0..q);
            b[i]= rng.gen_range(0..q);
        }
        let c_naive= naive_negacyclic_u32(&a, &b, q);

        let mut fa= a.clone();
        let mut fb= b.clone();
        table.forward_inplace(&mut fa);
        table.forward_inplace(&mut fb);

        pointwise_u32(&mut fa, &fb, q);

        table.backward_inplace(&mut fa);

        assert_eq!(fa, c_naive);
    }
}

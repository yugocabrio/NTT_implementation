#[path = "dft/mont_ntt.rs"]
pub mod mont_ntt;
#[path = "dft/field.rs"]
pub mod field;
pub mod mont_field;
pub mod shoup_field;
pub mod shoup_ntt;
pub mod goldilocks_field;
pub mod goldilocks_ntt;

pub trait DFT<O> {
    fn forward_inplace(&self, x: &mut [O]);
    fn backward_inplace(&self, x: &mut[O]);
}
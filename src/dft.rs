pub mod ntt;
#[path = "dft/field.rs"]
pub mod field;
pub mod mont_field;

pub trait DFT<O> {
    fn forward_inplace(&self, x: &mut [O]);
    fn forward_inplace_lazy(&self, x: &mut [O]);
    fn backward_inplace(&self, x: &mut[O]);
    fn backward_inplace_lazy(&self, x: &mut[O]);
}
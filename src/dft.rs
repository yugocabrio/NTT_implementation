pub mod barrett_field_32bit;
pub mod barrett_scalar_ntt;
pub mod barrett_vector_ntt;
#[path = "dft/field.rs"]
pub mod field;
pub mod goldilocks_field;
pub mod goldilocks_ntt;
pub mod mont_field;
#[path = "dft/mont_ntt.rs"]
pub mod mont_ntt;
pub mod shoup_field;
pub mod shoup_ntt;
pub mod util;

pub trait DFT<O> {
    fn forward_inplace(&self, x: &mut [O]);
    fn backward_inplace(&self, x: &mut [O]);
}

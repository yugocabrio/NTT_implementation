use crate::dft::DFT;

pub struct Table<O>{
    q:O, // NTT friendly prime modulus
    psi:O, // n-th root of unity
}

impl Table<u64> {
    pub fn new()->Self{
        Self{
            q:0x1fffffffffe00001u64, 
            psi:0x15eb043c7aa2b01fu64, //2^17th root of unity
        }
    }
}

impl DFT<u64> for Table<u64>{
    fn forward_inplace(&self, a: &mut [u64]){
        self.forward_inplace_core::<false>(a) 
    }

    fn forward_inplace_lazy(&self, a: &mut [u64]){
         self.forward_inplace_core::<true>(a) 
    }

    fn backward_inplace(&self, a: &mut [u64]){
        self.backward_inplace_core::<false>(a) 
    }

    fn backward_inplace_lazy(&self, a: &mut [u64]){
         self.backward_inplace_core::<true>(a) 
    }
}

impl Table<u64>{

    pub fn forward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {
    }

    pub fn backward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {
    }
}

use crate::dft::DFT;

pub struct Table<O> {
    /// NTT friendly prime modulus
    q: O,
    /// n-th root of unity
    psi: O,
}

impl Table<u64> {
    pub fn new() -> Self {
        Self {
            q: 0x1fffffffffe00001u64,
            psi: 0x15eb043c7aa2b01fu64, //2^17th root of unity
        }
    }
}

impl DFT<u64> for Table<u64> {
    /// NTT forward routine
    ///
    /// - `a`: vector with each element in range `[0, q)`
    fn forward_inplace(&self, a: &mut [u64]) {
        self.forward_inplace_core::<false>(a)
    }

    /// NTT forward lazy routine
    ///
    /// - `a`: vector with each element in range `[0, 2q)`
    fn forward_inplace_lazy(&self, a: &mut [u64]) {
        self.forward_inplace_core::<true>(a)
    }

    /// NTT backward routine
    ///
    /// - `a`: vector with each element in range `[0, q)`
    fn backward_inplace(&self, a: &mut [u64]) {
        self.backward_inplace_core::<false>(a)
    }

    /// NTT backward lazy routine
    ///
    /// - `a`: vector with each element in range `[0, 2q)`
    fn backward_inplace_lazy(&self, a: &mut [u64]) {
        self.backward_inplace_core::<true>(a)
    }
}

impl Table<u64> {
    pub fn forward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {}
    pub fn backward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {}
}

fn bitreverse(mut x: usize, log_n: usize) -> usize {
    let mut y = 0;
    for _ in 0..log_n {
        y = (y << 1) | (x & 1);
        x >>= 1;
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitreverse() {
        assert_eq!(bitreverse(6, 3), 3);
    }
}
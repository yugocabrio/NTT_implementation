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

/// (a + b) mod m
fn add(a: u64, b: u64, m: u64) -> u64 {
    let add_num = a + b;
    if add_num >= m { add_num - m } else { add_num }
}

/// (a - b) mod m
fn sub(a: u64, b: u64, m: u64) -> u64 {
    if a >= b { a - b } else { a + m - b }
}

/// (a * b) mod m
fn mul(a: u64, b: u64, m:u64) -> u64 {
    let mul_num = (a as u128) * (b as u128);
    (mul_num % (m as u128)) as u64
}

/// a^b mod m
fn exp(base: u64, exp: u64, m: u64) -> u64 {
    let mut exp_num  = 1u64;
    let mut current_base = base % m;
    let mut e = exp;
    while e > 0 {
        if (e & 1) == 1 {
            exp_num = mul(exp_num, current_base, m);
        }
        current_base = mul(current_base, current_base, m);
        e >>= 1;
    }
    exp_num
}

/// a^(m-2) mod m, Fermat's theorem
fn inv(a: u64, m: u64) -> Option<u64> {
    if a == 0 {
        None
    } else {
        Some(exp(a, m-2, m))
    }
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
    fn test_add() {
        let q = 5u64;
        assert_eq!(add(4, 1, q), 0);
    }

    #[test]
    fn test_sub() {
        let q = 5u64;
        assert_eq!(sub(2, 3, q), 4);
    }

    #[test]
    fn test_mul() {
        let q = 5u64;
        assert_eq!(mul(2, 3, q), 1);
    }

    #[test]
    fn test_exp() {
        let q = 5u64;
        assert_eq!(exp(2, 4, q), 1);
    }

    #[test]
    fn test_modinv() {
        let q = 5u64;
        assert_eq!(inv(4, q), Some(4));
        assert_eq!(inv(0, q), None);
    }

    #[test]
    fn test_bitreverse() {
        assert_eq!(bitreverse(6, 3), 3);
    }
}
use criterion::{
    black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion
};
use rand::Rng;
use app::dft::DFT;
use app::dft::mont_ntt::MontTable;
use app::dft::shoup_ntt::ShoupTable;
use concrete_ntt::prime64;
use app::dft::goldilocks_ntt::{GoldilocksNttTable};
use app::dft::goldilocks_field::{reduce128, GOLDILOCKS_P};

use plonky2_field::{
    fft::{fft, ifft},
    goldilocks_field::GoldilocksField,
    polynomial::{PolynomialCoeffs, PolynomialValues},
    types::Sample,
};

/// 61 bit prime (Mont, Shoup, concrete-ntt)
const PRIME: u64 = 0x1fffffffffe00001;
use app::dft::util::pointwise_u64;

/// forward benches
fn bench_ntt_compare(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_compare");

    for log_n in 11..=16 {
        let n = 1 << log_n;

        // Mont
        if let Some(mont_table) = MontTable::with_params(PRIME, n) {
            let bench_id = BenchmarkId::new("mont-forward_inplace", n);
            group.bench_with_input(bench_id, &n, |b, &_| {
                b.iter_batched(
                    || {
                        let mut rng = rand::thread_rng();
                        let mut data = vec![0u64; n];
                        for x in data.iter_mut() {
                            *x = rng.gen_range(0..PRIME);
                        }
                        data
                    },
                    |mut data| {
                        mont_table.forward_inplace(black_box(&mut data));
                    },
                    BatchSize::LargeInput,
                );
            });
        }

        // Shoup
        if let Some(shoup_table) = ShoupTable::with_params(PRIME, n) {
            let bench_id = BenchmarkId::new("shoup-forward_inplace", n);
            group.bench_with_input(bench_id, &n, |b, &_| {
                b.iter_batched(
                    || {
                        let mut rng = rand::thread_rng();
                        let mut data = vec![0u64; n];
                        for x in data.iter_mut() {
                            *x = rng.gen_range(0..PRIME);
                        }
                        data
                    },
                    |mut data| {
                        shoup_table.forward_inplace(black_box(&mut data));
                    },
                    BatchSize::LargeInput,
                );
            });
        }

        // concrete-ntt
        if let Some(plan) = prime64::Plan::try_new(n, PRIME) {
            let bench_id = BenchmarkId::new("concrete-forward", n);
            group.bench_with_input(bench_id, &n, |b, &_| {
                b.iter_batched(
                    || {
                        let mut rng = rand::thread_rng();
                        let mut data = vec![0u64; n];
                        for x in data.iter_mut() {
                            *x = rng.gen_range(0..PRIME);
                        }
                        data
                    },
                    |mut data| {
                        plan.fwd(black_box(&mut data));
                    },
                    BatchSize::LargeInput,
                );
            });
        }

        // Goldilocks
        if (GOLDILOCKS_P - 1) % (2*(n as u64)) == 0 {
            if let Some(goldi_table) = GoldilocksNttTable::with_params(n) {
                let bench_id = BenchmarkId::new("goldilocks-forward_inplace", n);
                group.bench_with_input(bench_id, &n, |b, &_| {
                    b.iter_batched(
                        || {
                            let mut rng = rand::thread_rng();
                            let mut data = vec![0u64; n];
                            for x in data.iter_mut() {
                                *x = rng.gen_range(0..GOLDILOCKS_P);
                            }
                            data
                        },
                        |mut data| {
                            goldi_table.forward_inplace(black_box(&mut data));
                        },
                        BatchSize::LargeInput,
                    );
                });
            }
        }

        // Plonky2
        {
            let bench_id = BenchmarkId::new("plonky2-forward", n);
            group.bench_with_input(bench_id, &n, |b, &_| {
                b.iter_batched(
                    || GoldilocksField::rand_vec(n),
                    |coeffs| {
                        let poly = PolynomialCoeffs { coeffs };
                        let values = fft(poly);
                        black_box(values);
                    },
                    BatchSize::LargeInput,
                );
            });
        }
    }
    group.finish();
}

/// poly mul
fn bench_ntt_polymul_compare(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_polymul_compare");

    for log_n in 11..=16 {
        let n = 1 << log_n;

        // Mont poly mul
        if let Some(mont_table) = MontTable::with_params(PRIME, n) {
            let bench_id = BenchmarkId::new("mont-polymul", n);
            group.bench_with_input(bench_id, &n, |b, &_| {
                b.iter_batched(
                    || {
                        let mut rng = rand::thread_rng();
                        let mut a = vec![0u64; n];
                        let mut b = vec![0u64; n];
                        for x in &mut a {
                            *x = rng.gen_range(0..PRIME);
                        }
                        for x in &mut b {
                            *x = rng.gen_range(0..PRIME);
                        }
                        (a, b)
                    },
                    |(mut a, mut b)| {
                        let mut fa = a.clone();
                        let mut fb = b.clone();
                        mont_table.forward_inplace(&mut fa);
                        mont_table.forward_inplace(&mut fb);

                        pointwise_u64(&mut fa, &fb, PRIME);

                        mont_table.backward_inplace(&mut fa);
                        black_box(&fa);
                    },
                    BatchSize::LargeInput,
                );
            });
        }

        // Shoup poly mul
        if let Some(shoup_table) = ShoupTable::with_params(PRIME, n) {
            let bench_id = BenchmarkId::new("shoup-polymul", n);
            group.bench_with_input(bench_id, &n, |b, &_| {
                b.iter_batched(
                    || {
                        let mut rng = rand::thread_rng();
                        let mut a = vec![0u64; n];
                        let mut b = vec![0u64; n];
                        for x in &mut a {
                            *x = rng.gen_range(0..PRIME);
                        }
                        for x in &mut b {
                            *x = rng.gen_range(0..PRIME);
                        }
                        (a, b)
                    },
                    |(mut a, mut b)| {
                        let mut fa = a.clone();
                        let mut fb = b.clone();
                        shoup_table.forward_inplace(&mut fa);
                        shoup_table.forward_inplace(&mut fb);

                        pointwise_u64(&mut fa, &fb, PRIME);

                        shoup_table.backward_inplace(&mut fa);
                        black_box(&fa);
                    },
                    BatchSize::LargeInput,
                );
            });
        }

        // concrete-ntt pol mul
        if let Some(plan) = prime64::Plan::try_new(n, PRIME) {
            let bench_id = BenchmarkId::new("concrete-polymul", n);
            group.bench_with_input(bench_id, &n, |b, &_| {
                b.iter_batched(
                    || {
                        let mut rng = rand::thread_rng();
                        let mut a = vec![0u64; n];
                        let mut b = vec![0u64; n];
                        for x in &mut a {
                            *x = rng.gen_range(0..PRIME);
                        }
                        for x in &mut b {
                            *x = rng.gen_range(0..PRIME);
                        }
                        (a, b)
                    },
                    |(mut a, mut b)| {
                        plan.fwd(&mut a);
                        plan.fwd(&mut b);

                        plan.mul_assign_normalize(&mut a, &b);
                        plan.inv(&mut a);
                        black_box(&a);
                    },
                    BatchSize::LargeInput,
                );
            });
        }

        // Goldilocks poly mul
        if (GOLDILOCKS_P - 1) % (2*(n as u64)) == 0 {
            if let Some(gold_table) = GoldilocksNttTable::with_params(n) {
                let bench_id= BenchmarkId::new("goldilocks-polymul", n);
                group.bench_with_input(bench_id, &n, |b, &_| {
                    b.iter_batched(
                        || {
                            let mut rng = rand::thread_rng();
                            let mut a = vec![0u64; n];
                            let mut b = vec![0u64; n];
                            for x in &mut a {
                                *x = rng.gen_range(0..GOLDILOCKS_P);
                            }
                            for x in &mut b {
                                *x = rng.gen_range(0..GOLDILOCKS_P);
                            }
                            (a,b)
                        },
                        |(mut a, mut b)| {
                            gold_table.forward_inplace(&mut a);
                            gold_table.forward_inplace(&mut b);

                            for i in 0..n {
                                let prod = (a[i] as u128) * (b[i] as u128);
                                a[i] = reduce128(prod);
                            }

                            gold_table.backward_inplace(&mut a);
                            black_box(&a);
                        },
                        BatchSize::LargeInput,
                    );
                });
            }
        }

        // Plonky2 poly mul
        {
            let bench_id = BenchmarkId::new("plonky2-polymul", n);
            group.bench_with_input(bench_id, &n, |b, &_| {
                b.iter_batched(
                    || {
                        let coeffs_a = GoldilocksField::rand_vec(n);
                        let coeffs_b = GoldilocksField::rand_vec(n);
                        (coeffs_a, coeffs_b)
                    },
                    |(coeffs_a, coeffs_b)| {
                        let poly_a = PolynomialCoeffs { coeffs: coeffs_a };
                        let poly_b = PolynomialCoeffs { coeffs: coeffs_b };

                        let vals_a = fft(poly_a);
                        let vals_b = fft(poly_b);

                        let mut vals_c = vals_a.values;
                        for i in 0..n {
                            vals_c[i] *= vals_b.values[i];
                        }

                        let result = ifft(PolynomialValues { values: vals_c });
                        black_box(result);
                    },
                    BatchSize::LargeInput,
                );
            });
        }
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_ntt_compare,
    bench_ntt_polymul_compare,
);
criterion_main!(benches);

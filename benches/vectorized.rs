use criterion::{
    black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion
};
use rand::Rng;

use app::dft::barrett_scalar_ntt::BarrettScalarNtt; 
use app::dft::barrett_vector_ntt::BarrettVectorNtt;
use concrete_ntt::prime32;
use app::dft::util::pointwise_u32;

const PRIME_32: u32 = 1062862849;

/// 1) forward のベンチマーク
fn bench_ntt_compare(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_forward_compare");

    // n=65536
    let log_n = 16;
    let n = 1 << log_n;

    // Barrett Scalar
    if let Some(table_scalar) = BarrettScalarNtt::with_params(PRIME_32, n) {
        let bench_id = BenchmarkId::new("32_scalar-forward_inplace", n);
        group.bench_with_input(bench_id, &n, |b, &_n| {
            b.iter_batched(
                || {
                    let mut rng = rand::thread_rng();
                    let mut data = vec![0u32; n];
                    for x in &mut data {
                        *x = rng.gen_range(0..PRIME_32);
                    }
                    data
                },
                |mut data| {
                    table_scalar.forward_inplace(black_box(&mut data));
                },
                BatchSize::LargeInput,
            );
        });
    }

    // Barrett Vector
    if let Some(table_vector) = BarrettVectorNtt::with_params(PRIME_32, n) {
        let bench_id = BenchmarkId::new("32_vector-forward_inplace", n);
        group.bench_with_input(bench_id, &n, |b, &_n| {
            b.iter_batched(
                || {
                    let mut rng = rand::thread_rng();
                    let mut data = vec![0u32; n];
                    for x in &mut data {
                        *x = rng.gen_range(0..PRIME_32);
                    }
                    data
                },
                |mut data| {
                    table_vector.forward_inplace(black_box(&mut data));
                },
                BatchSize::LargeInput,
            );
        });
    }

    // concrete 32 bit
    if let Some(plan) = prime32::Plan::try_new(n, PRIME_32) {
        let bench_id = BenchmarkId::new("concrete32-forward", n);
        group.bench_with_input(bench_id, &n, |b, &_n| {
            b.iter_batched(
                || {
                    let mut rng = rand::thread_rng();
                    let mut data = vec![0u32; n];
                    for x in &mut data {
                        *x = rng.gen_range(0..PRIME_32);
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

    group.finish();
}

/// polynomial multiplication
fn bench_ntt_polymul_compare(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_polymul_compare");

    let log_n = 16;
    let n = 1 << log_n;

    // Barrett Scalar
    if let Some(table_scalar) = BarrettScalarNtt::with_params(PRIME_32, n) {
        let bench_id = BenchmarkId::new("32_scalar-polymul", n);
        group.bench_with_input(bench_id, &n, |b, &_n| {
            b.iter_batched(
                || {
                    let mut rng = rand::thread_rng();
                    let mut a = vec![0u32; n];
                    let mut b = vec![0u32; n];
                    for x in &mut a {
                        *x = rng.gen_range(0..PRIME_32);
                    }
                    for x in &mut b {
                        *x = rng.gen_range(0..PRIME_32);
                    }
                    (a, b)
                },
                |(a, b)| {
                    let mut fa = a.clone();
                    let mut fb = b.clone();
                    table_scalar.forward_inplace(&mut fa);
                    table_scalar.forward_inplace(&mut fb);

                    pointwise_u32(&mut fa, &fb, PRIME_32);

                    table_scalar.backward_inplace(&mut fa);

                    black_box(&fa);
                },
                BatchSize::LargeInput,
            );
        });
    }

    // Barrett Vector
    if let Some(table_vector) = BarrettVectorNtt::with_params(PRIME_32, n) {
        let bench_id = BenchmarkId::new("32_neon-polymul", n);
        group.bench_with_input(bench_id, &n, |b, &_n| {
            b.iter_batched(
                || {
                    let mut rng = rand::thread_rng();
                    let mut a = vec![0u32; n];
                    let mut b = vec![0u32; n];
                    for x in &mut a {
                        *x = rng.gen_range(0..PRIME_32);
                    }
                    for x in &mut b {
                        *x = rng.gen_range(0..PRIME_32);
                    }
                    (a, b)
                },
                |(a, b)| {
                    let mut fa = a.clone();
                    let mut fb = b.clone();
                    table_vector.forward_inplace(&mut fa);
                    table_vector.forward_inplace(&mut fb);

                    pointwise_u32(&mut fa, &fb, PRIME_32);

                    table_vector.backward_inplace(&mut fa);
                    black_box(&fa);
                },
                BatchSize::LargeInput,
            );
        });
    }

    // concrete-ntt (32bit)
    if let Some(plan) = prime32::Plan::try_new(n, PRIME_32) {
        let bench_id = BenchmarkId::new("concrete32-polymul", n);
        group.bench_with_input(bench_id, &n, |b, &_n| {
            b.iter_batched(
                || {
                    let mut rng = rand::thread_rng();
                    let mut a = vec![0u32; n];
                    let mut b = vec![0u32; n];
                    for x in &mut a {
                        *x = rng.gen_range(0..PRIME_32);
                    }
                    for x in &mut b {
                        *x = rng.gen_range(0..PRIME_32);
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

    group.finish();
}

criterion_group!(
    benches,
    bench_ntt_compare,
    bench_ntt_polymul_compare,
);
criterion_main!(benches);

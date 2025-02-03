use criterion::{criterion_group, criterion_main, Criterion, black_box};
use rand::Rng;

use app::dft::field as naive;
use app::dft::mont_field::*;

const BIG_Q: u64 = 0x1fffffffffe00001;
const K: u32 = 62;

fn bench_naive_mul(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    c.bench_function("naive_mul_61bit", |b| {
        b.iter(|| {
            let a = black_box(rng.gen_range(0..BIG_Q));
            let b_ = black_box(rng.gen_range(0..BIG_Q));
            let _ = naive::mul(a, b_, BIG_Q);
        })
    });
}

fn bench_mont_mul(c: &mut Criterion) {
    let ctx = MontgomeryContext::new(BIG_Q, K);
    let mut rng = rand::thread_rng();

    c.bench_function("mont_mul_61bit", |b| {
        b.iter(|| {
            let a = rng.gen_range(0..BIG_Q);
            let b_ = rng.gen_range(0..BIG_Q);
            let a_mont = ctx.to_mont(a);
            let b_mont = ctx.to_mont(b_);
            let _ = mont_mul(a_mont, b_mont, &ctx);
        })
    });
}

fn bench_naive_exp(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    c.bench_function("naive_exp_61bit", |b| {
        b.iter(|| {
            let base = black_box(rng.gen_range(0..BIG_Q));
            let exp_ = black_box(rng.gen_range(0..10000)); 
            let _ = naive::exp(base, exp_, BIG_Q);
        })
    });
}

fn bench_mont_exp(c: &mut Criterion) {
    let ctx = MontgomeryContext::new(BIG_Q, K);
    let mut rng = rand::thread_rng();

    c.bench_function("mont_exp_61bit", |b| {
        b.iter(|| {
            let base = rng.gen_range(1..BIG_Q);
            let exp_ = rng.gen_range(0..10000);
            // to mont
            let base_m = ctx.to_mont(base);
            let _ = mont_exp(base_m, exp_, &ctx);
        })
    });
}

fn bench_naive_inv(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    c.bench_function("naive_inv_61bit", |b| {
        b.iter(|| {
            let a = black_box(rng.gen_range(1..BIG_Q));
            let _ = naive::inv(a, BIG_Q);
        })
    });
}

fn bench_mont_inv(c: &mut Criterion) {
    let ctx = MontgomeryContext::new(BIG_Q, K);
    let mut rng = rand::thread_rng();

    c.bench_function("mont_inv_61bit", |b| {
        b.iter(|| {
            let a = rng.gen_range(1..BIG_Q);
            let a_m = ctx.to_mont(a);
            let _ = mont_inv(a_m, &ctx);
        })
    });
}

fn criterion_benches(c: &mut Criterion) {
    bench_naive_mul(c);
    bench_mont_mul(c);
    bench_naive_exp(c);
    bench_mont_exp(c);
    bench_naive_inv(c);
    bench_mont_inv(c);
}

criterion_group!(field_benches, criterion_benches);
criterion_main!(field_benches);

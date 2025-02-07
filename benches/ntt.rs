use criterion::{
    black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion
};
use rand::Rng;
use app::dft::DFT;
use app::dft::mont_ntt::MontTable;
use app::dft::shoup_ntt::ShoupTable;
use concrete_ntt::prime64;

const PRIME: u64 = 0x1fffffffffe00001;

fn bench_ntt_compare(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_compare");

    for log_n in 16..=16 {
        let n = 1 << log_n;

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
    }

    group.finish();
}

fn bench_ntt_polymul_compare(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_polymul_compare");

    for log_n in 16..=16 {
        let n = 1 << log_n;

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
                        mont_table.forward_inplace(&mut a);
                        mont_table.forward_inplace(&mut b);

                        for i in 0..n {
                            let prod = (a[i] as u128 * b[i] as u128) % (PRIME as u128);
                            a[i] = prod as u64;
                        }

                        mont_table.backward_inplace(&mut a);

                        black_box(&a);
                    },
                    BatchSize::LargeInput,
                );
            });
        }

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
                        shoup_table.forward_inplace(&mut a);
                        shoup_table.forward_inplace(&mut b);

                        for i in 0..n {
                            let prod = (a[i] as u128 * b[i] as u128) % (PRIME as u128);
                            a[i] = prod as u64;
                        }

                        shoup_table.backward_inplace(&mut a);

                        black_box(&a);
                    },
                    BatchSize::LargeInput,
                );
            });
        }

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
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_ntt_compare,
    bench_ntt_polymul_compare,
);
criterion_main!(benches);

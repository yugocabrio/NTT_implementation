use criterion::{
    black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion
};
use rand::Rng;
use app::dft::DFT;
use app::dft::ntt::Table as MontTable;
use concrete_ntt::prime64;

const PRIME: u64 = 0x1fffffffffe00001;

fn bench_ntt_compare(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_compare");

    for log_n in 11..=16 {
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

criterion_group!(benches, bench_ntt_compare);
criterion_main!(benches);

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use app::dft::ntt::Table;
use app::dft::DFT;

fn forward_inplace(c: &mut Criterion) {
    fn runner(log_n:usize) -> Box<dyn FnMut()> {
        let ntt_table: Table<u64> = Table::<u64>::new();
        let mut a: Vec<u64> = vec![0; 1<<log_n as usize];
        for i in 0..a.len(){
            a[i] = i as u64;
        }
        Box::new(move || {
            ntt_table.forward_inplace(&mut a)
        })
    }

    let mut b: criterion::BenchmarkGroup<'_, criterion::measurement::WallTime> = c.benchmark_group("forward_inplace");
    for log_n in 11..17 {

        let runners: [(&str, Box<dyn FnMut()>); 1] = [
            ("prime", {
                runner(log_n)
            }),
        ];
        for (name, mut runner) in runners {
            let id = BenchmarkId::new(name, 1<<log_n);
            b.bench_with_input(id, &(), |b: &mut criterion::Bencher<'_>, _| b.iter(&mut runner));
        }
    }
}

criterion_group!(benches, forward_inplace);
criterion_main!(benches);

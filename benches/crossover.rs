use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nalgebra::Vector3;
use rayx::{intersect, Ray, Triangle};

fn bench_crossover(c: &mut Criterion) {
    let mut group = c.benchmark_group("crossover");

    let v1 = Vector3::new(0.0, 0.0, 0.0);
    let v2 = Vector3::new(1.0, 0.0, 0.0);
    let v3 = Vector3::new(0.0, 1.0, 0.0);

    let ray = Ray::new(Vector3::new(0.1, 0.1, 1.0), Vector3::new(0.0, 0.0, -1.0));

    let t_min = 0.0;
    let t_max = 10.0;

    // We vary the number of rays per triangle to find where BW (init + intersect)
    // outperforms MT (intersect only).
    for n in [1, 2, 3, 4, 5, 8, 10, 20, 50].iter() {
        group.bench_with_input(BenchmarkId::new("moller_trumbore", n), n, |b, &n| {
            b.iter(|| {
                let mut hits = 0;
                for _ in 0..n {
                    if intersect::intersect_moller_trumbore(
                        black_box(v1),
                        black_box(v2),
                        black_box(v3),
                        black_box(ray),
                        black_box(t_min),
                        black_box(t_max),
                    )
                    .is_some()
                    {
                        hits += 1;
                    }
                }
                hits
            })
        });

        group.bench_with_input(BenchmarkId::new("baldwin_weber_total", n), n, |b, &n| {
            b.iter(|| {
                let mut hits = 0;
                // BW Init
                let tri = Triangle::new(black_box(v1), black_box(v2), black_box(v3)).unwrap();
                // BW Intersect
                for _ in 0..n {
                    if tri
                        .intersect(black_box(ray), black_box(t_min), black_box(t_max))
                        .is_some()
                    {
                        hits += 1;
                    }
                }
                hits
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_crossover);
criterion_main!(benches);

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use nalgebra::Vector3;
use rayx::{Ray, Triangle, intersect};

fn bench_crossover(c: &mut Criterion) {
    // A new ID avoids comparing different ray workloads in CodSpeed history.
    let mut group = c.benchmark_group("crossover_diverse_rays");

    let v1 = Vector3::new(0.0, 0.0, 0.0);
    let v2 = Vector3::new(1.0, 0.0, 0.0);
    let v3 = Vector3::new(0.0, 1.0, 0.0);

    // Distinct rays prevent repeated identical intersections from dominating the
    // amortization comparison; both methods receive the exact same workload.
    let rays: Vec<_> = (0..50)
        .map(|i| {
            Ray::new(
                Vector3::new(
                    0.05 + (i % 10) as f64 * 0.1,
                    0.05 + (i / 10) as f64 * 0.1,
                    1.0,
                ),
                Vector3::new(
                    if i % 2 == 0 { 0.05 } else { -0.05 },
                    if i % 4 < 2 { 0.05 } else { -0.05 },
                    -1.0,
                ),
            )
        })
        .collect();

    let t_min = 0.0;
    let t_max = 10.0;

    // We vary the number of rays per triangle to find where BW (init + intersect)
    // outperforms MT (intersect only).
    for n in [1, 2, 3, 4, 5, 8, 10, 20, 50].iter() {
        group.bench_with_input(BenchmarkId::new("moller_trumbore", n), n, |b, &n| {
            b.iter(|| {
                let (v1, v2, v3) = black_box((v1, v2, v3));
                let mut hits = 0;
                for ray in black_box(&rays[..n]) {
                    if intersect::intersect_moller_trumbore(v1, v2, v3, *ray, t_min, t_max)
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
                let (v1, v2, v3) = black_box((v1, v2, v3));
                let mut hits = 0;
                // BW Init
                let tri = Triangle::new(v1, v2, v3).unwrap();
                // BW Intersect
                for ray in black_box(&rays[..n]) {
                    if tri.intersect(*ray, t_min, t_max).is_some() {
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

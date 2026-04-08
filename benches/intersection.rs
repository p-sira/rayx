use criterion::{Criterion, black_box, criterion_group, criterion_main};
use nalgebra::Vector3;
use rayx::{Ray, Triangle, intersect};

fn bench_intersection(c: &mut Criterion) {
    let v0 = Vector3::new(-1.0, 0.0, 0.5);
    let v1 = Vector3::new(2.0, 0.3, -0.2);
    let v2 = Vector3::new(0.2, 1.7, 0.1);
    let tri = Triangle::new(v0, v1, v2).unwrap();

    let ray = Ray::new(Vector3::new(0.1, 0.2, 2.0), Vector3::new(0.0, 0.0, -1.0));

    let mut group = c.benchmark_group("intersect");

    group.bench_function("baldwin_weber_precomputed", |b| {
        b.iter(|| tri.intersect(black_box(ray), black_box(0.0), black_box(1e9)))
    });

    group.bench_function("moller_trumbore", |b| {
        b.iter(|| {
            intersect::intersect_moller_trumbore(
                black_box(v0),
                black_box(v1),
                black_box(v2),
                black_box(ray),
                black_box(0.0),
                black_box(1e9),
            )
        })
    });

    group.finish();

    c.bench_function("baldwin_weber_init", |b| {
        b.iter(|| Triangle::new(black_box(v0), black_box(v1), black_box(v2)))
    });

    let mut group = c.benchmark_group("init_and_intersect");

    group.bench_function("baldwin_weber", |b| {
        b.iter(|| {
            let tri = Triangle::new(black_box(v0), black_box(v1), black_box(v2)).unwrap();
            tri.intersect(black_box(ray), black_box(0.0), black_box(1e9))
        })
    });

    group.bench_function("moller_trumbore", |b| {
        b.iter(|| {
            intersect::intersect_moller_trumbore(
                black_box(v0),
                black_box(v1),
                black_box(v2),
                black_box(ray),
                black_box(0.0),
                black_box(1e9),
            )
        })
    });

    group.finish();
}

criterion_group!(benches, bench_intersection);
criterion_main!(benches);

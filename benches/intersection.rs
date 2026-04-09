use criterion::{Criterion, criterion_group, criterion_main};

mod common;
use common::ModelBench;

fn bench_models(c: &mut Criterion) {
    let models = [
        (
            "suzanne",
            "tests/test_data/perfect-suzanne.stl",
            (-1.5, 1.5),
            (-1.0, 1.0),
            5.0,
            None,
        ),
        (
            "bunny_res2",
            "tests/test_data/bun_zipper_res2.stl",
            (-0.1, 0.07),
            (0.03, 0.19),
            1.0,
            Some(30),
        ),
        (
            "bunny_res3",
            "tests/test_data/bun_zipper_res3.stl",
            (-0.1, 0.07),
            (0.03, 0.19),
            1.0,
            Some(90),
        ),
        (
            "bunny_res4",
            "tests/test_data/bun_zipper_res4.stl",
            (-0.1, 0.07),
            (0.03, 0.19),
            1.0,
            Some(100),
        ),
        (
            "thai_statue_50",
            "tests/test_data/thai-statue-decim001.stl",
            (-75.0, 75.0),
            (-100.0, 100.0),
            200.0,
            Some(10),
        ),
        (
            "thai_statue_80",
            "tests/test_data/thai-statue-decim001.stl",
            (-50.0, 50.0),
            (-100.0, 100.0),
            200.0,
            Some(10),
        ),
    ];

    for (name, path, bounds_x, bounds_y, z_start, size) in models {
        {
            let bench = ModelBench::<f32>::new(name, path, bounds_x, bounds_y, z_start, size);
            bench.report_hit_rate();
            bench.run_init_bench(c, "f32");
            bench.run_intersect_bench(c, "f32");
        }

        {
            let bench = ModelBench::<f64>::new(name, path, bounds_x, bounds_y, z_start, size);
            bench.report_hit_rate();
            bench.run_init_bench(c, "f64");
            bench.run_intersect_bench(c, "f64");
        }
    }
}

criterion_group!(benches, bench_models);
criterion_main!(benches);

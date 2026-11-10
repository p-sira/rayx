use criterion::{Criterion, criterion_group, criterion_main};

mod common;
use common::{MODELS, ModelBench};

fn bench_models(c: &mut Criterion) {
    for (name, path, bounds_x, bounds_y, z_start, size) in MODELS {
        {
            let bench = ModelBench::<f32>::new(name, path, bounds_x, bounds_y, z_start, size);
            bench.report_hit_rate();
            bench.run_init_bench(c, "f32");
            bench.run_intersect_bench(c, "f32");
            if name == "thai_statue_50" {
                bench.run_chunked_bench(c, "f32");
            }
        }

        {
            let bench = ModelBench::<f64>::new(name, path, bounds_x, bounds_y, z_start, size);
            bench.report_hit_rate();
            bench.run_init_bench(c, "f64");
            bench.run_intersect_bench(c, "f64");
            if name == "thai_statue_50" {
                bench.run_chunked_bench(c, "f64");
            }
        }
    }
}

criterion_group!(benches, bench_models);
criterion_main!(benches);

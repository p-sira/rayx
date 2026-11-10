//! Local release-mode wall-clock measurements for the same mesh workloads as
//! `benches/intersection.rs`. Run with:
//! `taskset -c 2 cargo run --release --locked --example local_intersection`

use std::{hint::black_box, time::Instant};

use nalgebra::RealField;
use num_traits::Float;

#[allow(dead_code)]
#[path = "../benches/common.rs"]
mod common;

use common::{MODELS, ModelBench};

// Nine independently timed batches, with each batch calibrated to roughly
// 50 ms. Report the median wall time per complete workload invocation.
fn median_ms<F, R>(mut f: F) -> f64
where
    F: FnMut() -> R,
{
    for _ in 0..3 {
        black_box(f());
    }
    let start = Instant::now();
    black_box(f());
    let elapsed_ns = start.elapsed().as_nanos().max(1);
    let iterations = (50_000_000 / elapsed_ns).clamp(1, 100_000) as usize;

    let mut samples = [0.0; 9];
    for sample in &mut samples {
        let start = Instant::now();
        for _ in 0..iterations {
            black_box(f());
        }
        *sample = start.elapsed().as_secs_f64() * 1_000.0 / iterations as f64;
    }
    samples.sort_by(f64::total_cmp);
    samples[4]
}

fn benchmark_model<T: Float + RealField + Copy>(
    precision: &str,
    name: &str,
    path: &str,
    bounds_x: (f32, f32),
    bounds_y: (f32, f32),
    z_start: f32,
    size: Option<usize>,
) {
    let model = ModelBench::<T>::new(name, path, bounds_x, bounds_y, z_start, size);
    let init_ms = median_ms(|| model.initialize_observed());

    for (workload, rays) in [
        ("axis_64", &model.rays),
        ("oblique_16", &model.oblique_rays),
        ("miss_16", &model.miss_rays),
    ] {
        let bw_hits = model.count_baldwin_weber(rays);
        let mt_hits = model.count_moller_trumbore(rays);
        let bw_ms = median_ms(|| model.count_baldwin_weber(rays));
        let mt_ms = median_ms(|| model.count_moller_trumbore(rays));
        println!(
            "{name},{precision},{workload},{},{bw_hits},{mt_hits},{init_ms:.6},{bw_ms:.6},{mt_ms:.6}",
            rays.len()
        );
    }
}

fn main() {
    println!(
        "model,precision,workload,rays,bw_hits,mt_hits,bw_init_ms,bw_intersect_ms,mt_intersect_ms"
    );
    for (name, path, bounds_x, bounds_y, z_start, size) in MODELS {
        benchmark_model::<f32>("f32", name, path, bounds_x, bounds_y, z_start, size);
        benchmark_model::<f64>("f64", name, path, bounds_x, bounds_y, z_start, size);
    }
}

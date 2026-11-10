//! Compare ray-major traversal with the cache-friendly chunked batch API.
//! Run with `taskset -c 2 cargo run --release --locked --example local_chunked`.

use std::{hint::black_box, time::Instant};

use nalgebra::RealField;
use num_traits::Float;
use rayx::{Hit, intersect::intersect_baldwin_weber_chunked};

#[allow(dead_code)]
#[path = "../benches/common.rs"]
mod common;
use common::{MODELS, ModelBench};

fn median_ms(mut f: impl FnMut() -> usize) -> f64 {
    for _ in 0..3 {
        black_box(f());
    }
    let start = Instant::now();
    black_box(f());
    let count = (50_000_000 / start.elapsed().as_nanos().max(1)).clamp(1, 1_000) as usize;
    let mut samples = [0.0; 9];
    for sample in &mut samples {
        let start = Instant::now();
        for _ in 0..count {
            black_box(f());
        }
        *sample = start.elapsed().as_secs_f64() * 1_000.0 / count as f64;
    }
    samples.sort_by(f64::total_cmp);
    samples[4]
}

fn measure<T: Float + RealField + Copy>(precision: &str, bench: &ModelBench<T>) {
    for (workload, rays) in [
        ("axis_64", &bench.rays[..]),
        ("miss_16", &bench.miss_rays[..]),
        ("single_miss", &bench.miss_rays[..1]),
    ] {
        let expected = bench.count_baldwin_weber(rays);
        let baseline = median_ms(|| bench.count_baldwin_weber(rays));
        let mut scratch = vec![None; 512];
        let mut results: Vec<Option<(usize, Hit<T>)>> = vec![None; rays.len()];
        let mut batch = || {
            intersect_baldwin_weber_chunked(
                &bench.triangles_vertices,
                rays,
                &mut scratch,
                &mut results,
                T::zero(),
                T::from(1000.0).unwrap(),
            )
            .unwrap();
            results.iter().filter(|hit| hit.is_some()).count()
        };
        assert_eq!(batch(), expected);
        let chunked = median_ms(&mut batch);
        println!(
            "{},{precision},{workload},{baseline:.6},{chunked:.6}",
            bench.name
        );
    }
}

fn main() {
    println!("model,precision,workload,ray_major_ms,chunked_512_ms");
    for (name, path, bx, by, z, size) in MODELS {
        measure("f32", &ModelBench::<f32>::new(name, path, bx, by, z, size));
        measure("f64", &ModelBench::<f64>::new(name, path, bx, by, z, size));
    }
}

# RayX

`rayx` is a small Rust ray intersection library using **Baldwin-Weber's Fast Ray-Triangle Intersections by Coordinate Transformation** ([JCGT 2016](https://jcgt.org/published/0005/03/03/)).

<p>
    <a href="https://opensource.org/license/BSD-3-clause" style="text-decoration:none">
        <img src="https://img.shields.io/badge/License-BSD--3--Clause-brightgreen.svg" alt="License">
    </a>
    <a href="https://crates.io/crates/rayx" style="text-decoration:none">
        <img src="https://img.shields.io/crates/v/rayx" alt="Crate">
    </a>
    <a href="https://docs.rs/rayx" style="text-decoration:none">
        <img src="https://img.shields.io/badge/Docs-docs.rs-blue" alt="Documentation">
    </a>
    <a href="https://app.codspeed.io/p-sira/rayx?utm_source=badge" style="text-decoration:none">
        <img src="https://img.shields.io/endpoint?url=https://codspeed.io/badge.json" alt="CodSpeed"/>
    </a>
</p>

The precompute costs a little memory per triangle but reduces per-ray arithmetic.

## Usage

Install `rayx` with

```bash
cargo add rayx
```

Example:

```rust
use rayx::{Ray, Triangle};

let tri = Triangle::<f64>::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();

let ray = Ray::new([0.25, 0.25, 1.0], [0.0, 0.0, -1.0]);
let hit = tri.intersect(ray, 0.0, 10.0).unwrap();
```

## Performance

Local release benchmarks on an AMD Ryzen 5 4600H (rustc 1.99.1), pinned to one core. Each mesh uses **64 downward**, **16 oblique**, and **16 miss** rays; both algorithms scan triangles until the first hit. Values are medians of five runs.

| Model (triangles) | Type | BW init (µs) | Downward BW / MT (ms) | Downward | Oblique | Miss |
|---|---|---:|---:|---:|---:|---:|
| Suzanne (698) | f32 | 10.0 | 0.291 / 0.303 | +4.1% | +61.6% | +32.5% |
| Suzanne (698) | f64 | 10.6 | 0.284 / 0.334 | +17.5% | +71.9% | +39.5% |
| Bunny Low (908) | f32 | 16.5 | 0.445 / 0.486 | +9.2% | +55.8% | +29.2% |
| Bunny Low (908) | f64 | 17.5 | 0.446 / 0.536 | +20.2% | +81.0% | +32.7% |
| Bunny Medium (3,851) | f32 | 67.3 | 1.991 / 2.103 | +5.6% | +21.8% | +43.5% |
| Bunny Medium (3,851) | f64 | 75.1 | 2.122 / 2.223 | +4.8% | +4.0% | +87.3% |
| Bunny High (16,214) | f32 | 290.7 | 6.899 / 7.812 | +13.2% | +9.3% | +6.1% |
| Bunny High (16,214) | f64 | 339.7 | 7.404 / 7.602 | +2.7% | +2.9% | +21.2% |
| Thai Statue wide (100,000) | f32 | 2195.1 | 76.877 / 72.048 | -6.3% | -6.2% | -8.0% |
| Thai Statue wide (100,000) | f64 | 2706.0 | 87.658 / 83.513 | -4.7% | -13.9% | -11.2% |
| Thai Statue narrow (100,000) | f32 | 2255.3 | 65.705 / 61.806 | -5.9% | -6.8% | -10.1% |
| Thai Statue narrow (100,000) | f64 | 2496.4 | 77.373 / 70.611 | -8.7% | -8.1% | -9.5% |

**BW advantage** = `(MT / BW − 1) × 100%`; positive values favor Baldwin–Weber. Initialization time is separate from intersection time. The wide/narrow statue cases use different ray bounds on the same mesh.

The Bunny High `f32` downward workload has one differing hit (BW 33, MT 32), which affects its comparison. Other workloads have matching hit counts. Timings are local wall-clock results, not CodSpeed simulation estimates.

The table above measures the original ray-major traversal. For 64 downward rays on Thai Statue, `intersect::intersect_baldwin_weber_chunked` with 512-triangle chunks was **1.57× / 1.55× faster** (`f32` / `f64`) with wide bounds and **1.35× / 1.78× faster** with narrow bounds. These are medians of three separate local release runs; single-ray misses were slower with chunking.

To reproduce one run (repeat five times and take the median):

```bash
taskset -c 2 cargo run --release --locked --example local_intersection
```

For the chunked comparison, run `taskset -c 2 cargo run --release --locked --example local_chunked`.

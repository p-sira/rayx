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
</p>

The precompute costs a little memory per triangle but reduces per-ray arithmetic.

## Performance

On a benchmark at **f32** precision on a single core [1] (2,500 rays vs meshes):

| Model | Triangles | Mem (BW/MT; KB) | Hit | BW Init | BW Intersect | MT Intersect | Speed Up |
|-------|-----------|-----|-----|---------|--------------|--------------|----------|
| Suzanne | 698 | 33.5/25.1 | 43.5% | 1.7 µs | 7.0 ms | 10.2 ms | 45.7% |
| Bunny (Low) | 908 | 43.6/32.7 | 51.0% | 2.3 µs | 9.4 ms | 13.3 ms | 41.5% |
| Bunny (Med) | 3,851 | 185/139 | 52.5% | 9.6 µs | 48.4 ms | 54.5 ms | 12.6% |
| Bunny (High) | 16,214 | 778/584 | 53.1% | 39.9 µs | 240.6 ms | 299.0 ms | 24.3% |
| Thai Statue | 100,000 | 4800/3600 | 53.7% | 260 µs | 2.22 s | 2.50 s | 12.6% |
| Thai Statue | 100,000 | 4800/3600 | 78.6% | 286 µs | 1.77 s | 1.97 s | 11.3% |

On a benchmark at **f64** precision:

| Model | Triangles | Mem (BW/MT; KB) | Hit | BW Init | BW Intersect | MT Intersect | Speed Up |
|-------|-----------|-----|-----|---------|--------------|--------------|----------|
| Suzanne | 698 | 67.0/50.2 | 43.5% | 1.5 µs | 8.2 ms | 9.7 ms | 18.3% |
| Bunny (Low) | 908 | 87.2/65.4 | 51.0% | 2.0 µs | 10.8 ms | 11.7 ms | 8.3% |
| Bunny (Med) | 3,851 | 370/278 | 52.5% | 8.8 µs | 54.0 ms | 61.1 ms | 13.1% |
| Bunny (High) | 16,214 | 1556/1168 | 53.1% | 36.5 µs | 264.4 ms | 288.5 ms | 9.1% |
| Thai Statue | 100,000 | 9600/7200\* | 53.7% | 435 µs | 2.43 s | 2.43 s | 0% |
| Thai Statue | 100,000 | 9600/7200\* | 78.6% | 500 µs | 1.92 s | 1.92 s | 0% |

**Baldwin-Weber consistently outperforms Moller-Trumbore in intersection time across various model complexities.** While it requires a one-time precomputation phase, the cost is minimal and quickly amortized over many ray intersections. The performance gap is most pronounced at lower precisions (f32) where the reduced arithmetic complexity yields up to 45% speedup.

\* Notice the Baldwin-Weber implementation spills into the RAM, while Moller-Trumbore implementation stays within the L3 cache. This is due to the larger memory footprint of the Baldwin-Weber implementation.

[1]: AMD Ryzen 5 4600H with Radeon Graphics @3.0 GHz running x86_64-unknown-linux-gnu with rustc 1.90.0, L2 cache 3 MiB, L3 cache 8 MiB, RAM 16 GB, SSD SATA 600 (6 Gbps).

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

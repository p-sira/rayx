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

## Usage

Install `rayx` with

```bash
cargo add rayx
```

Example:

```rust
use nalgebra::Vector3;
use rayx::{Ray, Triangle};

let tri = Triangle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();

let ray = Ray::new([0.25, 0.25, 1.0], [0.0, 0.0, -1.0]);
let hit = tri.intersect(ray, 0.0, 10.0).unwrap();

assert!((hit.t - 1.0).abs() < 1e-12);
assert!((hit.u - 0.25).abs() < 1e-12);
assert!((hit.v - 0.25).abs() < 1e-12);
```

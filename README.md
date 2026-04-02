# RayX

`rayx` is a small Rust ray intersection library using **Baldwin-Weber's Fast Ray-Triangle Intersections by Coordinate Transformation** ([JCGT 2016] (https://jcgt.org/published/0005/03/03/)).

## Usage

Install `rayx` with

```bash
cargo add rayx
```

Example:

```rust
use nalgebra::Vector3;
use rayx::{Ray, Triangle};

fn main() {
    let tri = Triangle::new(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    )
    .unwrap();

    let ray = Ray::new(Vector3::new(0.25, 0.25, 1.0), Vector3::new(0.0, 0.0, -1.0));
    let hit = tri.intersect(ray, 0.0, 10.0).unwrap();

    assert!((hit.t - 1.0).abs() < 1e-6);
    assert!((hit.b1 - 0.25).abs() < 1e-6);
    assert!((hit.b2 - 0.25).abs() < 1e-6);
}
```

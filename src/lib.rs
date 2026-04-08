//! RayX: a small, fast ray intersection library.
//!
//! `rayx` is a small Rust ray intersection library using **Baldwin-Weber's Fast Ray-Triangle Intersections by Coordinate Transformation** ([JCGT 2016](https://jcgt.org/published/0005/03/03/)).
//!
//! The precompute costs a little memory per triangle but reduces per-ray arithmetic.
//!
//! ## Performance
//!
//! On a benchmark at **f32** precision on a single core [1] (2,500 rays vs meshes):
//!
//! | Model | Triangles | Mem (BW/MT; KB) | Hit | BW Init | BW Intersect | MT Intersect | Speed Up |
//! |-------|-----------|-----|-----|---------|--------------|--------------|----------|
//! | Suzanne | 698 | 33.5/25.1 | 43.5% | 1.7 µs | 7.0 ms | 10.2 ms | 45.7% |
//! | Bunny (Low) | 908 | 43.6/32.7 | 51.0% | 2.3 µs | 9.4 ms | 13.3 ms | 41.5% |
//! | Bunny (Med) | 3,851 | 185/139 | 52.5% | 9.6 µs | 48.4 ms | 54.5 ms | 12.6% |
//! | Bunny (High) | 16,214 | 778/584 | 53.1% | 39.9 µs | 240.6 ms | 299.0 ms | 24.3% |
//! | Thai Statue | 100,000 | 4800/3600 | 53.7% | 260 µs | 2.22 s | 2.50 s | 12.6% |
//! | Thai Statue | 100,000 | 4800/3600 | 78.6% | 286 µs | 1.77 s | 1.97 s | 11.3% |
//!
//! On a benchmark at **f64** precision:
//!
//! | Model | Triangles | Mem (BW/MT; KB) | Hit | BW Init | BW Intersect | MT Intersect | Speed Up |
//! |-------|-----------|-----|-----|---------|--------------|--------------|----------|
//! | Suzanne | 698 | 67.0/50.2 | 43.5% | 1.5 µs | 8.2 ms | 9.7 ms | 18.3% |
//! | Bunny (Low) | 908 | 87.2/65.4 | 51.0% | 2.0 µs | 10.8 ms | 11.7 ms | 8.3% |
//! | Bunny (Med) | 3,851 | 370/278 | 52.5% | 8.8 µs | 54.0 ms | 61.1 ms | 13.1% |
//! | Bunny (High) | 16,214 | 1556/1168 | 53.1% | 36.5 µs | 264.4 ms | 288.5 ms | 9.1% |
//! | Thai Statue | 100,000 | 9600/7200\* | 53.7% | 435 µs | 2.43 s | 2.43 s | 0% |
//! | Thai Statue | 100,000 | 9600/7200\* | 78.6% | 500 µs | 1.92 s | 1.92 s | 0% |
//!
//! **Baldwin-Weber consistently outperforms Moller-Trumbore in intersection time across various model complexities.**
//! While it requires a one-time precomputation phase, the cost is minimal and quickly amortized over many ray intersections.
//! The performance gap is most pronounced at lower precisions (f32) where the reduced arithmetic complexity yields up to 45% speedup.
//!
//! \* Notice the Baldwin-Weber implementation spills into the RAM, while Moller-Trumbore implementation stays within the L3 cache. This is due to the larger memory footprint of the Baldwin-Weber implementation.
//!
//! [1]: AMD Ryzen 5 4600H with Radeon Graphics @3.0 GHz running x86_64-unknown-linux-gnu with rustc 1.90.0,
//! RAM 16 GB, L2 cache 3 MiB, L3 cache 8 MiB, SSD SATA 600 (6 Gbps).
//!
//! ## Usage
//!
//! ```
//! use rayx::{Ray, Triangle};
//!
//! let tri = Triangle::<f64>::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
//!
//! let ray = Ray::new([0.25, 0.25, 1.0], [0.0, 0.0, -1.0]);
//! let hit = tri.intersect(ray, 0.0, 10.0).unwrap();
//! # assert!((hit.t - 1.0).abs() < 1e-12);
//! # assert!((hit.u - 0.25).abs() < 1e-12);
//! # assert!((hit.v - 0.25).abs() < 1e-12);
//! ```

#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(docsrs, feature(doc_cfg))]

use core::fmt;

mod base;
pub mod intersect;

pub use base::{Hit, Ray};

use nalgebra::{RealField, Vector3};
use num_traits::Float;

/// Triangle with precomputed transform coefficients for fast intersection tests.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Triangle<T: Float + RealField = f32> {
    m: nalgebra::Matrix3<T>,
    v: nalgebra::Vector3<T>,
}

impl<T: Float + RealField> Triangle<T> {
    /// Build a triangle and store the precomputed transform.
    ///
    /// Note: The vertices are not stored.
    #[inline]
    pub fn new(
        v1: impl Into<Vector3<T>>,
        v2: impl Into<Vector3<T>>,
        v3: impl Into<Vector3<T>>,
    ) -> Result<Self, TriangleError> {
        let m = intersect::compute_baldwin_weber_transform(v1.into(), v2.into(), v3.into())?;
        Ok(Self {
            m: m.fixed_view::<3, 3>(0, 0).into(),
            v: m.column(3).into(),
        })
    }

    /// Intersect ray with this triangle using the precomputed transform.
    ///
    /// Returns `Some(Hit)` if the ray intersects within `[t_min, t_max]`.
    #[inline]
    pub fn intersect(&self, ray: Ray<T>, t_min: T, t_max: T) -> Option<Hit<T>> {
        intersect::intersect_baldwin_weber(&self.m, &self.v, ray, t_min, t_max)
    }

    /// Reconstruct the original triangle vertices from the precomputed transform.
    #[inline]
    pub fn reconstruct_vertices(&self) -> [Vector3<T>; 3] {
        let inv_m = self
            .m
            .try_inverse()
            .expect("Cannot invert the triangle matrix");

        let v1 = inv_m * (-self.v);
        let v2 = v1 + inv_m.column(0);
        let v3 = v1 + inv_m.column(1);

        [v1, v2, v3]
    }

    #[inline]
    pub fn m(&self) -> nalgebra::Matrix3<T> {
        self.m
    }

    #[inline]
    pub fn v(&self) -> nalgebra::Vector3<T> {
        self.v
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TriangleError {
    Degenerate,
}

impl fmt::Display for TriangleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TriangleError::Degenerate => write!(f, "Triangle is degenerate (zero surface area)"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for TriangleError {}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx<T: Float>(a: T, b: T, tol: T) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn degenerate_triangle_rejected() {
        let zeros = Vector3::new(0.0, 0.0, 0.0);
        assert_eq!(
            Triangle::new(zeros, zeros, zeros),
            Err(TriangleError::Degenerate)
        );
    }

    #[test]
    fn basic_hit_matches_expected_barycentrics() {
        // Canonical triangle in z=0 plane: (0,0,0), (1,0,0), (0,1,0)
        let tri = Triangle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();

        let ray = Ray::new([0.25, 0.25, 1.0], [0.0, 0.0, -1.0]);
        let hit = tri.intersect(ray, 0.0, 10.0).unwrap();
        assert!(approx(hit.t, 1.0, 1e-12));
        assert!(approx(hit.u, 0.25, 1e-12));
        assert!(approx(hit.v, 0.25, 1e-12));
        assert!(approx(hit.w(), 0.5, 1e-12));
    }

    #[test]
    fn miss_outside_triangle() {
        let tri = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        )
        .unwrap();

        let ray = Ray::new(Vector3::new(0.8, 0.8, 1.0), Vector3::new(0.0, 0.0, -1.0));
        assert_eq!(tri.intersect(ray, 0.0, 10.0), None);
    }

    #[test]
    fn agrees_with_moller_trumbore() {
        let v0 = Vector3::new(-1.0, 0.0, 0.5);
        let v1 = Vector3::new(2.0, 0.3, -0.2);
        let v2 = Vector3::new(0.2, 1.7, 0.1);
        let tri = Triangle::new(v0, v1, v2).unwrap();

        let rays = [
            Ray::new(Vector3::new(0.1, 0.2, 2.0), Vector3::new(0.0, 0.0, -1.0)),
            Ray::new(
                Vector3::new(10.0, 10.0, 10.0),
                Vector3::new(-1.0, -1.0, -1.0),
            ),
            Ray::new(Vector3::new(0.0, -2.0, 0.0), Vector3::new(0.2, 1.0, 0.1)),
            Ray::new(Vector3::new(0.3, 0.4, 0.5), Vector3::new(1.0, 0.1, -0.2)),
        ];

        for ray in rays {
            let a = tri.intersect(ray, 0.0, 1e9);
            let b = intersect::intersect_moller_trumbore(v0, v1, v2, ray, 0.0, 1e9);
            match (a, b) {
                (None, None) => {}
                (Some(ha), Some(hb)) => {
                    assert!(approx(ha.t, hb.t, 1e-12), "t mismatch: {ha:?} vs {hb:?}");
                    assert!(approx(ha.u, hb.u, 1e-12), "b1 mismatch: {ha:?} vs {hb:?}");
                    assert!(approx(ha.v, hb.v, 1e-12), "b2 mismatch: {ha:?} vs {hb:?}");
                }
                _ => panic!("disagreement: precomputed={a:?}, mt={b:?}"),
            }
        }
    }

    #[test]
    fn reconstruct_vertices_matches_input() {
        let v0 = Vector3::new(-1.2, 0.5, 3.1);
        let v1 = Vector3::new(2.4, -1.1, 0.8);
        let v2 = Vector3::new(0.3, 4.2, -1.5);

        let tri = Triangle::new(v0, v1, v2).unwrap();
        let reconstructed = tri.reconstruct_vertices();

        assert!(approx(reconstructed[0].x, v0.x, 1e-12));
        assert!(approx(reconstructed[0].y, v0.y, 1e-12));
        assert!(approx(reconstructed[0].z, v0.z, 1e-12));

        assert!(approx(reconstructed[1].x, v1.x, 1e-12));
        assert!(approx(reconstructed[1].y, v1.y, 1e-12));
        assert!(approx(reconstructed[1].z, v1.z, 1e-12));

        assert!(approx(reconstructed[2].x, v2.x, 1e-12));
        assert!(approx(reconstructed[2].y, v2.y, 1e-12));
        assert!(approx(reconstructed[2].z, v2.z, 1e-12));
    }
}

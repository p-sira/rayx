//! RayX: a small, fast ray intersection library.
//!
//! `rayx` is a small Rust ray intersection library using **Baldwin-Weber's Fast Ray-Triangle Intersections by Coordinate Transformation** ([JCGT 2016] (https://jcgt.org/published/0005/03/03/)).
//!
//! The precompute costs a little memory per triangle but reduces per-ray arithmetic.

#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(docsrs, feature(doc_cfg))]

use core::fmt;

mod base;
pub mod intersect;

pub use base::{Hit, Ray};

use nalgebra::{RealField, Vector3};
use num_traits::Float;

/// Triangle with precomputed transform coefficients for fast intersection tests.
///
/// Stores the top 3 rows of the global→barycentric transform matrix (12 coefficients).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Triangle<T: Float + RealField> {
    // Row-major 3x4 matrix, applied as:
    // [x'; y'; z'] = M * [x y z 1]^T
    pub m: [T; 12],
}

impl<T: Float + RealField> Triangle<T> {
    /// Build a triangle and store the precomputed transform.
    ///
    /// Note: The vertices are not stored.
    pub fn new(v1: Vector3<T>, v2: Vector3<T>, v3: Vector3<T>) -> Result<Self, TriangleError> {
        let m = intersect::compute_baldwin_weber_transform(v1, v2, v3)?;
        Ok(Self { m })
    }

    /// Intersect ray with this triangle using the precomputed transform.
    ///
    /// Returns `Some(Hit)` if the ray intersects within `[t_min, t_max]`.
    pub fn intersect(&self, ray: Ray<T>, t_min: T, t_max: T) -> Option<Hit<T>> {
        intersect::intersect_baldwin_weber(&self.m, ray, t_min, t_max)
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
        let tri = Triangle::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        )
        .unwrap();

        let ray = Ray::new(Vector3::new(0.25, 0.25, 1.0), Vector3::new(0.0, 0.0, -1.0));
        let hit = tri.intersect(ray, 0.0, 10.0).unwrap();
        assert!(approx(hit.t, 1.0, 1e-6));
        assert!(approx(hit.b1, 0.25, 1e-6));
        assert!(approx(hit.b2, 0.25, 1e-6));
        assert!(approx(hit.b0(), 0.5, 1e-6));
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
                    assert!(approx(ha.t, hb.t, 1e-4), "t mismatch: {ha:?} vs {hb:?}");
                    assert!(approx(ha.b1, hb.b1, 1e-3), "b1 mismatch: {ha:?} vs {hb:?}");
                    assert!(approx(ha.b2, hb.b2, 1e-3), "b2 mismatch: {ha:?} vs {hb:?}");
                }
                _ => panic!("disagreement: precomputed={a:?}, mt={b:?}"),
            }
        }
    }
}

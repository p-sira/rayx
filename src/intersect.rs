use nalgebra::{RealField, Vector3};
use num_traits::Float;

use crate::{Hit, Ray, TriangleError};

/// Computes the Baldwin-Weber transform matrix (top 3 rows) for a triangle.
/// Returns the 12 matrix coefficients.
#[inline]
pub fn compute_baldwin_weber_transform<T: Float + RealField>(
    v1: Vector3<T>,
    v2: Vector3<T>,
    v3: Vector3<T>,
) -> Result<[T; 12], TriangleError> {
    let e1 = v2 - v1;
    let e2 = v3 - v1;
    let n = e1.cross(&e2);

    if n.dot(&n) < T::epsilon() {
        return Err(TriangleError::Degenerate);
    }

    let ax = Float::abs(n.x);
    let ay = Float::abs(n.y);
    let az = Float::abs(n.z);

    // Paper defines three cases for inverse of T based on the major axis
    let m = if ax >= ay && ax >= az {
        let inv = T::one() / n.x;
        let v3xv1 = v3.cross(&v1).x;
        let v2xv1 = v2.cross(&v1).x;
        [
            T::zero(),
            e2.z * inv,
            (-e2.y) * inv,
            v3xv1 * inv,
            T::zero(),
            (-e1.z) * inv,
            e1.y * inv,
            (-v2xv1) * inv,
            T::one(),
            n.y * inv,
            n.z * inv,
            (-n.dot(&v1)) * inv,
        ]
    } else if ay >= az {
        let inv = T::one() / n.y;
        let v3xv1 = v3.cross(&v1).y;
        let v2xv1 = v2.cross(&v1).y;
        [
            (-e2.z) * inv,
            T::zero(),
            e2.x * inv,
            v3xv1 * inv,
            e1.z * inv,
            T::zero(),
            (-e1.x) * inv,
            (-v2xv1) * inv,
            n.x * inv,
            T::one(),
            n.z * inv,
            (-n.dot(&v1)) * inv,
        ]
    } else {
        let inv = T::one() / n.z;
        let v3xv1 = v3.cross(&v1).z;
        let v2xv1 = v2.cross(&v1).z;
        [
            e2.y * inv,
            (-e2.x) * inv,
            T::zero(),
            v3xv1 * inv,
            (-e1.y) * inv,
            e1.x * inv,
            T::zero(),
            (-v2xv1) * inv,
            n.x * inv,
            n.y * inv,
            T::one(),
            (-n.dot(&v1)) * inv,
        ]
    };

    Ok(m)
}

/// Baldwin-Weber ray–triangle intersection algorithm.
#[inline]
pub fn intersect_baldwin_weber<T: Float + RealField>(
    m: &[T; 12],
    ray: Ray<T>,
    t_min: T,
    t_max: T,
) -> Option<Hit<T>> {
    // Transform origin as point (w=1) and direction as vector (w=0).
    let ox = m[0] * ray.origin.x + m[1] * ray.origin.y + m[2] * ray.origin.z + m[3];
    let oy = m[4] * ray.origin.x + m[5] * ray.origin.y + m[6] * ray.origin.z + m[7];
    let oz = m[8] * ray.origin.x + m[9] * ray.origin.y + m[10] * ray.origin.z + m[11];

    let dx = m[0] * ray.dir.x + m[1] * ray.dir.y + m[2] * ray.dir.z;
    let dy = m[4] * ray.dir.x + m[5] * ray.dir.y + m[6] * ray.dir.z;
    let dz = m[8] * ray.dir.x + m[9] * ray.dir.y + m[10] * ray.dir.z;

    // Epsilon check prevents division by extremely small numbers or exact zeroes
    // yielding unreliable infinities.
    if Float::abs(dz) < T::epsilon() {
        return None;
    }

    let t = -oz / dz;
    if !(t_min..=t_max).contains(&t) {
        return None;
    }

    let b1 = ox + t * dx;
    if b1 < T::zero() || b1 > T::one() {
        return None;
    }

    let b2 = oy + t * dy;
    if b2 < T::zero() || b1 + b2 > T::one() {
        return None;
    }

    Some(Hit { t, u: b1, v: b2 })
}

/// Möller–Trumbore ray–triangle intersection algorithm.
pub fn intersect_moller_trumbore<T: Float + RealField>(
    v1: Vector3<T>,
    v2: Vector3<T>,
    v3: Vector3<T>,
    ray: Ray<T>,
    t_min: T,
    t_max: T,
) -> Option<Hit<T>> {
    let eps = T::epsilon() * T::from(16.0).unwrap();
    let e1 = v2 - v1;
    let e2 = v3 - v1;
    let p = ray.dir.cross(&e2);
    let det = e1.dot(&p);
    if Float::abs(det) < eps {
        return None;
    }
    let inv_det = T::one() / det;
    let tvec = ray.origin - v1;
    let u = tvec.dot(&p) * inv_det;
    if u < T::zero() || u > T::one() {
        return None;
    }
    let q = tvec.cross(&e1);
    let v = ray.dir.dot(&q) * inv_det;
    if v < T::zero() || u + v > T::one() {
        return None;
    }
    let t = e2.dot(&q) * inv_det;
    if !(t_min..=t_max).contains(&t) {
        return None;
    }
    Some(Hit { t, u, v })
}

use nalgebra::{Matrix3x4, RealField, Vector3};
use num_traits::Float;

use crate::{Hit, Ray, TriangleError};

/// Computes the Baldwin-Weber transform matrix (top 3 rows) for a triangle.
/// Returns the 12 matrix coefficients.
#[inline]
pub fn compute_baldwin_weber_transform<T: Float + RealField>(
    v1: Vector3<T>,
    v2: Vector3<T>,
    v3: Vector3<T>,
) -> Result<Matrix3x4<T>, TriangleError> {
    let e1 = v2 - v1;
    let e2 = v3 - v1;
    let n = e1.cross(&e2);

    if n.dot(&n) < T::epsilon() {
        return Err(TriangleError::Degenerate);
    }

    let a = n.abs();

    // Paper defines three cases for inverse of T based on the major axis
    let m = if a.x >= a.y && a.x >= a.z {
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
    } else if a.y >= a.z {
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

    Ok(Matrix3x4::from_row_slice(&m))
}

/// Baldwin-Weber ray–triangle intersection algorithm.
#[inline]
pub fn intersect_baldwin_weber<T: Float + RealField>(
    m: &Matrix3x4<T>,
    ray: Ray<T>,
    t_min: T,
    t_max: T,
) -> Option<Hit<T>> {
    let d = m.fixed_view::<3, 3>(0, 0) * ray.dir;

    // Epsilon check prevents division by extremely small numbers or exact zeroes
    // yielding unreliable infinities.
    if Float::abs(d.z) < T::epsilon() {
        return None;
    }

    let o = m.fixed_view::<3, 3>(0, 0) * ray.origin + m.column(3);

    let t = -o.z / d.z;
    if !(t_min..=t_max).contains(&t) {
        return None;
    }

    let b1 = o.x + t * d.x;
    if b1 < T::zero() || b1 > T::one() {
        return None;
    }

    let b2 = o.y + t * d.y;
    if b2 < T::zero() || b1 + b2 > T::one() {
        return None;
    }

    Some(Hit { t, u: b1, v: b2 })
}

/// Möller–Trumbore ray–triangle intersection algorithm.
#[inline]
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

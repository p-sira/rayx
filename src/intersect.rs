use nalgebra::{Matrix3x4, RealField, Vector3};
use num_traits::Float;

use crate::{Hit, Ray, Triangle, TriangleError};

/// Invalid buffers supplied to [`intersect_baldwin_weber_chunked`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ChunkedIntersectError {
    /// Each ray needs one result slot.
    ResultLengthMismatch,
    /// A nonempty ray batch and mesh require at least one scratch slot.
    EmptyScratch,
}

impl core::fmt::Display for ChunkedIntersectError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ResultLengthMismatch => f.write_str("one result slot is required for each ray"),
            Self::EmptyScratch => f.write_str("triangle scratch buffer cannot be empty"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for ChunkedIntersectError {}

/// Intersect a batch of rays against a mesh, precomputing only one chunk at a time.
///
/// Process each triangle chunk against all rays before moving to the next chunk,
/// so its Baldwin–Weber transforms can be reused while cache-resident. The
/// caller supplies reusable buffers; this function does not allocate and also
/// works in `no_std` environments.
///
/// Each result contains the **first hit in input triangle order**, together
/// with its triangle index. It is not necessarily the closest hit. Degenerate
/// triangles are skipped. Results are reset on every successful call.
///
/// A scratch buffer of 2,048 slots is a useful starting point for benchmarking,
/// but the best size depends on the CPU and batch size.
///
/// ```
/// use rayx::{Ray, intersect::intersect_baldwin_weber_chunked};
/// let triangles = [[
///     [0.0, 0.0, 0.0].into(),
///     [1.0, 0.0, 0.0].into(),
///     [0.0, 1.0, 0.0].into(),
/// ]];
/// let rays = [Ray::new([0.25, 0.25, 1.0], [0.0, 0.0, -1.0])];
/// let mut scratch = [None; 2];
/// let mut hits = [None; 1];
/// intersect_baldwin_weber_chunked(&triangles, &rays, &mut scratch, &mut hits, 0.0, 10.0)
///     .unwrap();
/// assert_eq!(hits[0].unwrap().0, 0);
/// ```
pub fn intersect_baldwin_weber_chunked<T: Float + RealField>(
    vertices: &[[Vector3<T>; 3]],
    rays: &[Ray<T>],
    scratch: &mut [Option<Triangle<T>>],
    results: &mut [Option<(usize, Hit<T>)>],
    t_min: T,
    t_max: T,
) -> Result<(), ChunkedIntersectError> {
    if results.len() != rays.len() {
        return Err(ChunkedIntersectError::ResultLengthMismatch);
    }
    if scratch.is_empty() && !vertices.is_empty() && !rays.is_empty() {
        return Err(ChunkedIntersectError::EmptyScratch);
    }

    results.fill(None);
    if vertices.is_empty() || rays.is_empty() {
        return Ok(());
    }

    let mut remaining = rays.len();
    for (chunk_index, chunk) in vertices.chunks(scratch.len()).enumerate() {
        for (slot, v) in scratch.iter_mut().zip(chunk) {
            *slot = Triangle::new(v[0], v[1], v[2]).ok();
        }

        for (ray, result) in rays.iter().zip(results.iter_mut()) {
            if result.is_some() {
                continue;
            }
            for (index, triangle) in scratch[..chunk.len()].iter().enumerate() {
                if let Some(hit) = triangle
                    .as_ref()
                    .and_then(|tri| tri.intersect(*ray, t_min, t_max))
                {
                    *result = Some((chunk_index * scratch.len() + index, hit));
                    remaining -= 1;
                    break;
                }
            }
        }

        if remaining == 0 {
            break;
        }
    }
    Ok(())
}

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

    // Threshold for degeneracy should be small enough to allow small valid triangles,
    // especially in f32 where epsilon() is relatively large (1e-7).
    let area_eps = T::from(1e-20).unwrap_or_else(T::epsilon);
    if n.dot(&n) < area_eps {
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
    m: &nalgebra::Matrix3<T>,
    v: &nalgebra::Vector3<T>,
    ray: Ray<T>,
    t_min: T,
    t_max: T,
) -> Option<Hit<T>> {
    let d = m * ray.dir;

    // Epsilon check prevents division by extremely small numbers or exact zeroes
    // yielding unreliable infinities.
    if Float::abs(d.z) < T::epsilon() {
        return None;
    }

    let o = m * ray.origin + v;

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

#[cfg(test)]
mod chunked_tests {
    use super::*;

    fn check<T: Float + RealField + core::fmt::Debug>() {
        let z = T::zero();
        let one = T::one();
        let half = one / (one + one);
        let vertices = [
            [
                Vector3::new(z, z, z),
                Vector3::new(one, z, z),
                Vector3::new(z, one, z),
            ],
            // The degenerate triangle should not change triangle indices.
            [Vector3::new(z, z, z); 3],
            [
                Vector3::new(z, z, one),
                Vector3::new(one, z, one),
                Vector3::new(z, one, one),
            ],
            [
                Vector3::new(z, z, one + one),
                Vector3::new(one, z, one + one),
                Vector3::new(z, one, one + one),
            ],
        ];
        let rays = [
            Ray::new(
                Vector3::new(half / (one + one), half / (one + one), one + one + one),
                Vector3::new(z, z, -one),
            ),
            Ray::new(
                Vector3::new(half / (one + one), half / (one + one), -one),
                Vector3::new(z, z, one),
            ),
            Ray::new(
                Vector3::new(one + one, one + one, one),
                Vector3::new(z, z, -one),
            ),
        ];
        let expected: [Option<(usize, Hit<T>)>; 3] = core::array::from_fn(|i| {
            vertices.iter().enumerate().find_map(|(j, v)| {
                Triangle::new(v[0], v[1], v[2]).ok().and_then(|tri| {
                    tri.intersect(rays[i], z, one + one + one + one)
                        .map(|hit| (j, hit))
                })
            })
        });
        for chunk_size in 1..=vertices.len() + 1 {
            let mut scratch = vec![None; chunk_size];
            let mut results = vec![None; rays.len()];
            intersect_baldwin_weber_chunked(
                &vertices,
                &rays,
                &mut scratch,
                &mut results,
                z,
                one + one + one + one,
            )
            .unwrap();
            assert_eq!(results, expected, "chunk size {chunk_size}");
            // The downward ray hits z=2 first geometrically, but the z=0
            // triangle appears first in the original mesh order.
            assert_eq!(results[0].map(|p| p.0), Some(0));
            assert_eq!(results[1].map(|p| p.0), Some(0));
            assert!(results[2].is_none());
        }
    }

    #[test]
    fn matches_sequential_first_hit_f32() {
        check::<f32>();
    }

    #[test]
    fn matches_sequential_first_hit_f64() {
        check::<f64>();
    }

    #[test]
    fn rejects_invalid_buffers_without_mutating_results() {
        let v = [[
            Vector3::new(0., 0., 0.),
            Vector3::new(1., 0., 0.),
            Vector3::new(0., 1., 0.),
        ]];
        let rays = [Ray::new([0.2, 0.2, 1.], [0., 0., -1.])];
        let mut no_scratch = [];
        let mut results = [Some((
            7,
            Hit {
                t: 1.,
                u: 0.,
                v: 0.,
            },
        ))];
        assert_eq!(
            intersect_baldwin_weber_chunked(&v, &rays, &mut no_scratch, &mut results, 0., 2.),
            Err(ChunkedIntersectError::EmptyScratch)
        );
        assert_eq!(results[0].unwrap().0, 7);
        let mut no_results = [];
        assert_eq!(
            intersect_baldwin_weber_chunked(&v, &rays, &mut no_scratch, &mut no_results, 0., 2.),
            Err(ChunkedIntersectError::ResultLengthMismatch)
        );
        intersect_baldwin_weber_chunked(&v, &[], &mut no_scratch, &mut no_results, 0., 2.).unwrap();
        intersect_baldwin_weber_chunked(&[], &rays, &mut no_scratch, &mut results, 0., 2.).unwrap();
        assert_eq!(results, [None]);
    }
}

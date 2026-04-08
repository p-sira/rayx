use nalgebra::{RealField, Vector3};
use num_traits::Float;

// MARK: Ray

/// A ray with parametric form `o + t*d`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ray<T: Float + RealField = f32> {
    pub origin: Vector3<T>,
    pub dir: Vector3<T>,
}

impl<T: Float + RealField> Ray<T> {
    #[inline]
    pub fn new(origin: impl Into<Vector3<T>>, dir: impl Into<Vector3<T>>) -> Self {
        Self {
            origin: origin.into(),
            dir: dir.into(),
        }
    }

    #[inline]
    pub fn at(self, t: T) -> Vector3<T> {
        self.origin + self.dir * t
    }
}

// MARK: Hit

/// Ray–triangle intersection result with the distance to intersection
/// `t` at the barycentric coordinates (u, v).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit<T: Float + RealField = f32> {
    pub t: T,
    pub u: T,
    pub v: T,
}

impl<T: Float + RealField> Hit<T> {
    #[inline]
    pub fn w(self) -> T {
        T::one() - self.u - self.v
    }
}

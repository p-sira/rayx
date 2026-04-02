use nalgebra::{RealField, Vector3};
use num_traits::Float;

// MARK: Ray

/// A ray with parametric form `o + t*d`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ray<T: Float + RealField> {
    pub origin: Vector3<T>,
    pub dir: Vector3<T>,
}

impl<T: Float + RealField> Ray<T> {
    #[inline]
    pub fn new(origin: Vector3<T>, dir: Vector3<T>) -> Self {
        Self { origin, dir }
    }

    #[inline]
    pub fn at(self, t: T) -> Vector3<T> {
        self.origin + self.dir * t
    }
}

// MARK: Hit

/// Ray–triangle intersection result (barycentric coords (b1, b2)).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit<T: Float + RealField> {
    pub t: T,
    pub b1: T,
    pub b2: T,
}

impl<T: Float + RealField> Hit<T> {
    #[inline]
    pub fn b0(self) -> T {
        T::one() - self.b1 - self.b2
    }
}

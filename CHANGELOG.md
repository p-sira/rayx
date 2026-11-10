# 0.3

## 0.3.0

**Performance**

- Add an allocation-free, `no_std`-compatible chunked Baldwin–Weber batch intersection API that reuses transforms across rays and preserves triangle-order first hits ([#3](https://github.com/p-sira/rayx/pull/3)).

**Testing**

- Add representative downward, oblique, and miss-ray benchmarks, plus reproducible local performance comparisons ([#3](https://github.com/p-sira/rayx/pull/3)).

# 0.2

## 0.2.0

- Add method to reconstruct vertices from `m` and `v`.
- Remove public access to the `Triangle`'s fields.
- Add benchmarks.

# 0.1

## 0.1.2

- Implement `Eq` for `Triangle`.

## 0.1.1

- Upgrade `nalgebra` to v0.34.
- Accept types that implements `Into<Vector3>` for structs.

## 0.1.0

- First Release

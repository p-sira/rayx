use criterion::{Criterion, black_box};
use nalgebra::{RealField, Vector3};
use num_traits::{Float, cast};
use rayx::{Ray, Triangle, intersect};
use std::fs::File;

pub type ModelSpec = (
    &'static str,
    &'static str,
    (f32, f32),
    (f32, f32),
    f32,
    Option<usize>,
);

pub const MODELS: [ModelSpec; 6] = [
    (
        "suzanne",
        "tests/test_data/perfect-suzanne.stl",
        (-1.5, 1.5),
        (-1.0, 1.0),
        5.0,
        None,
    ),
    (
        "bunny_res2",
        "tests/test_data/bun_zipper_res2.stl",
        (-0.1, 0.07),
        (0.03, 0.19),
        1.0,
        Some(30),
    ),
    (
        "bunny_res3",
        "tests/test_data/bun_zipper_res3.stl",
        (-0.1, 0.07),
        (0.03, 0.19),
        1.0,
        Some(90),
    ),
    (
        "bunny_res4",
        "tests/test_data/bun_zipper_res4.stl",
        (-0.1, 0.07),
        (0.03, 0.19),
        1.0,
        Some(100),
    ),
    (
        "thai_statue_50",
        "tests/test_data/thai-statue-decim001.stl",
        (-75.0, 75.0),
        (-100.0, 100.0),
        200.0,
        Some(10),
    ),
    (
        "thai_statue_80",
        "tests/test_data/thai-statue-decim001.stl",
        (-50.0, 50.0),
        (-100.0, 100.0),
        200.0,
        Some(10),
    ),
];

pub struct ModelBench<T: Float + RealField + Copy> {
    pub name: String,
    pub triangles_vertices: Vec<[Vector3<T>; 3]>,
    pub triangles_precomputed: Vec<Triangle<T>>,
    pub rays: Vec<Ray<T>>,
    pub oblique_rays: Vec<Ray<T>>,
    pub miss_rays: Vec<Ray<T>>,
    pub size: Option<usize>,
}

impl<T: Float + RealField + Copy> ModelBench<T> {
    pub fn new(
        name: &str,
        path: &str,
        bounds_x: (f32, f32),
        bounds_y: (f32, f32),
        z_start: f32,
        size: Option<usize>,
    ) -> Self {
        let mut file = File::open(path).expect("Failed to open model");
        let mesh = stl_io::read_stl(&mut file).expect("Failed to read model");

        let mut degenerate_count = 0usize;
        let triangles_data: Vec<([Vector3<T>; 3], Triangle<T>)> = mesh
            .faces
            .iter()
            .filter_map(|f| {
                let v = [
                    Vector3::new(
                        cast(mesh.vertices[f.vertices[0]][0]).unwrap(),
                        cast(mesh.vertices[f.vertices[0]][1]).unwrap(),
                        cast(mesh.vertices[f.vertices[0]][2]).unwrap(),
                    ),
                    Vector3::new(
                        cast(mesh.vertices[f.vertices[1]][0]).unwrap(),
                        cast(mesh.vertices[f.vertices[1]][1]).unwrap(),
                        cast(mesh.vertices[f.vertices[1]][2]).unwrap(),
                    ),
                    Vector3::new(
                        cast(mesh.vertices[f.vertices[2]][0]).unwrap(),
                        cast(mesh.vertices[f.vertices[2]][1]).unwrap(),
                        cast(mesh.vertices[f.vertices[2]][2]).unwrap(),
                    ),
                ];
                match Triangle::new(v[0], v[1], v[2]) {
                    Ok(tri) => Some((v, tri)),
                    Err(_) => {
                        degenerate_count += 1;
                        None
                    }
                }
            })
            .collect();

        let (triangles_vertices, triangles_precomputed): (Vec<[Vector3<T>; 3]>, Vec<Triangle<T>>) =
            triangles_data.into_iter().unzip();

        if degenerate_count > 0 {
            eprintln!(
                "{path}: {} triangles (skipped {} degenerate)",
                triangles_vertices.len() + degenerate_count,
                degenerate_count
            );
        } else {
            eprintln!("{path}: {} triangles", triangles_vertices.len());
        }

        let mut rays = Vec::new();
        // Keep the CodSpeed simulation bounded even on the 100k-triangle mesh.
        // Eight samples per axis provide both hits and misses at each mesh size.
        let n = 8;
        let bx_0: T = cast(bounds_x.0).unwrap();
        let bx_1: T = cast(bounds_x.1).unwrap();
        let by_0: T = cast(bounds_y.0).unwrap();
        let by_1: T = cast(bounds_y.1).unwrap();
        let zs: T = cast(z_start).unwrap();
        let neg_one: T = cast(-1.0).unwrap();

        for i in 0..n {
            for j in 0..n {
                let x = (cast::<usize, T>(i).unwrap() / cast::<usize, T>(n).unwrap())
                    * (bx_1 - bx_0)
                    + bx_0;
                let y = (cast::<usize, T>(j).unwrap() / cast::<usize, T>(n).unwrap())
                    * (by_1 - by_0)
                    + by_0;
                rays.push(Ray::new(
                    Vector3::new(x, y, zs),
                    Vector3::new(T::zero(), T::zero(), neg_one),
                ));
            }
        }

        // Smaller, deterministic workloads cover tilted rays and full misses
        // without requiring millions of instrumented triangle tests.
        let mut oblique_rays = Vec::new();
        let mut miss_rays = Vec::new();
        let n = 4;
        let mesh_max_x = mesh
            .vertices
            .iter()
            .map(|v| v[0])
            .fold(f32::NEG_INFINITY, f32::max);
        let mesh_min_x = mesh
            .vertices
            .iter()
            .map(|v| v[0])
            .fold(f32::INFINITY, f32::min);
        let outside_x: T = cast(mesh_max_x + (mesh_max_x - mesh_min_x).max(1.0)).unwrap();
        let slope_x: T = cast(0.15 * (bounds_x.1 - bounds_x.0) / z_start).unwrap();
        let slope_y: T = cast(0.15 * (bounds_y.1 - bounds_y.0) / z_start).unwrap();

        for i in 0..n {
            for j in 0..n {
                let u: T = cast::<usize, T>(i).unwrap() / cast::<usize, T>(n).unwrap();
                let v: T = cast::<usize, T>(j).unwrap() / cast::<usize, T>(n).unwrap();
                let x = u * (bx_1 - bx_0) + bx_0;
                let y = v * (by_1 - by_0) + by_0;
                let dx = if i % 2 == 0 { -slope_x } else { slope_x };
                let dy = if j % 2 == 0 { -slope_y } else { slope_y };
                oblique_rays.push(Ray::new(
                    Vector3::new(x, y, zs),
                    Vector3::new(dx, dy, neg_one),
                ));
                // The origin lies beyond the mesh's actual x extent, so these
                // downward rays traverse every triangle without hitting one.
                miss_rays.push(Ray::new(
                    Vector3::new(outside_x + u * (bx_1 - bx_0), y, zs),
                    Vector3::new(T::zero(), T::zero(), neg_one),
                ));
            }
        }

        Self {
            name: name.to_string(),
            triangles_vertices,
            triangles_precomputed,
            rays,
            oblique_rays,
            miss_rays,
            size,
        }
    }

    pub fn report_hit_rate(&self) {
        let t_min: T = cast(0.0).unwrap();
        let t_max: T = cast(1000.0).unwrap();
        for (workload, rays) in [
            ("axis_64", &self.rays),
            ("oblique_16", &self.oblique_rays),
            ("miss_16", &self.miss_rays),
        ] {
            let mut hits = 0;
            for ray in rays {
                for tri in &self.triangles_precomputed {
                    if tri.intersect(*ray, t_min, t_max).is_some() {
                        hits += 1;
                        break;
                    }
                }
            }
            let hit_rate = hits as f32 / rays.len() as f32;
            eprintln!(
                "{} {workload} hit rate ({}): {:.2}% ({}/{})",
                self.name,
                if std::mem::size_of::<T>() == 4 {
                    "f32"
                } else {
                    "f64"
                },
                hit_rate * 100.0,
                hits,
                rays.len()
            );
        }
    }

    pub fn run_intersect_bench(&self, c: &mut Criterion, precision_label: &str) {
        for (workload, rays) in [
            ("axis_64", &self.rays),
            ("oblique_16", &self.oblique_rays),
            ("miss_16", &self.miss_rays),
        ] {
            // Ray counts changed; retain historical benchmark data under its
            // original IDs and publish these as separate CodSpeed cases.
            let mut group =
                c.benchmark_group(format!("{}/{}/{}", self.name, precision_label, workload));
            if let Some(size) = self.size {
                group.sample_size(size);
            }

            group.bench_function("baldwin_weber", |b| {
                b.iter(|| self.count_baldwin_weber(rays))
            });

            group.bench_function("moller_trumbore", |b| {
                b.iter(|| self.count_moller_trumbore(rays))
            });

            group.finish();
        }
    }

    pub fn run_init_bench(&self, c: &mut Criterion, precision_label: &str) {
        let mut group = c.benchmark_group(format!("{}/{}", self.name, precision_label));
        group.sample_size(10);

        // This ID differs from the old benchmark, which discarded the result
        // and could have optimized the construction away.
        group.bench_function("baldwin_weber_init_observed", |b| {
            b.iter(|| self.initialize_observed())
        });

        group.finish();
    }

    /// Shared kernels for the CodSpeed and local wall-clock harnesses.
    pub fn count_baldwin_weber(&self, rays: &[Ray<T>]) -> usize {
        let mut hits = 0;
        let t_min: T = cast(0.0).unwrap();
        let t_max: T = cast(1000.0).unwrap();
        for ray in black_box(rays) {
            for tri in black_box(&self.triangles_precomputed) {
                if tri.intersect(*ray, t_min, t_max).is_some() {
                    hits += 1;
                    break;
                }
            }
        }
        hits
    }

    pub fn count_moller_trumbore(&self, rays: &[Ray<T>]) -> usize {
        let mut hits = 0;
        let t_min: T = cast(0.0).unwrap();
        let t_max: T = cast(1000.0).unwrap();
        for ray in black_box(rays) {
            for v in black_box(&self.triangles_vertices) {
                if intersect::intersect_moller_trumbore(v[0], v[1], v[2], *ray, t_min, t_max)
                    .is_some()
                {
                    hits += 1;
                    break;
                }
            }
        }
        hits
    }

    pub fn initialize_observed(&self) -> usize {
        let mut count = 0;
        for v in black_box(&self.triangles_vertices) {
            black_box(Triangle::new(v[0], v[1], v[2]).unwrap());
            count += 1;
        }
        count
    }
}

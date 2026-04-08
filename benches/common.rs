use criterion::{Criterion, black_box};
use nalgebra::{RealField, Vector3};
use num_traits::{Float, cast};
use rayx::{Ray, Triangle, intersect};
use std::fs::File;

pub struct ModelBench<T: Float + RealField + Copy> {
    pub name: String,
    pub triangles_vertices: Vec<[Vector3<T>; 3]>,
    pub triangles_precomputed: Vec<Triangle<T>>,
    pub rays: Vec<Ray<T>>,
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
        let n = 50;
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

        Self {
            name: name.to_string(),
            triangles_vertices,
            triangles_precomputed,
            rays,
            size,
        }
    }

    pub fn report_hit_rate(&self) {
        let mut hits = 0;
        let t_min: T = cast(0.0).unwrap();
        let t_max: T = cast(1000.0).unwrap();
        for ray in &self.rays {
            for tri in &self.triangles_precomputed {
                if tri.intersect(*ray, t_min, t_max).is_some() {
                    hits += 1;
                    break;
                }
            }
        }
        let hit_rate = hits as f32 / self.rays.len() as f32;
        eprintln!(
            "{} hit rate ({}): {:.2}% ({}/{})",
            self.name,
            if std::mem::size_of::<T>() == 4 {
                "f32"
            } else {
                "f64"
            },
            hit_rate * 100.0,
            hits,
            self.rays.len()
        );
    }

    pub fn run_intersect_bench(&self, c: &mut Criterion, precision_label: &str) {
        let mut group = c.benchmark_group(format!("{}/{}", self.name, precision_label));
        if let Some(size) = self.size {
            group.sample_size(size);
        }

        let t_min: T = cast(0.0).unwrap();
        let t_max: T = cast(1000.0).unwrap();

        group.bench_function("baldwin_weber", |b| {
            b.iter(|| {
                let mut hits = 0;
                for ray in black_box(&self.rays) {
                    for tri in black_box(&self.triangles_precomputed) {
                        if tri.intersect(*ray, t_min, t_max).is_some() {
                            hits += 1;
                            break;
                        }
                    }
                }
                hits
            })
        });

        group.bench_function("moller_trumbore", |b| {
            b.iter(|| {
                let mut hits = 0;
                for ray in black_box(&self.rays) {
                    for v in black_box(&self.triangles_vertices) {
                        if intersect::intersect_moller_trumbore(
                            v[0], v[1], v[2], *ray, t_min, t_max,
                        )
                        .is_some()
                        {
                            hits += 1;
                            break;
                        }
                    }
                }
                hits
            })
        });

        group.finish();
    }

    pub fn run_init_bench(&self, c: &mut Criterion, precision_label: &str) {
        let mut group = c.benchmark_group(format!("{}/{}", self.name, precision_label));
        group.sample_size(10);

        group.bench_function("baldwin_weber_init", |b| {
            b.iter(|| {
                for v in black_box(&self.triangles_vertices) {
                    let _ = Triangle::new(v[0], v[1], v[2]).unwrap();
                }
            })
        });

        group.finish();
    }
}

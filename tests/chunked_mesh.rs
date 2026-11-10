use nalgebra::{RealField, Vector3};
use num_traits::Float;
use rayx::{Hit, Ray, Triangle, intersect::intersect_baldwin_weber_chunked};

fn check<T: Float + RealField + core::fmt::Debug>(vertices: &[[Vector3<T>; 3]], rays: &[Ray<T>]) {
    let end = T::from(1000.0).unwrap();
    let transformed: Vec<_> = vertices
        .iter()
        .map(|v| Triangle::new(v[0], v[1], v[2]).ok())
        .collect();
    let expected: Vec<Option<(usize, Hit<T>)>> = rays
        .iter()
        .map(|ray| {
            transformed.iter().enumerate().find_map(|(index, v)| {
                v.as_ref()
                    .and_then(|tri| tri.intersect(*ray, T::zero(), end).map(|hit| (index, hit)))
            })
        })
        .collect();

    for size in [512, 2048] {
        let mut scratch = vec![None; size];
        let mut found = vec![None; rays.len()];
        intersect_baldwin_weber_chunked(vertices, rays, &mut scratch, &mut found, T::zero(), end)
            .unwrap();
        assert_eq!(found, expected, "chunk size {size}");
    }
}

#[test]
fn real_meshes_match_sequential_first_hit() {
    for (file_name, height) in [
        ("perfect-suzanne.stl", 5.0_f32),
        ("bun_zipper_res2.stl", 1.0),
        ("thai-statue-decim001.stl", 200.0),
    ] {
        let path = std::path::Path::new("tests/test_data").join(file_name);
        if !path.exists() {
            continue;
        }
        let mesh = stl_io::read_stl(&mut std::fs::File::open(path).unwrap()).unwrap();
        let vertices: Vec<[Vector3<f32>; 3]> = mesh
            .faces
            .iter()
            .map(|f| {
                f.vertices.map(|i| {
                    let v = mesh.vertices[i];
                    Vector3::new(v[0], v[1], v[2])
                })
            })
            .collect();
        let rays = [
            Ray::new([0.0, 0.0, height], [0.0, 0.0, -1.0]),
            Ray::new([0.0, 0.1, height], [0.0, 0.0, -1.0]),
            Ray::new([0.5, -0.1, height], [0.1, 0.0, -1.0]),
            Ray::new([10_000.0, 0.0, height], [0.0, 0.0, -1.0]),
        ];
        check(&vertices, &rays);

        let vertices64: Vec<_> = vertices
            .iter()
            .map(|v| v.map(|x| x.map(f64::from)))
            .collect();
        let rays64: Vec<_> = rays
            .iter()
            .map(|r| Ray::new(r.origin.map(f64::from), r.dir.map(f64::from)))
            .collect();
        check(&vertices64, &rays64);
    }
}

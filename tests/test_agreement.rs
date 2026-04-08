use nalgebra::Vector3;
use rayx::{Ray, Triangle, intersect};
use std::fs::File;
use std::path::Path;

fn approx(a: f32, b: f32, tol: f32) -> bool {
    (a - b).abs() <= tol
}

fn run_stl_agreement_test(file_name: &str, rtol: f32) {
    let test_data_dir = Path::new("tests/test-data");
    let path = test_data_dir.join(file_name);
    if !path.exists() {
        println!(
            "File {} not found, skipping. Please download STL files from the github repo: https://github.com/p-sira/rayx/tests/test-data",
            file_name
        );
        return;
    }

    let mut file = File::open(&path).expect("Failed to open model");
    let mesh = stl_io::read_stl(&mut file).expect("Failed to read model");

    // Load triangles
    let mut triangles = Vec::new();
    for f in &mesh.faces {
        let v = [
            Vector3::new(
                mesh.vertices[f.vertices[0]][0],
                mesh.vertices[f.vertices[0]][1],
                mesh.vertices[f.vertices[0]][2],
            ),
            Vector3::new(
                mesh.vertices[f.vertices[1]][0],
                mesh.vertices[f.vertices[1]][1],
                mesh.vertices[f.vertices[1]][2],
            ),
            Vector3::new(
                mesh.vertices[f.vertices[2]][0],
                mesh.vertices[f.vertices[2]][1],
                mesh.vertices[f.vertices[2]][2],
            ),
        ];
        if let Ok(tri) = Triangle::<f32>::new(v[0], v[1], v[2]) {
            triangles.push((v, tri));
        }
    }

    // Generate some rays aimed at the origin where models are usually centered
    let rays = [
        Ray::new(Vector3::new(0.0, 0.0, 10.0), Vector3::new(0.0, 0.0, -1.0)),
        Ray::new(Vector3::new(0.1, 0.1, 10.0), Vector3::new(0.0, 0.0, -1.0)),
        Ray::new(Vector3::new(-0.1, -0.1, 10.0), Vector3::new(0.0, 0.0, -1.0)),
        Ray::new(Vector3::new(10.0, 0.0, 0.0), Vector3::new(-1.0, 0.0, 0.0)),
        Ray::new(Vector3::new(0.0, 10.0, 0.0), Vector3::new(0.0, -1.0, 0.0)),
    ];

    let t_min = 0.0;
    let t_max = 1000.0;

    for ray in rays {
        for (v, tri) in &triangles {
            let res_bw = tri.intersect(ray, t_min, t_max);
            let res_mt = intersect::intersect_moller_trumbore(v[0], v[1], v[2], ray, t_min, t_max);

            match (res_bw, res_mt) {
                (Some(hb), Some(hm)) => {
                    assert!(
                        approx(hb.t, hm.t, rtol),
                        "t mismatch in {}: BW={:?} vs MT={:?}",
                        file_name,
                        hb,
                        hm
                    );
                    assert!(
                        approx(hb.u, hm.u, rtol),
                        "u mismatch in {}: BW={:?} vs MT={:?}",
                        file_name,
                        hb,
                        hm
                    );
                    assert!(
                        approx(hb.v, hm.v, rtol),
                        "v mismatch in {}: BW={:?} vs MT={:?}",
                        file_name,
                        hb,
                        hm
                    );
                }
                (None, None) => {}
                _ => {
                    panic!(
                        "Disagreement in {}: Ray {:?}, v0={:?}, v1={:?}, v2={:?}\nBW Result: {:?}\nMT Result: {:?}",
                        file_name, ray, v[0], v[1], v[2], res_bw, res_mt
                    );
                }
            }
        }
    }
}

macro_rules! test_stl_agreement {
    ($($(#[$meta:meta])* $name:ident: $rtol:expr),* $(,)?) => {
        $(
            #[test]
            $(#[$meta])*
            fn $name() {
                run_stl_agreement_test(concat!(stringify!($name), ".stl"), $rtol);
            }
        )*
    };
}

test_stl_agreement! {
    bun_zipper_res2: 1e-4,
    bun_zipper_res3: 1e-4,
    bun_zipper_res4: 1e-4,
    perfect_suzanne: 1e-4,
    #[ignore]
    thai_statue: 1e-4,
}

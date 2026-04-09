use crate::common::{approx, indexed_mesh_to_triangle_vertices, try_load_mesh};
use rayx::Triangle;

fn run_reconstruct_vertices_test(file_name: &str, rtol: f32) {
    let mesh = match try_load_mesh(file_name) {
        Some(mesh) => mesh,
        None => return,
    };

    let triangle_vertices = indexed_mesh_to_triangle_vertices(&mesh);
    let triangles = triangle_vertices
        .iter()
        .map(|v| Triangle::new(v[0], v[1], v[2]).unwrap())
        .collect::<Vec<_>>();

    for (v, tri) in triangle_vertices.iter().zip(triangles.iter()) {
        let reconstructed = tri.reconstruct_vertices();
        for i in 0..3 {
            assert!(approx(v[i].x, reconstructed[i].x, rtol));
            assert!(approx(v[i].y, reconstructed[i].y, rtol));
            assert!(approx(v[i].z, reconstructed[i].z, rtol));
        }
    }
}

macro_rules! test_reconstruct_vertices {
    ($($(#[$meta:meta])* $name:ident),* $(,)?) => {
        $(
            #[test]
            $(#[$meta])*
            fn $name() {
                run_reconstruct_vertices_test(concat!(stringify!($name), ".stl"), 1e-5);
            }
        )*
    };
}

test_reconstruct_vertices! {
    bun_zipper_res2,
    bun_zipper_res3,
    bun_zipper_res4,
    perfect_suzanne,
    #[ignore]
    thai_statue,
}

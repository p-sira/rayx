use std::{fs::File, path::Path};

use nalgebra::Vector3;
use stl_io::{IndexedMesh, IndexedTriangle};

pub fn approx(a: f32, b: f32, tol: f32) -> bool {
    (a - b).abs() <= tol
}

pub fn try_load_mesh(file_name: &str) -> Option<IndexedMesh> {
    let test_data_dir = Path::new("tests/test_data");
    let path = test_data_dir.join(file_name);
    if !path.exists() {
        println!(
            "File {} not found, skipping. Please download STL files from the github repo: https://github.com/p-sira/rayx/tests/test_data",
            file_name
        );
        return None;
    }

    let mut file = File::open(&path).expect("Failed to open model");
    let mesh = stl_io::read_stl(&mut file).expect("Failed to read model");
    Some(mesh)
}

pub fn indexed_mesh_to_triangle_vertices(mesh: &IndexedMesh) -> Vec<[Vector3<f32>; 3]> {
    let idx =
        |m: &IndexedMesh, i: usize, f: &IndexedTriangle, j: usize| m.vertices[f.vertices[j]][i];

    mesh.faces
        .iter()
        .map(|f| {
            [
                Vector3::new(idx(mesh, 0, f, 0), idx(mesh, 1, f, 0), idx(mesh, 2, f, 0)),
                Vector3::new(idx(mesh, 0, f, 1), idx(mesh, 1, f, 1), idx(mesh, 2, f, 1)),
                Vector3::new(idx(mesh, 0, f, 2), idx(mesh, 1, f, 2), idx(mesh, 2, f, 2)),
            ]
        })
        .collect()
}

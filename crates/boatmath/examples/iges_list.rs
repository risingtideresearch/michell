// Scratch: list an IGES file's surfaces.
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let f = hullgeom::iges::parse(&std::fs::read_to_string(path).unwrap()).unwrap();
    println!("units {} entities {:?}", f.units_scale, f.entity_counts);
    let mut rows = Vec::new();
    for (i, s) in f.surfaces.iter().enumerate() {
        let b = |k: usize| {
            s.ctrl
                .iter()
                .map(|p| p[k])
                .fold((f64::INFINITY, f64::NEG_INFINITY), |m, v| {
                    (m.0.min(v), m.1.max(v))
                })
        };
        let (x, y, z) = (b(0), b(1), b(2));
        let w = s
            .weights
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |m, &v| {
                (m.0.min(v), m.1.max(v))
            });
        rows.push(format!("{i:3} deg {}x{} n {:3}x{:3} rational {:5} w {:.3}..{:.3} trim {} x {:7.3}..{:7.3} y {:7.3}..{:7.3} z {:7.3}..{:7.3}",
            s.degree_u, s.degree_v, s.n_ctrl_u, s.n_ctrl_v, !s.is_polynomial(), w.0, w.1, s.trim_uv.is_some(), x.0, x.1, y.0, y.1, z.0, z.1));
    }
    for r in rows {
        println!("{r}");
    }
}

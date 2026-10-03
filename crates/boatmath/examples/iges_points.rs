// Scratch: sample every surface of an IGES file (weights honoured), as "id x y z" lines.
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let f = hullgeom::iges::parse(&std::fs::read_to_string(path).unwrap()).unwrap();
    for (i, s) in f.surfaces.iter().enumerate() {
        let ((u0, u1), (v0, v1)) = match s.trim_uv {
            Some([a, b, c, d]) => ((a, b), (c, d)),
            None => (s.u_domain(), s.v_domain()),
        };
        let (nu, nv) = (120, 120);
        for a in 0..=nu {
            for b in 0..=nv {
                let (u, v) = (
                    u0 + (u1 - u0) * a as f64 / nu as f64,
                    v0 + (v1 - v0) * b as f64 / nv as f64,
                );
                // Rational: weighted sum by hand via the homogeneous trick on the ctrl net.
                let p = s.point(u, v);
                println!("{i} {:.4} {:.4} {:.4}", p[0], p[1], p[2]);
            }
        }
    }
}

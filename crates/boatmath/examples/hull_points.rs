// Scratch: sample a hull record's patches as "x y z" lines.
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let r: serde_json::Value = serde_json::from_str(
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .next()
            .unwrap(),
    )
    .unwrap();
    for p in r["geometry"]["hulls"][0]["patches"].as_array().unwrap() {
        let s = boatmath::native::patch_from_value(p).unwrap();
        let ((u0, u1), (v0, v1)) = match s.trim_uv {
            Some([a, b, c, d]) => ((a, b), (c, d)),
            None => (s.u_domain(), s.v_domain()),
        };
        let n = 200;
        for a in 0..=n {
            for b in 0..=n {
                let q = s.point(
                    u0 + (u1 - u0) * a as f64 / n as f64,
                    v0 + (v1 - v0) * b as f64 / n as f64,
                );
                println!("{:.4} {:.4} {:.4}", q[0], q[1], q[2]);
            }
        }
    }
}

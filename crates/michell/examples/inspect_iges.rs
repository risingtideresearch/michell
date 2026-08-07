//! Inventory the B-spline surfaces in an IGES file (debugging aid).
//!
//! Usage: cargo run --example inspect_iges -- file.igs

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: inspect_iges <file.igs>");
    let text = std::fs::read_to_string(&path).expect("read file");
    let file = michell::iges::parse(&text).expect("parse");
    println!("units scale: {}", file.units_scale);
    println!("entities: {:?}", file.entity_counts);
    println!("surfaces: {}", file.surfaces.len());
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    let mut rational = 0usize;
    let mut degrees = std::collections::BTreeMap::<(usize, usize), usize>::new();
    for s in &file.surfaces {
        if !s.is_polynomial() {
            rational += 1;
        }
        *degrees.entry((s.degree_u, s.degree_v)).or_default() += 1;
        for p in &s.ctrl {
            for c in 0..3 {
                lo[c] = lo[c].min(p[c]);
                hi[c] = hi[c].max(p[c]);
            }
        }
    }
    println!("rational: {rational}");
    println!("degrees (u,v) -> count: {degrees:?}");
    println!(
        "control-net bbox [m]: x {:.3}..{:.3}  y {:.3}..{:.3}  z {:.3}..{:.3}",
        lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]
    );
    // Knot-span structure of the first few patches.
    for (i, s) in file.surfaces.iter().take(3).enumerate() {
        println!(
            "patch {i}: ctrl {}x{}, knots-u {:?}, knots-v {:?}",
            s.n_ctrl_u, s.n_ctrl_v, s.knots_u, s.knots_v
        );
    }
}

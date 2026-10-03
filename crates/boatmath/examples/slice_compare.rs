// Scratch: two hull records sliced square to x, as "which x y z" lines; the second shifted by dx.
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let load = |p: &str| -> serde_json::Value {
        serde_json::from_str(std::fs::read_to_string(p).unwrap().lines().next().unwrap()).unwrap()
    };
    let (t, f) = (load(&a[1]), load(&a[2]));
    let dx: f64 = a[3].parse().unwrap();
    let tt = boatmath::camber::fit::Target::from_geometry(&t["geometry"]).unwrap();
    let ff = boatmath::camber::fit::Target::from_geometry(&f["geometry"]).unwrap();
    for x in a[4].split(',').map(|v| v.parse::<f64>().unwrap()) {
        for (w, g, s) in [("hull", &tt, 0.0), ("fit", &ff, dx)] {
            for seg in g.slice([x - s, 0.0], [1.0, 0.0]) {
                for p in seg {
                    println!("{w} {x} {:.5} {:.5}", -p[0], p[1]);
                }
            }
        }
    }
}

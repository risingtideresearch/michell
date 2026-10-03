// Scratch: score a camber document against a hull record's geometry.
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let hull: serde_json::Value = serde_json::from_str(
        std::fs::read_to_string(&a[1])
            .unwrap()
            .lines()
            .next()
            .unwrap(),
    )
    .unwrap();
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&a[2]).unwrap()).unwrap();
    let d = boatmath::camber::Document::parse(&v).unwrap();
    let z0: f64 = a[3].parse().unwrap();
    let x0: f64 = a[4].parse().unwrap();
    println!(
        "{:?}",
        boatmath::camber::fit::score(&hull["geometry"], &d, z0, x0).unwrap()
    );
}

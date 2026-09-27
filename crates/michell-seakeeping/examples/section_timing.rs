//! Cost of one section heave solve at a few panel counts.
use michell_seakeeping::section2d::Section;

fn main() {
    for n in [12, 20, 32] {
        let s = Section::semicircle(0.15, n);
        let t = std::time::Instant::now();
        let reps = 50;
        for i in 0..reps {
            let _ = s.heave(3.0 + 0.1 * i as f64, 9.81, 1000.0).unwrap();
        }
        println!("{n} panels: {:.2} ms per solve", t.elapsed().as_secs_f64() * 1e3 / reps as f64);
    }
}

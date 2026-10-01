//! Gauss–Legendre nodes and weights on [-1, 1].

/// Nodes and weights of the n-point Gauss–Legendre rule on [-1, 1].
///
/// Newton iteration on the Legendre polynomial recurrence; accurate to
/// machine precision for the modest n used in this crate.
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    assert!(n >= 1);
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    let m = n.div_ceil(2);
    for i in 0..m {
        // Tricomi initial guess for the i-th root (descending order).
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            // Evaluate P_n(x) and P_{n-1}(x) by recurrence.
            let mut p0 = 1.0;
            let mut p1 = x;
            for k in 2..=n {
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            // p1 = P_n, p0 = P_{n-1}
            dp = (n as f64) * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() <= 1e-16 * x.abs().max(1.0) {
                // One extra polish iteration, then stop.
                let mut q0 = 1.0;
                let mut q1 = x;
                for k in 2..=n {
                    let kf = k as f64;
                    let q2 = ((2.0 * kf - 1.0) * x * q1 - (kf - 1.0) * q0) / kf;
                    q0 = q1;
                    q1 = q2;
                }
                dp = (n as f64) * (x * q1 - q0) / (x * x - 1.0);
                x -= q1 / dp;
                break;
            }
        }
        let w = 2.0 / ((1.0 - x * x) * dp * dp);
        nodes[i] = -x;
        nodes[n - 1 - i] = x;
        weights[i] = w;
        weights[n - 1 - i] = w;
    }
    if n % 2 == 1 {
        // The middle node of an odd rule is exactly zero.
        nodes[n / 2] = 0.0;
    }
    (nodes, weights)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrates_polynomials_exactly() {
        // n-point GL is exact for degree 2n-1.
        for n in [1usize, 2, 3, 5, 8, 16, 24] {
            let (x, w) = gauss_legendre(n);
            for deg in 0..(2 * n) {
                let approx: f64 = x
                    .iter()
                    .zip(&w)
                    .map(|(&xi, &wi)| wi * xi.powi(deg as i32))
                    .sum();
                let exact = if deg % 2 == 0 {
                    2.0 / (deg as f64 + 1.0)
                } else {
                    0.0
                };
                assert!(
                    (approx - exact).abs() < 1e-13,
                    "n={n} deg={deg}: {approx} vs {exact}"
                );
            }
        }
    }
}

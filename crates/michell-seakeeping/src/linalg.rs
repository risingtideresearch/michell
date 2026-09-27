//! Dense complex linear solves (panel systems are small: tens of unknowns).

use michell_geometry::C64;

/// Solve `A x = b` for square `A` (row-major, `n × n`) by LU with partial
/// pivoting; `None` if `A` is singular to working precision.
pub fn solve(mut a: Vec<C64>, mut b: Vec<C64>) -> Option<Vec<C64>> {
    let n = b.len();
    debug_assert_eq!(a.len(), n * n);
    for k in 0..n {
        let (piv, mag) = (k..n)
            .map(|i| (i, a[i * n + k].abs()))
            .fold((k, -1.0), |best, c| if c.1 > best.1 { c } else { best });
        if !(mag > 0.0) {
            return None;
        }
        if piv != k {
            for j in 0..n {
                a.swap(k * n + j, piv * n + j);
            }
            b.swap(k, piv);
        }
        let inv = a[k * n + k].recip();
        for i in k + 1..n {
            let l = a[i * n + k] * inv;
            if l == C64::ZERO {
                continue;
            }
            for j in k + 1..n {
                let v = a[k * n + j];
                a[i * n + j] = a[i * n + j] - l * v;
            }
            b[i] = b[i] - l * b[k];
        }
    }
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in i + 1..n {
            s = s - a[i * n + j] * b[j];
        }
        b[i] = s / a[i * n + i];
    }
    Some(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_a_small_complex_system() {
        let a = vec![
            C64::new(0.0, 1.0), C64::new(2.0, 0.0), C64::new(1.0, -1.0),
            C64::new(3.0, 0.0), C64::new(0.0, 0.0), C64::new(1.0, 0.0),
            C64::new(1.0, 1.0), C64::new(-1.0, 0.5), C64::new(0.0, 2.0),
        ];
        let x = vec![C64::new(1.0, 2.0), C64::new(-0.5, 0.0), C64::new(0.3, -1.0)];
        let b: Vec<C64> = (0..3)
            .map(|i| (0..3).fold(C64::ZERO, |s, j| s + a[i * 3 + j] * x[j]))
            .collect();
        let got = solve(a, b).unwrap();
        for (g, w) in got.iter().zip(&x) {
            assert!((*g - *w).abs() < 1e-12);
        }
    }
}

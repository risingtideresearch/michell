//! Closed-form moment integrals against oscillatory and exponential kernels.
//!
//! These are the primitives that make the Michell inner integrals exact for
//! piecewise-polynomial (B-spline) hulls: on every knot span the integrand is
//! a polynomial times `exp(i k x)` in x and a polynomial times `exp(-κ z)` in
//! z, and both families of moments below have closed forms.

use std::ops::{Add, Mul, Sub};

/// Minimal complex number (kept local to avoid any dependency).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

impl C64 {
    pub const ZERO: C64 = C64 { re: 0.0, im: 0.0 };

    #[inline]
    pub fn new(re: f64, im: f64) -> C64 {
        C64 { re, im }
    }

    /// e^{iθ}
    #[inline]
    pub fn cis(theta: f64) -> C64 {
        let (s, c) = theta.sin_cos();
        C64 { re: c, im: s }
    }

    #[inline]
    pub fn scale(self, s: f64) -> C64 {
        C64 {
            re: self.re * s,
            im: self.im * s,
        }
    }

    #[inline]
    pub fn abs_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    #[inline]
    pub fn abs(self) -> f64 {
        self.abs_sq().sqrt()
    }

    /// Complex exponential `e^{z}`.
    #[inline]
    pub fn exp(self) -> C64 {
        let r = self.re.exp();
        let (s, c) = self.im.sin_cos();
        C64 {
            re: r * c,
            im: r * s,
        }
    }

    /// Reciprocal `1/z = conj(z)/|z|²`.
    #[inline]
    pub fn recip(self) -> C64 {
        let d = self.abs_sq();
        C64 {
            re: self.re / d,
            im: -self.im / d,
        }
    }
}

impl Add for C64 {
    type Output = C64;
    #[inline]
    fn add(self, o: C64) -> C64 {
        C64::new(self.re + o.re, self.im + o.im)
    }
}

impl Sub for C64 {
    type Output = C64;
    #[inline]
    fn sub(self, o: C64) -> C64 {
        C64::new(self.re - o.re, self.im - o.im)
    }
}

impl Mul for C64 {
    type Output = C64;
    #[inline]
    fn mul(self, o: C64) -> C64 {
        C64::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}

/// Switch between power series (small argument, avoids cancellation) and the
/// closed-form recurrence (large argument, exact).
const SERIES_THRESHOLD: f64 = 4.0;

/// Oscillatory moments `M_a = ∫_0^h t^a e^{i k t} dt` for `a = 0..=a_max`.
///
/// For |kh| below [`SERIES_THRESHOLD`] uses the entire series
/// `M_a = h^{a+1} Σ_m (ikh)^m / (m! (a+m+1))`; above it, the stable upward
/// recurrence `M_a = (h^a e^{ikh} - a M_{a-1}) / (ik)`.
pub fn osc_moments(k: f64, h: f64, a_max: usize, out: &mut Vec<C64>) {
    debug_assert!(h > 0.0);
    out.clear();
    let kh = k * h;
    if kh.abs() <= SERIES_THRESHOLD {
        let ikh = C64::new(0.0, kh);
        let mut h_pow = h;
        for a in 0..=a_max {
            let mut term = C64::new(1.0, 0.0); // (ikh)^m / m!
            let mut sum = C64::ZERO;
            for m in 0..80 {
                let contrib = term.scale(1.0 / (a as f64 + m as f64 + 1.0));
                sum = sum + contrib;
                if contrib.abs() <= 1e-18 * sum.abs() {
                    break;
                }
                term = term * ikh.scale(1.0 / (m as f64 + 1.0));
            }
            out.push(sum.scale(h_pow));
            h_pow *= h;
        }
    } else {
        let e = C64::cis(kh);
        let inv_ik = C64::new(0.0, -1.0 / k); // 1/(ik)
        let mut prev = (e - C64::new(1.0, 0.0)) * inv_ik;
        out.push(prev);
        let mut h_pow = h;
        for a in 1..=a_max {
            let cur = (e.scale(h_pow) - prev.scale(a as f64)) * inv_ik;
            out.push(cur);
            prev = cur;
            h_pow *= h;
        }
    }
}

/// Exponential moments `N_b = ∫_0^h t^b e^{-κ t} dt` for `b = 0..=b_max`, κ ≥ 0.
///
/// Small κh uses the entire series `N_b = h^{b+1} Σ_m (-κh)^m / (m! (b+m+1))`;
/// large κh the upward recurrence `N_b = (b N_{b-1} - h^b e^{-κh}) / κ`.
pub fn exp_moments(kappa: f64, h: f64, b_max: usize, out: &mut Vec<f64>) {
    debug_assert!(h > 0.0);
    debug_assert!(kappa >= 0.0);
    out.clear();
    let x = kappa * h;
    if x <= SERIES_THRESHOLD {
        let mut h_pow = h;
        for b in 0..=b_max {
            let mut term = 1.0f64; // (-x)^m / m!
            let mut sum = 0.0f64;
            for m in 0..80 {
                let contrib = term / (b as f64 + m as f64 + 1.0);
                sum += contrib;
                if contrib.abs() <= 1e-18 * sum.abs() {
                    break;
                }
                term *= -x / (m as f64 + 1.0);
            }
            out.push(sum * h_pow);
            h_pow *= h;
        }
    } else {
        let e = (-x).exp();
        let mut prev = (1.0 - e) / kappa;
        out.push(prev);
        let mut h_pow = h;
        for b in 1..=b_max {
            let cur = (b as f64 * prev - h_pow * e) / kappa;
            out.push(cur);
            prev = cur;
            h_pow *= h;
        }
    }
}

/// Exponential moments `N_b = ∫_0^h t^b e^{-κ t} dt` for a **complex** decay
/// `κ` with `Re(κ) >= 0` — the same series/recurrence as [`exp_moments`], in
/// `C64`. Used only by the heeled-hull kernel, where the tilt makes the
/// vertical decay complex (`κ = νλ²cosφ + i νλ√(λ²−1) sinφ`); the upright path
/// keeps the real routine untouched.
pub fn exp_moments_complex(kappa: C64, h: f64, b_max: usize, out: &mut Vec<C64>) {
    debug_assert!(h > 0.0);
    debug_assert!(kappa.re >= 0.0);
    out.clear();
    let x = kappa.scale(h); // κh
    if x.abs() <= SERIES_THRESHOLD {
        let mut h_pow = h;
        for b in 0..=b_max {
            let mut term = C64::new(1.0, 0.0); // (-x)^m / m!
            let mut sum = C64::ZERO;
            for m in 0..80 {
                let contrib = term.scale(1.0 / (b as f64 + m as f64 + 1.0));
                sum = sum + contrib;
                if contrib.abs() <= 1e-18 * sum.abs() {
                    break;
                }
                term = term * x.scale(-1.0 / (m as f64 + 1.0));
            }
            out.push(sum.scale(h_pow));
            h_pow *= h;
        }
    } else {
        let e = x.scale(-1.0).exp(); // e^{-x}
        let inv_k = kappa.recip();
        let mut prev = (C64::new(1.0, 0.0) - e) * inv_k;
        out.push(prev);
        let mut h_pow = h;
        for b in 1..=b_max {
            let cur = (prev.scale(b as f64) - e.scale(h_pow)) * inv_k;
            out.push(cur);
            prev = cur;
            h_pow *= h;
        }
    }
}

/// Ordered pair moments
/// `T_{b,b'} = ∫_0^h dZ ∫_0^Z dζ  Z^b ζ^{b'} e^{-κ(Z-ζ)}` for `b, b' = 0..=b_max`,
/// κ ≥ 0, laid out row-major: `out[b*(b_max+1) + b'] = T_{b,b'}`.
///
/// This is the single-span primitive for the non-separable kernel
/// `sgn(z-ζ) e^{-κ|z-ζ|}` — the unbounded-fluid (Munk-type) part of the
/// thin-ship pitching moment. The caller wants the antisymmetric combination
/// `T_{b,b'} - T_{b',b}`; the symmetric one `T_{b,b'} + T_{b',b}` is the
/// unordered `∫_0^h∫_0^h Z^b ζ^{b'} e^{-κ|Z-ζ|}`, which the tests lean on.
///
/// Small κh expands `e^{-κ(Z-ζ)} = Σ_m (-κ)^m (Z-ζ)^m / m!` and uses the Beta
/// integral `∫_0^Z ζ^{b'} (Z-ζ)^m dζ = b'! m! Z^{b'+m+1} / (b'+m+1)!`, giving
/// the entire series `T = h^{b+b'+2} Σ_m (-κh)^m b'! / ((b'+m+1)! (b+b'+m+2))`.
/// Large κh integrates the inner ζ integral by parts,
/// `∫_0^Z ζ^{b'} e^{κζ} dζ = e^{κZ} Σ_{j≤b'} (-1)^j (b'!/(b'-j)!) Z^{b'-j}/κ^{j+1}
///   + (-1)^{b'+1} b'!/κ^{b'+1}`,
/// after which the outer integral is elementary except for the last term,
/// which is `N_b` from [`exp_moments`]:
/// `T = Σ_{j≤b'} (-1)^j (b'!/(b'-j)!) h^{b+b'-j+1} / (κ^{j+1} (b+b'-j+1))
///   + (-1)^{b'+1} b'! N_b / κ^{b'+1}`.
/// The closed form cancels catastrophically as κh → 0 (every term is
/// O(h^{b+b'+1}/κ) while `T` is O(h^{b+b'+2})), hence the usual switch at
/// [`SERIES_THRESHOLD`]. Underflow of `e^{-κh}` at extreme κ is absorbed by
/// [`exp_moments`] and the result stays finite (`T → h^{b+b'+1}/(κ(b+b'+1))`).
#[allow(dead_code)] // wired up by the pitching-moment kernel; not yet called
pub fn exp_pair_moments_ordered(kappa: f64, h: f64, b_max: usize, out: &mut Vec<f64>) {
    debug_assert!(h > 0.0);
    debug_assert!(kappa >= 0.0);
    if kappa * h <= SERIES_THRESHOLD {
        exp_pair_moments_ordered_series(kappa, h, b_max, out);
    } else {
        exp_pair_moments_ordered_closed(kappa, h, b_max, out);
    }
}

/// Series branch of [`exp_pair_moments_ordered`]; valid for any κh, accurate
/// (and cheap) only for small κh.
fn exp_pair_moments_ordered_series(kappa: f64, h: f64, b_max: usize, out: &mut Vec<f64>) {
    out.clear();
    let x = kappa * h;
    for b in 0..=b_max {
        for bp in 0..=b_max {
            let mut term = 1.0 / (bp as f64 + 1.0); // (-x)^m b'! / (b'+m+1)!
            let mut sum = 0.0f64;
            for m in 0..80 {
                let contrib = term / ((b + bp + m) as f64 + 2.0);
                sum += contrib;
                if contrib.abs() <= 1e-18 * sum.abs() {
                    break;
                }
                term *= -x / ((bp + m) as f64 + 2.0);
            }
            out.push(sum * h.powi((b + bp) as i32 + 2));
        }
    }
}

/// Closed-form branch of [`exp_pair_moments_ordered`]; exact for any κ > 0
/// but cancels as κh → 0.
fn exp_pair_moments_ordered_closed(kappa: f64, h: f64, b_max: usize, out: &mut Vec<f64>) {
    out.clear();
    let mut n_b = Vec::new();
    exp_moments(kappa, h, b_max, &mut n_b);
    let inv_k = 1.0 / kappa;
    for b in 0..=b_max {
        for bp in 0..=b_max {
            let mut fall = 1.0f64; // b'!/(b'-j)!
            let mut h_pow = h.powi((b + bp) as i32 + 1); // h^{b+b'-j+1}
            let mut k_pow = inv_k; // κ^{-(j+1)}
            let mut sign = 1.0f64; // (-1)^j
            let mut sum = 0.0f64;
            for j in 0..=bp {
                if j > 0 {
                    fall *= (bp + 1 - j) as f64;
                    h_pow /= h;
                    k_pow *= inv_k;
                    sign = -sign;
                }
                sum += sign * fall * k_pow * h_pow / ((b + bp - j) as f64 + 1.0);
            }
            // Leaving the loop: fall = b'!, k_pow = κ^{-(b'+1)}, sign = (-1)^{b'}.
            sum -= sign * fall * k_pow * n_b[b];
            out.push(sum);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Composite Simpson reference for ∫_0^h f(t) dt.
    fn simpson<F: Fn(f64) -> f64>(f: F, h: f64, n: usize) -> f64 {
        let n = n + n % 2;
        let dt = h / n as f64;
        let mut s = f(0.0) + f(h);
        for i in 1..n {
            let w = if i % 2 == 1 { 4.0 } else { 2.0 };
            s += w * f(i as f64 * dt);
        }
        s * dt / 3.0
    }

    #[test]
    fn osc_moments_match_quadrature() {
        let mut out = Vec::new();
        // Spans both the series and recurrence branches, incl. near-threshold.
        for &(k, h) in &[
            (1e-9, 1.0),
            (0.5, 1.0),
            (3.9, 1.0),
            (4.1, 1.0),
            (25.0, 0.7),
            (-13.0, 2.3),
            (400.0, 0.05),
        ] {
            osc_moments(k, h, 5, &mut out);
            #[allow(clippy::needless_range_loop)]
            for a in 0..=5usize {
                let re = simpson(|t| t.powi(a as i32) * (k * t).cos(), h, 20000);
                let im = simpson(|t| t.powi(a as i32) * (k * t).sin(), h, 20000);
                let scale = h.powi(a as i32 + 1) / (a as f64 + 1.0);
                assert!(
                    (out[a].re - re).abs() < 1e-9 * scale && (out[a].im - im).abs() < 1e-9 * scale,
                    "k={k} h={h} a={a}: got ({}, {}), want ({re}, {im})",
                    out[a].re,
                    out[a].im
                );
            }
        }
    }

    #[test]
    fn exp_moments_match_quadrature() {
        let mut out = Vec::new();
        for &(kappa, h) in &[
            (0.0, 1.0),
            (1e-8, 2.0),
            (0.5, 1.0),
            (3.9, 1.0),
            (4.1, 1.0),
            (30.0, 0.7),
            (900.0, 0.5),
        ] {
            exp_moments(kappa, h, 5, &mut out);
            #[allow(clippy::needless_range_loop)]
            for b in 0..=5usize {
                let want = simpson(|t| t.powi(b as i32) * (-kappa * t).exp(), h, 20000);
                let scale = (h.powi(b as i32 + 1) / (b as f64 + 1.0)).max(want.abs());
                assert!(
                    (out[b] - want).abs() < 1e-9 * scale,
                    "kappa={kappa} h={h} b={b}: got {}, want {want}",
                    out[b]
                );
            }
        }
    }

    #[test]
    fn exp_moments_complex_match_quadrature() {
        let mut out = Vec::new();
        let mut real = Vec::new();
        for &(kr, ki, h) in &[
            (0.5, 0.0, 1.0),   // real ⇒ must match exp_moments
            (0.3, 0.7, 1.0),   // series branch
            (2.0, -3.0, 1.2),  // near/over threshold
            (20.0, 15.0, 0.7), // recurrence branch
            (5.0, 40.0, 0.3),  // strongly oscillatory decay
        ] {
            let kappa = C64::new(kr, ki);
            exp_moments_complex(kappa, h, 5, &mut out);
            #[allow(clippy::needless_range_loop)]
            for b in 0..=5usize {
                let re = simpson(
                    |t| t.powi(b as i32) * (-kr * t).exp() * (ki * t).cos(),
                    h,
                    40000,
                );
                let im = simpson(
                    |t| t.powi(b as i32) * (-kr * t).exp() * -(ki * t).sin(),
                    h,
                    40000,
                );
                let scale = h.powi(b as i32 + 1) / (b as f64 + 1.0);
                assert!(
                    (out[b].re - re).abs() < 1e-9 * scale && (out[b].im - im).abs() < 1e-9 * scale,
                    "κ=({kr},{ki}) h={h} b={b}: got ({}, {}), want ({re}, {im})",
                    out[b].re,
                    out[b].im
                );
            }
            // A real decay must reproduce the real routine bit-for-bit closely.
            if ki == 0.0 {
                exp_moments(kr, h, 5, &mut real);
                for b in 0..=5usize {
                    assert!((out[b].re - real[b]).abs() < 1e-14 && out[b].im.abs() < 1e-300);
                }
            }
        }
    }

    #[test]
    fn extreme_decay_underflows_gracefully() {
        let mut out = Vec::new();
        exp_moments(1e6, 1.0, 3, &mut out);
        // N_b -> b!/κ^{b+1}
        assert!((out[0] - 1e-6).abs() < 1e-18);
        assert!((out[1] - 1e-12).abs() < 1e-24);
        assert!(out.iter().all(|v| v.is_finite()));
    }

    /// Nested composite-Simpson reference for the pair moments, all `(b, b')`
    /// at once: `∫_0^h dZ ∫_0^Z dζ Z^b ζ^{b'} e^{-κ(Z-ζ)}` when `ordered`, else
    /// the full square `∫_0^h dZ ∫_0^h dζ Z^b ζ^{b'} e^{-κ|Z-ζ|}` with the inner
    /// integral split at the kink ζ = Z so each piece is smooth.
    fn pair_quadrature(kappa: f64, h: f64, b_max: usize, ordered: bool, n: usize) -> Vec<f64> {
        let nb = b_max + 1;
        let n = n + n % 2;
        let weight = |i: usize| -> f64 {
            if i == 0 || i == n {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            }
        };
        // ∫_lo^hi ζ^{b'} e^{-κ|z-ζ|} dζ for every b', accumulated into `acc`.
        let inner = |z: f64, lo: f64, hi: f64, acc: &mut [f64]| {
            let dz = (hi - lo) / n as f64;
            for i in 0..=n {
                let zeta = lo + i as f64 * dz;
                let e = weight(i) * (-kappa * (z - zeta).abs()).exp() * dz / 3.0;
                for bp in 0..nb {
                    acc[bp] += e * zeta.powi(bp as i32);
                }
            }
        };
        let mut out = vec![0.0; nb * nb];
        let mut row = vec![0.0; nb];
        let dz = h / n as f64;
        for i in 0..=n {
            let z = i as f64 * dz;
            row.iter_mut().for_each(|v| *v = 0.0);
            inner(z, 0.0, z, &mut row);
            if !ordered {
                inner(z, z, h, &mut row);
            }
            for b in 0..nb {
                let f = weight(i) * z.powi(b as i32) * dz / 3.0;
                for bp in 0..nb {
                    out[b * nb + bp] += f * row[bp];
                }
            }
        }
        out
    }

    #[test]
    fn exp_pair_moments_ordered_match_quadrature() {
        let mut out = Vec::new();
        // Spans both branches, incl. either side of the threshold.
        for &(kappa, h) in &[
            (0.0, 1.0),
            (1e-8, 2.0),
            (0.5, 1.0),
            (3.9, 1.0),
            (4.1, 1.0),
            (20.0, 0.7),
            (200.0, 0.05),
        ] {
            exp_pair_moments_ordered(kappa, h, 3, &mut out);
            let want = pair_quadrature(kappa, h, 3, true, 2000);
            for b in 0..=3usize {
                for bp in 0..=3usize {
                    let i = b * 4 + bp;
                    assert!(
                        (out[i] - want[i]).abs() < 1e-9 * want[i].abs(),
                        "kappa={kappa} h={h} b={b} b'={bp}: got {}, want {}",
                        out[i],
                        want[i]
                    );
                }
            }
        }
    }

    #[test]
    fn exp_pair_moments_symmetric_sum_is_unordered_integral() {
        let mut out = Vec::new();
        for &(kappa, h) in &[(0.0, 1.0), (1.5, 1.0), (4.0, 1.0), (12.0, 0.8)] {
            exp_pair_moments_ordered(kappa, h, 3, &mut out);
            let want = pair_quadrature(kappa, h, 3, false, 2000);
            for b in 0..=3usize {
                for bp in 0..=3usize {
                    let i = b * 4 + bp;
                    let got = out[i] + out[bp * 4 + b];
                    assert!(
                        (got - want[i]).abs() < 1e-9 * want[i].abs(),
                        "kappa={kappa} h={h} b={b} b'={bp}: got {got}, want {}",
                        want[i]
                    );
                    // The unordered integral is symmetric in (b, b') (to
                    // quadrature accuracy — the reference's inner grids differ).
                    assert!((want[i] - want[bp * 4 + b]).abs() < 1e-9 * want[i].abs());
                }
            }
        }
    }

    #[test]
    fn exp_pair_moments_zero_decay_closed_form() {
        let mut out = Vec::new();
        for &h in &[0.3, 1.0, 2.5] {
            exp_pair_moments_ordered(0.0, h, 4, &mut out);
            for b in 0..=4usize {
                for bp in 0..=4usize {
                    // ∫_0^h Z^b · Z^{b'+1}/(b'+1) dZ
                    let want =
                        h.powi((b + bp) as i32 + 2) / ((bp as f64 + 1.0) * ((b + bp) as f64 + 2.0));
                    let got = out[b * 5 + bp];
                    assert!(
                        (got - want).abs() < 1e-14 * want,
                        "h={h} b={b} b'={bp}: got {got}, want {want}"
                    );
                    // With ζ < Z on the domain, T_{b,b'} > T_{b',b} whenever b > b'.
                    if b > bp {
                        assert!(got > out[bp * 5 + b]);
                    }
                }
            }
        }
    }

    #[test]
    fn exp_pair_moments_branches_agree_at_threshold() {
        let mut series = Vec::new();
        let mut closed = Vec::new();
        for &(kappa, h) in &[(4.0, 1.0), (2.0, 2.0), (3.9, 1.0), (4.1, 1.0), (8.0, 0.5)] {
            exp_pair_moments_ordered_series(kappa, h, 3, &mut series);
            exp_pair_moments_ordered_closed(kappa, h, 3, &mut closed);
            for i in 0..16 {
                assert!(
                    (series[i] - closed[i]).abs() < 1e-13 * series[i].abs(),
                    "kappa={kappa} h={h} i={i}: series {}, closed {}",
                    series[i],
                    closed[i]
                );
            }
        }
    }

    #[test]
    fn exp_pair_moments_extreme_decay_underflows_gracefully() {
        let mut out = Vec::new();
        exp_pair_moments_ordered(1e6, 1.0, 3, &mut out);
        // T_{0,0} = (κh - 1 + e^{-κh})/κ² -> 1/κ - 1/κ².
        assert!((out[0] - (1e-6 - 1e-12)).abs() < 1e-18);
        // T_{1,0} = h²/(2κ) - N_1/κ -> 1/(2κ) - 1/κ³.
        assert!((out[4] - (0.5e-6 - 1e-18)).abs() < 1e-18);
        assert!(out.iter().all(|v| v.is_finite() && *v > 0.0));
    }
}

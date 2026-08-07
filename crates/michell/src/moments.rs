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
}

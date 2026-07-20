//! Clamped, polynomial (non-rational) tensor-product B-spline surfaces.
//!
//! This is the canonical hull representation: `y = S(x, z)` with `y` the local
//! half-beam. Rational surfaces (NURBS with non-unit weights) are deliberately
//! not supported — polynomial spans are what make the Michell inner integrals
//! evaluable in closed form (see [`crate::moments`]).

use crate::error::{Error, Result};

#[derive(Debug, Clone)]
pub struct BSplineSurface {
    degree_x: usize,
    degree_z: usize,
    knots_x: Vec<f64>,
    knots_z: Vec<f64>,
    /// Row-major: `control[i * n_ctrl_z + j]` for x-index `i`, z-index `j`.
    control: Vec<f64>,
    n_ctrl_x: usize,
    n_ctrl_z: usize,
}

impl BSplineSurface {
    /// Build and validate a clamped tensor-product B-spline surface.
    ///
    /// Requirements per direction (degree `p`, `n` control points):
    /// - `p >= 1`, `n >= p + 1`, `knots.len() == n + p + 1`;
    /// - knots finite and non-decreasing, positive-length domain;
    /// - clamped: first and last knot each appear exactly `p + 1` times;
    /// - interior knot multiplicity at most `p` (multiplicity `p` gives a
    ///   C⁰ crease — how chines are represented).
    ///
    /// `control` is row-major with the z index fastest:
    /// `control[i * n_ctrl_z + j]`.
    pub fn new(
        degree_x: usize,
        degree_z: usize,
        knots_x: Vec<f64>,
        knots_z: Vec<f64>,
        control: Vec<f64>,
    ) -> Result<Self> {
        let n_ctrl_x = validate_knots(&knots_x, degree_x, "x")?;
        let n_ctrl_z = validate_knots(&knots_z, degree_z, "z")?;
        if control.len() != n_ctrl_x * n_ctrl_z {
            return Err(Error::InvalidSpline(format!(
                "control net has {} entries, expected {} x {} = {}",
                control.len(),
                n_ctrl_x,
                n_ctrl_z,
                n_ctrl_x * n_ctrl_z
            )));
        }
        if control.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidSpline(
                "control net contains a non-finite value".into(),
            ));
        }
        Ok(BSplineSurface {
            degree_x,
            degree_z,
            knots_x,
            knots_z,
            control,
            n_ctrl_x,
            n_ctrl_z,
        })
    }

    pub fn degree_x(&self) -> usize {
        self.degree_x
    }

    pub fn degree_z(&self) -> usize {
        self.degree_z
    }

    pub fn knots_x(&self) -> &[f64] {
        &self.knots_x
    }

    pub fn knots_z(&self) -> &[f64] {
        &self.knots_z
    }

    pub fn control(&self) -> &[f64] {
        &self.control
    }

    pub fn n_ctrl_x(&self) -> usize {
        self.n_ctrl_x
    }

    pub fn n_ctrl_z(&self) -> usize {
        self.n_ctrl_z
    }

    /// Parametric domain in x: `[knots_x[p], knots_x[n]]`.
    pub fn x_domain(&self) -> (f64, f64) {
        (self.knots_x[self.degree_x], self.knots_x[self.n_ctrl_x])
    }

    /// Parametric domain in z: `[knots_z[q], knots_z[n]]`.
    pub fn z_domain(&self) -> (f64, f64) {
        (self.knots_z[self.degree_z], self.knots_z[self.n_ctrl_z])
    }

    /// Surface value. Arguments outside the domain are clamped to it.
    pub fn eval(&self, x: f64, z: f64) -> f64 {
        self.eval_deriv(x, z, 0, 0)
    }

    /// Mixed partial derivative `∂^{dx+dz} S / ∂x^{dx} ∂z^{dz}`.
    /// Arguments outside the domain are clamped to it.
    pub fn eval_deriv(&self, x: f64, z: f64, dx: usize, dz: usize) -> f64 {
        let (x0, x1) = self.x_domain();
        let (z0, z1) = self.z_domain();
        let x = x.clamp(x0, x1);
        let z = z.clamp(z0, z1);
        let sx = find_span(&self.knots_x, self.degree_x, self.n_ctrl_x, x);
        let sz = find_span(&self.knots_z, self.degree_z, self.n_ctrl_z, z);
        let ndx = ders_basis(&self.knots_x, self.degree_x, sx, x, dx);
        let ndz = ders_basis(&self.knots_z, self.degree_z, sz, z, dz);
        if dx > self.degree_x || dz > self.degree_z {
            return 0.0;
        }
        let mut sum = 0.0;
        for (i, &bx) in ndx[dx].iter().enumerate() {
            let ci = sx - self.degree_x + i;
            let row = &self.control[ci * self.n_ctrl_z..];
            let mut inner = 0.0;
            for (j, &bz) in ndz[dz].iter().enumerate() {
                inner += bz * row[sz - self.degree_z + j];
            }
            sum += bx * inner;
        }
        sum
    }

    /// Indices `s` of non-empty knot spans `[knots_x[s], knots_x[s+1])`.
    pub(crate) fn x_span_indices(&self) -> Vec<usize> {
        span_indices(&self.knots_x, self.degree_x, self.n_ctrl_x)
    }

    pub(crate) fn z_span_indices(&self) -> Vec<usize> {
        span_indices(&self.knots_z, self.degree_z, self.n_ctrl_z)
    }

    /// All mixed partials `D[a][b] = ∂^{a+b} S / ∂x^a ∂z^b` for
    /// `a = 0..=degree_x`, `b = 0..=degree_z`, evaluated at the lower-left
    /// corner of the given (non-empty) span pair. Together with Taylor's
    /// theorem this yields the exact local polynomial on the span rectangle.
    pub(crate) fn corner_partials(&self, span_x: usize, span_z: usize) -> Vec<Vec<f64>> {
        let p = self.degree_x;
        let q = self.degree_z;
        let x0 = self.knots_x[span_x];
        let z0 = self.knots_z[span_z];
        let ndx = ders_basis(&self.knots_x, p, span_x, x0, p);
        let ndz = ders_basis(&self.knots_z, q, span_z, z0, q);
        let mut d = vec![vec![0.0; q + 1]; p + 1];
        for (a, row_a) in d.iter_mut().enumerate() {
            for (b, dab) in row_a.iter_mut().enumerate() {
                let mut sum = 0.0;
                for (i, &bx) in ndx[a].iter().enumerate() {
                    let ci = span_x - p + i;
                    let row = &self.control[ci * self.n_ctrl_z..];
                    let mut inner = 0.0;
                    for (j, &bz) in ndz[b].iter().enumerate() {
                        inner += bz * row[span_z - q + j];
                    }
                    sum += bx * inner;
                }
                *dab = sum;
            }
        }
        d
    }
}

/// Validate one knot vector; returns the control-point count `n`.
fn validate_knots(knots: &[f64], degree: usize, dir: &str) -> Result<usize> {
    if degree < 1 {
        return Err(Error::InvalidSpline(format!(
            "{dir}: degree must be at least 1"
        )));
    }
    if knots.len() < 2 * (degree + 1) {
        return Err(Error::InvalidSpline(format!(
            "{dir}: need at least {} knots for degree {degree}, got {}",
            2 * (degree + 1),
            knots.len()
        )));
    }
    let n = knots.len() - degree - 1;
    if knots.iter().any(|k| !k.is_finite()) {
        return Err(Error::InvalidSpline(format!(
            "{dir}: knot vector contains a non-finite value"
        )));
    }
    if knots.windows(2).any(|w| w[1] < w[0]) {
        return Err(Error::InvalidSpline(format!(
            "{dir}: knot vector is not non-decreasing"
        )));
    }
    if knots[degree] >= knots[n] {
        return Err(Error::InvalidSpline(format!(
            "{dir}: knot domain has zero length"
        )));
    }
    // Clamped: first and last values appear exactly degree + 1 times.
    let first = knots[0];
    let last = *knots.last().unwrap();
    let first_mult = knots.iter().take_while(|&&k| k == first).count();
    let last_mult = knots.iter().rev().take_while(|&&k| k == last).count();
    if first_mult != degree + 1 || last_mult != degree + 1 {
        return Err(Error::InvalidSpline(format!(
            "{dir}: knot vector must be clamped (first and last knot each \
             repeated exactly degree+1 = {} times; found {first_mult} and {last_mult})",
            degree + 1
        )));
    }
    // Interior multiplicity at most `degree`.
    let interior = &knots[degree + 1..n];
    let mut run = 1usize;
    for w in interior.windows(2) {
        if w[1] == w[0] {
            run += 1;
            if run > degree {
                return Err(Error::InvalidSpline(format!(
                    "{dir}: interior knot {} has multiplicity greater than degree {degree}",
                    w[0]
                )));
            }
        } else {
            run = 1;
        }
    }
    Ok(n)
}

/// Non-zero B-spline basis values and first derivatives at `u`: returns the
/// index of the first non-zero basis function and two rows of `degree + 1`
/// entries — `rows[0]` the basis values, `rows[1]` their first derivatives.
pub(crate) fn basis_rows1(
    knots: &[f64],
    degree: usize,
    n_ctrl: usize,
    u: f64,
) -> (usize, Vec<Vec<f64>>) {
    let span = find_span(knots, degree, n_ctrl, u);
    let ders = ders_basis(knots, degree, span, u, 1);
    (span - degree, ders)
}

/// Index `s` such that `knots[s] <= u < knots[s+1]` within the domain,
/// with the right end mapped into the last non-empty span.
pub(crate) fn find_span(knots: &[f64], degree: usize, n_ctrl: usize, u: f64) -> usize {
    if u >= knots[n_ctrl] {
        // Last non-empty span.
        let mut s = n_ctrl - 1;
        while s > degree && knots[s + 1] <= knots[s] {
            s -= 1;
        }
        return s;
    }
    if u <= knots[degree] {
        return degree;
    }
    let (mut lo, mut hi) = (degree, n_ctrl);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if u < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo
}

fn span_indices(knots: &[f64], degree: usize, n_ctrl: usize) -> Vec<usize> {
    (degree..n_ctrl)
        .filter(|&s| knots[s + 1] > knots[s])
        .collect()
}

/// Non-zero basis functions and derivatives at `u` (Piegl & Tiller A2.3).
///
/// Returns `ders[k][j]` = k-th derivative of basis function `N_{span-p+j, p}`
/// at `u`, for `k = 0..=min(n, p)` and `j = 0..=p`.
pub(crate) fn ders_basis(knots: &[f64], p: usize, span: usize, u: f64, n: usize) -> Vec<Vec<f64>> {
    let n = n.min(p);
    let mut ndu = vec![vec![0.0f64; p + 1]; p + 1];
    ndu[0][0] = 1.0;
    let mut left = vec![0.0f64; p + 1];
    let mut right = vec![0.0f64; p + 1];
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            // Lower triangle: knot differences.
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            // Upper triangle: basis values.
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let mut ders = vec![vec![0.0f64; p + 1]; n + 1];
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }
    if n == 0 {
        return ders;
    }
    let mut a = [vec![0.0f64; p + 1], vec![0.0f64; p + 1]];
    for r in 0..=p {
        let mut s1 = 0usize;
        let mut s2 = 1usize;
        a[0].iter_mut().for_each(|v| *v = 0.0);
        a[1].iter_mut().for_each(|v| *v = 0.0);
        a[0][0] = 1.0;
        for k in 1..=n {
            let mut d = 0.0;
            let rk = r as isize - k as isize;
            let pk = p - k;
            if r >= k {
                a[s2][0] = a[s1][0] / ndu[pk + 1][(rk) as usize];
                d = a[s2][0] * ndu[rk as usize][pk];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if r as isize - 1 <= pk as isize {
                k - 1
            } else {
                p - r
            };
            for j in j1..=j2 {
                a[s2][j] =
                    (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][(rk + j as isize) as usize];
                d += a[s2][j] * ndu[(rk + j as isize) as usize][pk];
            }
            if r <= pk {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut factor = p as f64;
    for (k, row) in ders.iter_mut().enumerate().skip(1) {
        for v in row.iter_mut() {
            *v *= factor;
        }
        factor *= (p - k) as f64;
    }
    ders
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Wigley-style biquadratic: exactly representable at degree (2, 2).
    fn wigley_surface(l: f64, b: f64, t: f64) -> BSplineSurface {
        let a = l / 2.0;
        let knots_x = vec![-a, -a, -a, a, a, a];
        let knots_z = vec![0.0, 0.0, 0.0, t, t, t];
        let gx = [0.0, 2.0, 0.0];
        let hz = [1.0, 1.0, 0.0];
        let mut control = Vec::with_capacity(9);
        for gi in gx {
            for hj in hz {
                control.push(b / 2.0 * gi * hj);
            }
        }
        BSplineSurface::new(2, 2, knots_x, knots_z, control).unwrap()
    }

    fn wigley_exact(l: f64, b: f64, t: f64, x: f64, z: f64, dx: usize, dz: usize) -> f64 {
        // f = (b/2) g(x) h(z), g = 1-(2x/l)^2, h = 1-(z/t)^2
        let g = match dx {
            0 => 1.0 - (2.0 * x / l).powi(2),
            1 => -8.0 * x / (l * l),
            2 => -8.0 / (l * l),
            _ => 0.0,
        };
        let h = match dz {
            0 => 1.0 - (z / t).powi(2),
            1 => -2.0 * z / (t * t),
            2 => -2.0 / (t * t),
            _ => 0.0,
        };
        b / 2.0 * g * h
    }

    #[test]
    fn reproduces_biquadratic_and_derivatives() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let s = wigley_surface(l, b, t);
        for i in 0..=10 {
            for j in 0..=10 {
                let x = -l / 2.0 + l * i as f64 / 10.0;
                let z = t * j as f64 / 10.0;
                for dx in 0..=3 {
                    for dz in 0..=3 {
                        let got = s.eval_deriv(x, z, dx, dz);
                        let want = wigley_exact(l, b, t, x, z, dx, dz);
                        assert!(
                            (got - want).abs() < 1e-10 * (1.0 + want.abs()),
                            "x={x} z={z} dx={dx} dz={dz}: got {got}, want {want}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn corner_partials_match_eval_deriv() {
        // Multi-span cubic x quadratic surface with a smooth control net.
        let knots_x = vec![0.0, 0.0, 0.0, 0.0, 2.5, 5.0, 7.5, 10.0, 10.0, 10.0, 10.0];
        let knots_z = vec![0.0, 0.0, 0.0, 0.6, 1.2, 1.2, 1.2];
        let (nx, nz) = (7usize, 4usize);
        let mut control = vec![0.0; nx * nz];
        for i in 0..nx {
            for j in 0..nz {
                control[i * nz + j] =
                    (1.0 + (i as f64 * 0.9).sin().abs()) * (1.5 - 0.3 * j as f64);
            }
        }
        let s = BSplineSurface::new(3, 2, knots_x, knots_z, control).unwrap();
        for &sx in &s.x_span_indices() {
            for &sz in &s.z_span_indices() {
                let d = s.corner_partials(sx, sz);
                let x0 = s.knots_x()[sx];
                let z0 = s.knots_z()[sz];
                for (a, row) in d.iter().enumerate() {
                    for (b, &dab) in row.iter().enumerate() {
                        let want = s.eval_deriv(x0, z0, a, b);
                        assert!(
                            (dab - want).abs() < 1e-9 * (1.0 + want.abs()),
                            "span ({sx},{sz}) a={a} b={b}: {dab} vs {want}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn taylor_reconstruction_matches_surface() {
        // The corner partials must reconstruct the surface anywhere in the span.
        let knots_x = vec![0.0, 0.0, 0.0, 0.0, 2.5, 5.0, 7.5, 10.0, 10.0, 10.0, 10.0];
        let knots_z = vec![0.0, 0.0, 0.0, 0.6, 1.2, 1.2, 1.2];
        let (nx, nz) = (7usize, 4usize);
        let mut control = vec![0.0; nx * nz];
        for i in 0..nx {
            for j in 0..nz {
                control[i * nz + j] = (0.3 + i as f64).sqrt() * (2.0 - 0.4 * j as f64);
            }
        }
        let s = BSplineSurface::new(3, 2, knots_x, knots_z, control).unwrap();
        let mut fact = [1.0f64; 8];
        for i in 1..8 {
            fact[i] = fact[i - 1] * i as f64;
        }
        for &sx in &s.x_span_indices() {
            for &sz in &s.z_span_indices() {
                let d = s.corner_partials(sx, sz);
                let x0 = s.knots_x()[sx];
                let x1 = s.knots_x()[sx + 1];
                let z0 = s.knots_z()[sz];
                let z1 = s.knots_z()[sz + 1];
                for fx in [0.13, 0.5, 0.97] {
                    for fz in [0.21, 0.5, 0.88] {
                        let x = x0 + fx * (x1 - x0);
                        let z = z0 + fz * (z1 - z0);
                        let mut taylor = 0.0;
                        for (a, row) in d.iter().enumerate() {
                            for (b, &dab) in row.iter().enumerate() {
                                taylor += dab / (fact[a] * fact[b])
                                    * (x - x0).powi(a as i32)
                                    * (z - z0).powi(b as i32);
                            }
                        }
                        let want = s.eval(x, z);
                        assert!(
                            (taylor - want).abs() < 1e-9 * (1.0 + want.abs()),
                            "span ({sx},{sz}) x={x} z={z}: {taylor} vs {want}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn rejects_bad_input() {
        // Unclamped knot vector.
        assert!(BSplineSurface::new(
            2,
            2,
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0; 9],
        )
        .is_err());
        // Decreasing knots.
        assert!(BSplineSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, -1.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0; 12],
        )
        .is_err());
        // Wrong control-net size.
        assert!(BSplineSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0; 8],
        )
        .is_err());
        // Interior multiplicity above degree (C^{-1} split).
        assert!(BSplineSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0; 18],
        )
        .is_err());
        // Non-finite control point.
        assert!(BSplineSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 0.0, f64::NAN, 0.0, 0.0, 0.0, 0.0],
        )
        .is_err());
    }

    #[test]
    fn chine_via_repeated_interior_knot() {
        // Degree-2 with an interior knot of multiplicity 2: C0 crease allowed.
        let s = BSplineSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
            vec![1.0; 3 * 5],
        );
        assert!(s.is_ok());
    }
}

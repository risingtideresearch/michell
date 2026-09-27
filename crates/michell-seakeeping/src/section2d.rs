//! Two-dimensional **section hydrodynamics** in heave by the Frank
//! close-fit source method (Frank 1967).
//!
//! A section is its starboard half-curve from the waterline to the keel,
//! cut into straight panels and mirrored to port; a constant source strength
//! on each panel (equal on the mirror, heave being symmetric) is found so
//! that the normal velocity matches the body's at every panel midpoint. The
//! source is the pulsating free-surface source of [`crate::green`]: its
//! logarithms are integrated over each panel in closed form, its bounded
//! wave part by Gauss quadrature.
//!
//! Conventions: section plane `(y, z)`, `z` up, fluid below `z = 0`; time
//! dependence `e^{−iωt}`; normals point **into the fluid**. The heave
//! potential is per unit **velocity** (`∂ψ/∂n = n_z`), so the pressure of a
//! heave velocity `V` is `iωρVψ`, and the force `−∫pn_z ds = (ω²a + iωb)η`
//! gives the added mass and damping per unit length
//! `a = −ρ Re∫ψ n_z ds`, `b = −ωρ Im∫ψ n_z ds` (both sides).
//!
//! Like every source method it fails at the **irregular frequencies** of
//! the interior (the sloshing modes of the fluid the section would enclose,
//! the first near `ν ≈ (π/B) coth(πT/B)` for a box of beam `B`, draft `T`).
//! They are narrow, and betray themselves: the damping stops matching the
//! energy the far field carries. [`Section::heave`] checks that balance on
//! every solve and bridges a failure by interpolating across it.

use crate::green::{far_factor, wave_part};
use crate::linalg::solve;
use michell_geometry::quadrature::gauss_legendre;
use michell_geometry::C64;
use std::f64::consts::PI;

/// A section's starboard half, as panels from the waterline to the keel.
#[derive(Debug, Clone)]
pub struct Section {
    /// Panel end points `(y, z)`, `z` up, from the waterline to the keel.
    nodes: Vec<[f64; 2]>,
}

/// Panel geometry derived once per section.
#[derive(Debug, Clone, Copy)]
struct Panel {
    a: [f64; 2],
    b: [f64; 2],
    mid: [f64; 2],
    len: f64,
    /// Unit normal into the fluid.
    n: [f64; 2],
}

/// A section's heave radiation solution at one frequency.
#[derive(Debug, Clone)]
pub struct HeaveSolution {
    /// Frequency `ω` [rad/s].
    pub omega: f64,
    /// Added mass per unit length `a₃₃` [kg/m].
    pub added_mass: f64,
    /// Damping per unit length `b₃₃` [kg/(m·s)].
    pub damping: f64,
    /// Far-field amplitude: `ψ → C e^{νz} e^{iν|y|}` per unit heave
    /// velocity as `|y| → ∞`.
    pub far: C64,
    /// Whether this solution was interpolated across an irregular
    /// frequency (see [`Section::heave`]).
    pub interpolated: bool,
    density: f64,
    panels: Vec<Panel>,
    /// `ψ` at the panel midpoints (starboard; port is equal).
    psi: Vec<C64>,
}

/// Largest relative disagreement between the near-field damping and the
/// radiated energy `ρω|C|²` a solution may have before it is taken to be
/// polluted by an irregular frequency. Converged solutions away from those
/// agree to about a percent.
pub const ENERGY_TOLERANCE: f64 = 0.05;

impl Section {
    /// A section from its curve as the sectional hull samples it —
    /// `(half-beam, depth below the waterline)` from the waterline round to
    /// the keel — resampled into `panels` panels, clustered (cosine) toward
    /// both ends. `None` for a section with no wetted girth.
    pub fn from_curve(curve: &[(f64, f64)], panels: usize) -> Option<Section> {
        if curve.len() < 2 || panels < 2 {
            return None;
        }
        let pts: Vec<[f64; 2]> = curve
            .iter()
            .map(|&(y, d)| [y.max(0.0), -d.max(0.0)])
            .collect();
        let mut s = vec![0.0];
        for w in pts.windows(2) {
            let d = ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2)).sqrt();
            s.push(s.last().unwrap() + d);
        }
        let total = *s.last().unwrap();
        if !(total > 0.0) {
            return None;
        }
        let at = |t: f64| -> [f64; 2] {
            let target = t * total;
            let i = match s.iter().position(|&v| v >= target) {
                Some(0) => 1,
                Some(i) => i,
                None => s.len() - 1,
            };
            let (s0, s1) = (s[i - 1], s[i]);
            let f = if s1 > s0 {
                (target - s0) / (s1 - s0)
            } else {
                0.0
            };
            [
                pts[i - 1][0] + f * (pts[i][0] - pts[i - 1][0]),
                pts[i - 1][1] + f * (pts[i][1] - pts[i - 1][1]),
            ]
        };
        let nodes = (0..=panels)
            .map(|i| at(0.5 * (1.0 - (PI * i as f64 / panels as f64).cos())))
            .collect();
        Some(Section { nodes })
    }

    /// A section from its panel end points `(y, z)`, `z` up, from the
    /// waterline to the keel.
    pub fn from_nodes(nodes: Vec<[f64; 2]>) -> Section {
        Section { nodes }
    }

    /// A semicircle of radius `r` (a heaving half-immersed cylinder).
    pub fn semicircle(r: f64, panels: usize) -> Section {
        let nodes = (0..=panels)
            .map(|i| {
                let th = 0.5 * PI * i as f64 / panels as f64;
                [r * th.cos(), -r * th.sin()]
            })
            .collect();
        Section { nodes }
    }

    /// A rectangle of half-beam `b` and draft `t`, `panels` split between
    /// side and bottom in proportion to their lengths.
    pub fn rectangle(b: f64, t: f64, panels: usize) -> Section {
        let ns = ((panels as f64 * t / (t + b)).round() as usize).clamp(1, panels - 1);
        let nb = panels - ns;
        let mut nodes = Vec::new();
        for i in 0..=ns {
            let u = 0.5 * (1.0 - (PI * i as f64 / ns as f64).cos());
            nodes.push([b, -t * u]);
        }
        for i in 1..=nb {
            let u = 0.5 * (1.0 - (PI * i as f64 / nb as f64).cos());
            nodes.push([b * (1.0 - u), -t]);
        }
        Section { nodes }
    }

    /// The panel end points `(y, z)`, `z` up, from the waterline to the keel.
    pub fn nodes(&self) -> &[[f64; 2]] {
        &self.nodes
    }

    /// Half-beam at the waterline.
    pub fn half_beam(&self) -> f64 {
        self.nodes[0][0]
    }

    /// Draft: the depth of the deepest node.
    pub fn draft(&self) -> f64 {
        self.nodes.iter().map(|p| -p[1]).fold(0.0, f64::max)
    }

    fn panels(&self) -> Vec<Panel> {
        self.nodes
            .windows(2)
            .filter_map(|w| {
                let (a, b) = (w[0], w[1]);
                let (dy, dz) = (b[0] - a[0], b[1] - a[1]);
                let len = (dy * dy + dz * dz).sqrt();
                (len > 0.0).then(|| Panel {
                    a,
                    b,
                    mid: [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])],
                    len,
                    n: [-dz / len, dy / len],
                })
            })
            .collect()
    }

    /// The heave radiation problem at frequency `omega`.
    ///
    /// Near an irregular frequency the source method's solution is wrong
    /// but still finite; the tell is that its damping no longer matches the
    /// energy its far field carries. When the two disagree by more than
    /// [`ENERGY_TOLERANCE`], the problem is re-solved a little below and
    /// above `omega` and the two interpolated — the irregular frequencies
    /// of a section are narrow — and [`HeaveSolution::interpolated`] says so.
    pub fn heave(&self, omega: f64, gravity: f64, density: f64) -> Option<HeaveSolution> {
        let sol = self.heave_at(omega, gravity, density)?;
        if sol.energy_error() <= ENERGY_TOLERANCE {
            return Some(sol);
        }
        for eps in [0.02, 0.04, 0.08] {
            let lo = self.heave_at(omega * (1.0 - eps), gravity, density)?;
            let hi = self.heave_at(omega * (1.0 + eps), gravity, density)?;
            if lo.energy_error() <= ENERGY_TOLERANCE && hi.energy_error() <= ENERGY_TOLERANCE {
                let mix = |a: C64, b: C64| (a + b).scale(0.5);
                return Some(HeaveSolution {
                    omega,
                    added_mass: 0.5 * (lo.added_mass + hi.added_mass),
                    // Damping scales with frequency through b = ρω|C|²:
                    // interpolate |C|² and ψ, rebuild b at ω.
                    damping: 0.5 * (lo.damping / lo.omega + hi.damping / hi.omega) * omega,
                    far: mix(lo.far, hi.far),
                    psi: lo
                        .psi
                        .iter()
                        .zip(&hi.psi)
                        .map(|(&a, &b)| mix(a, b))
                        .collect(),
                    panels: lo.panels,
                    interpolated: true,
                    density,
                });
            }
        }
        Some(sol)
    }

    fn heave_at(&self, omega: f64, gravity: f64, density: f64) -> Option<HeaveSolution> {
        let panels = self.panels();
        let nu = omega * omega / gravity;
        let (s, d) = influence(&panels, nu);
        let n = panels.len();
        let rhs: Vec<C64> = panels.iter().map(|p| C64::new(p.n[1], 0.0)).collect();
        let sigma = solve(d, rhs)?;
        let psi = apply(&s, &sigma, n);
        let force: C64 = panels.iter().zip(&psi).fold(C64::ZERO, |acc, (p, &f)| {
            acc + f.scale(2.0 * p.n[1] * p.len)
        });
        let far = far_amplitude(&panels, &sigma, nu);
        Some(HeaveSolution {
            omega,
            added_mass: -density * force.re,
            damping: -omega * density * force.im,
            far,
            panels,
            psi,
            interpolated: false,
            density,
        })
    }

    /// The head-seas diffraction force per unit length and unit wave
    /// amplitude, `−∫ p_D n_z ds`, by solving the diffraction problem itself
    /// at the wave's own frequency (`ν = k`) — the reference the Haskind
    /// form [`HeaveSolution::diffraction`] is checked against.
    pub fn diffraction_direct(&self, omega: f64, gravity: f64, density: f64) -> Option<C64> {
        let panels = self.panels();
        let nu = omega * omega / gravity;
        let (s, d) = influence(&panels, nu);
        let n = panels.len();
        // φ_I = −i(g/ω) e^{kz}; ∂φ_D/∂n = −∂φ_I/∂n = i(gk/ω) e^{kz} n_z.
        let rhs: Vec<C64> = panels
            .iter()
            .map(|p| C64::new(0.0, gravity * nu / omega * (nu * p.mid[1]).exp() * p.n[1]))
            .collect();
        let sigma = solve(d, rhs)?;
        let phi = apply(&s, &sigma, n);
        let integral = panels.iter().zip(&phi).fold(C64::ZERO, |acc, (p, &f)| {
            acc + f.scale(2.0 * p.n[1] * p.len)
        });
        Some(C64::new(0.0, -omega * density) * integral)
    }
}

impl HeaveSolution {
    /// Relative disagreement between the damping and the radiated energy,
    /// `|b/(ρω|C|²) − 1|`.
    pub fn energy_error(&self) -> f64 {
        let far = self.density * self.omega * self.far.abs_sq();
        if far > 0.0 {
            (self.damping / far - 1.0).abs()
        } else if self.damping.abs() > 0.0 {
            f64::INFINITY
        } else {
            0.0
        }
    }

    /// Radiated wave amplitude per unit heave amplitude,
    /// `|ζ/η| = ω²|C|/g`.
    pub fn wave_ratio(&self, gravity: f64) -> f64 {
        self.omega * self.omega * self.far.abs() / gravity
    }

    /// The diffraction heave force per unit length and unit wave amplitude
    /// of a wave of wavenumber `k` at heading `heading` (`π` head seas), by
    /// the Haskind relation on this radiation solution:
    ///
    /// ```text
    /// f_D = ρ g k ∫ ψ e^{kz} [n_z cos(k y sin β) − sin β n_y sin(k y sin β)] ds
    /// ```
    ///
    /// exact (within strip theory) when this solution's frequency is the
    /// wave's own; strip theory at speed evaluates `ψ` at the encounter
    /// frequency instead (see [`crate::strip`]). Phase relative to a crest
    /// over the section's centreplane.
    pub fn diffraction(&self, k: f64, heading: f64, gravity: f64, density: f64) -> C64 {
        let sb = heading.sin();
        let sum = self
            .panels
            .iter()
            .zip(&self.psi)
            .fold(C64::ZERO, |acc, (p, &f)| {
                let (y, z) = (p.mid[0], p.mid[1]);
                let arg = k * y * sb;
                let w = (k * z).exp() * (p.n[1] * arg.cos() - sb * p.n[0] * arg.sin());
                acc + f.scale(2.0 * w * p.len)
            });
        sum.scale(density * gravity * k)
    }
}

/// `ψ = S σ`.
fn apply(s: &[C64], sigma: &[C64], n: usize) -> Vec<C64> {
    (0..n)
        .map(|i| (0..n).fold(C64::ZERO, |acc, j| acc + s[i * n + j] * sigma[j]))
        .collect()
}

/// The far-field amplitude `C = −i Σ_j σ_j ∫_j 2 e^{νζ} cos νη ds`.
fn far_amplitude(panels: &[Panel], sigma: &[C64], nu: f64) -> C64 {
    let (gx, gw) = gauss_legendre(8);
    let mut sum = C64::ZERO;
    for (p, &sg) in panels.iter().zip(sigma) {
        let mut f = C64::ZERO;
        for (&t, &w) in gx.iter().zip(&gw) {
            let u = 0.5 * (1.0 + t);
            let q = [
                p.a[0] + u * (p.b[0] - p.a[0]),
                p.a[1] + u * (p.b[1] - p.a[1]),
            ];
            f = f
                + (far_factor(nu, q[0], q[1], 1.0) + far_factor(nu, q[0], q[1], -1.0))
                    .scale(0.5 * w * p.len);
        }
        sum = sum + sg * f;
    }
    C64::new(0.0, -1.0) * sum
}

/// Distance from `p` to the segment `a → b`.
fn seg_distance(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let (dy, dz) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dy * dy + dz * dz;
    let t = if l2 > 0.0 {
        (((p[0] - a[0]) * dy + (p[1] - a[1]) * dz) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    ((p[0] - a[0] - t * dy).powi(2) + (p[1] - a[1] - t * dz).powi(2)).sqrt()
}

/// `∫ ln|p − q| ds_q` over the segment `a → b`, and its gradient in `p`.
/// `on_panel` puts `p` on the segment's own line exactly (a panel's own
/// midpoint): its normal offset is then taken as zero rather than as
/// whatever roundoff leaves, whose sign would otherwise flip the principal
/// value of the normal derivative by `±π`.
fn log_panel(p: [f64; 2], a: [f64; 2], b: [f64; 2], on_panel: bool) -> (f64, [f64; 2]) {
    let (dy, dz) = (b[0] - a[0], b[1] - a[1]);
    let len = (dy * dy + dz * dz).sqrt();
    let t = [dy / len, dz / len];
    let nl = [-t[1], t[0]];
    let rel = [p[0] - a[0], p[1] - a[1]];
    let x0 = rel[0] * t[0] + rel[1] * t[1];
    let y0 = if on_panel {
        0.0
    } else {
        rel[0] * nl[0] + rel[1] * nl[1]
    };
    let (u1, u2) = (-x0, len - x0);
    let f = |u: f64| -> f64 {
        let q = u * u + y0 * y0;
        let log_term = if u == 0.0 { 0.0 } else { 0.5 * u * q.ln() };
        let atan_term = if y0 == 0.0 { 0.0 } else { y0 * (u / y0).atan() };
        log_term - u + atan_term
    };
    let value = f(u2) - f(u1);
    let gx = -0.5 * (u2 * u2 + y0 * y0).ln() + 0.5 * (u1 * u1 + y0 * y0).ln();
    let gy = if y0 == 0.0 {
        0.0
    } else {
        (u2 / y0).atan() - (u1 / y0).atan()
    };
    (value, [gx * t[0] + gy * nl[0], gx * t[1] + gy * nl[1]])
}

/// The influence matrices of the panels and their port mirrors (the source
/// of strength 1 on panel `j` and its mirror, normalised `1/2π`): `S` the
/// potential and `D` the normal velocity at each midpoint, `D` carrying the
/// self term `½` of the fluid-side limit.
fn influence(panels: &[Panel], nu: f64) -> (Vec<C64>, Vec<C64>) {
    let n = panels.len();
    let (gx, gw) = gauss_legendre(8);
    let mut s = vec![C64::ZERO; n * n];
    let mut d = vec![C64::ZERO; n * n];
    let inv2pi = 1.0 / (2.0 * PI);
    let mirror = |q: [f64; 2]| [-q[0], q[1]];
    let image = |q: [f64; 2]| [q[0], -q[1]];
    for (i, pi) in panels.iter().enumerate() {
        let p = pi.mid;
        for (j, pj) in panels.iter().enumerate() {
            // Logarithms: the panel and its mirror, and both their images.
            // The images' integrals `∫ ln r₁` also carry the logarithmic
            // singularity of the wave part's vertical gradient (below).
            let mut val = 0.0;
            let mut grad = [0.0; 2];
            let mut image_logs = 0.0;
            for (m, (a, b)) in [
                (pj.a, pj.b),
                (mirror(pj.b), mirror(pj.a)),
                (image(pj.a), image(pj.b)),
                (image(mirror(pj.b)), image(mirror(pj.a))),
            ]
            .into_iter()
            .enumerate()
            {
                let (v, g) = log_panel(p, a, b, m == 0 && i == j);
                val += v;
                grad[0] += g[0];
                grad[1] += g[1];
                if m >= 2 {
                    image_logs += v;
                }
            }
            // The wave part. Its value and horizontal gradient are bounded;
            // its vertical gradient `−2ν Re P` goes as `2ν ln(ν r₁)` where the
            // panel meets the waterline, so that logarithm is subtracted
            // here and added back in closed form from `∫ ln r₁`. Panels are
            // subdivided in proportion to their length over their distance.
            let near = seg_distance(p, pj.a, pj.b)
                .min(seg_distance(p, mirror(pj.a), mirror(pj.b)))
                .max(1e-12 * pj.len);
            let subs = ((2.0 * pj.len / near).ceil() as usize).clamp(1, 64);
            let mut wv = C64::ZERO;
            let mut wg = [C64::ZERO; 2];
            for sub in 0..subs {
                for (&t, &w) in gx.iter().zip(&gw) {
                    let u = (sub as f64 + 0.5 * (1.0 + t)) / subs as f64;
                    let q = [
                        pj.a[0] + u * (pj.b[0] - pj.a[0]),
                        pj.a[1] + u * (pj.b[1] - pj.a[1]),
                    ];
                    let ds = 0.5 * w * pj.len / subs as f64;
                    for qy in [q[0], -q[0]] {
                        let (dy, zs) = (p[0] - qy, p[1] + q[1]);
                        let (v, gy, gz) = wave_part(nu, dy, zs);
                        let r1 = (dy * dy + zs * zs).sqrt();
                        let gz_reg = if r1 > 0.0 {
                            gz - C64::new(2.0 * nu * (nu * r1).ln(), 0.0)
                        } else {
                            gz
                        };
                        wv = wv + v.scale(ds);
                        wg[0] = wg[0] + gy.scale(ds);
                        wg[1] = wg[1] + gz_reg.scale(ds);
                    }
                }
            }
            // ∫ 2ν ln(ν r₁) over the panel and its mirror.
            wg[1] = wg[1] + C64::new(2.0 * nu * (image_logs + 2.0 * pj.len * nu.ln()), 0.0);
            s[i * n + j] = (C64::new(val, 0.0) + wv).scale(inv2pi);
            let dn = C64::new(grad[0] * pi.n[0] + grad[1] * pi.n[1], 0.0)
                + wg[0].scale(pi.n[0])
                + wg[1].scale(pi.n[1]);
            d[i * n + j] = dn.scale(inv2pi)
                + if i == j {
                    C64::new(0.5, 0.0)
                } else {
                    C64::ZERO
                };
        }
    }
    (s, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    const G: f64 = 9.81;
    const RHO: f64 = 1025.0;

    fn omega_for(nu_r: f64, r: f64) -> f64 {
        (nu_r / r * G).sqrt()
    }

    /// Damping from the near-field pressure equals the energy the far-field
    /// waves carry away: `b = ρω|C|²`.
    #[test]
    fn damping_balances_the_radiated_energy() {
        let r = 1.0;
        for &nur in &[0.2, 0.6, 1.0, 1.4] {
            let sol = Section::semicircle(r, 48)
                .heave(omega_for(nur, r), G, RHO)
                .unwrap();
            let far = RHO * sol.omega * sol.far.abs_sq();
            assert!(
                (sol.damping - far).abs() < 5e-3 * far,
                "νR {nur}: b {} vs far {far}",
                sol.damping
            );
        }
    }

    /// Sections resampled from sampled curves — the path real hulls take —
    /// with panels graded to millimetres at the waterline and keel: the
    /// same energy balance, and agreement with the exact semicircle. (A
    /// roundoff-signed self term once flipped these tiny panels' diagonal.)
    #[test]
    fn curve_sections_balance_energy() {
        let r = 0.15;
        let semi: Vec<(f64, f64)> = (0..=128)
            .map(|i| {
                let th = 0.5 * PI * i as f64 / 128.0;
                (r * th.cos(), r * th.sin())
            })
            .collect();
        let (b, t) = (0.15, 0.1875);
        let para: Vec<(f64, f64)> = (0..=64)
            .map(|i| {
                let z = t * i as f64 / 64.0;
                (b * (1.0 - (z / t).powi(2)), z)
            })
            .collect();
        for &nu in &[2.0, 4.0, 8.0] {
            let w = (nu * G).sqrt();
            let exact = Section::semicircle(r, 64).heave(w, G, RHO).unwrap();
            for (name, curve) in [("semicircle", &semi), ("parabola", &para)] {
                for n in [16, 24, 40] {
                    let sol = Section::from_curve(curve, n)
                        .unwrap()
                        .heave(w, G, RHO)
                        .unwrap();
                    let far = RHO * w * sol.far.abs_sq();
                    assert!(sol.damping > 0.0, "{name} {n} ν {nu}: b {}", sol.damping);
                    assert!(
                        (sol.damping - far).abs() < 3e-2 * far,
                        "{name} {n} ν {nu}: b {} vs {far}",
                        sol.damping
                    );
                    if name == "semicircle" {
                        assert!(
                            (sol.added_mass - exact.added_mass).abs() < 3e-2 * exact.added_mass,
                            "{n} ν {nu}: a {} vs {}",
                            sol.added_mass,
                            exact.added_mass
                        );
                        assert!(
                            (sol.damping - exact.damping).abs() < 3e-2 * exact.damping,
                            "{n} ν {nu}: b {} vs {}",
                            sol.damping,
                            exact.damping
                        );
                    }
                }
            }
        }
    }

    /// A Wigley midship section's first irregular frequency (ν ≈ 12/m for
    /// B = 0.3 m, T = 0.1875 m) is detected by its energy mismatch and
    /// bridged: the result lies between its neighbours.
    #[test]
    fn an_irregular_frequency_is_bridged() {
        let (b, t) = (0.15, 0.1875);
        let para: Vec<(f64, f64)> = (0..=64)
            .map(|i| {
                let z = t * i as f64 / 64.0;
                (b * (1.0 - (z / t).powi(2)), z)
            })
            .collect();
        let sec = Section::from_curve(&para, 24).unwrap();
        let at = |nu: f64| sec.heave((nu * G).sqrt(), G, RHO).unwrap();
        let raw = sec.heave_at((12.0 * G).sqrt(), G, RHO).unwrap();
        assert!(
            raw.energy_error() > ENERGY_TOLERANCE,
            "the raw solve should be polluted: {}",
            raw.energy_error()
        );
        let (lo, mid, hi) = (at(11.0), at(12.0), at(13.0));
        assert!(mid.interpolated && !lo.interpolated && !hi.interpolated);
        let between =
            |a: f64, m: f64, c: f64| m > a.min(c) - 0.05 * a.abs() && m < a.max(c) + 0.05 * a.abs();
        assert!(
            between(lo.added_mass, mid.added_mass, hi.added_mass),
            "{} {} {}",
            lo.added_mass,
            mid.added_mass,
            hi.added_mass
        );
        let bw = |s: &HeaveSolution| s.damping / s.omega;
        assert!(
            between(bw(&lo), bw(&mid), bw(&hi)),
            "{} {} {}",
            bw(&lo),
            bw(&mid),
            bw(&hi)
        );
    }

    #[test]
    fn heave_coefficients_converge_with_panels() {
        let r = 1.0;
        let w = omega_for(0.8, r);
        let coarse = Section::semicircle(r, 32).heave(w, G, RHO).unwrap();
        let fine = Section::semicircle(r, 128).heave(w, G, RHO).unwrap();
        let m = RHO * PI * r * r / 2.0;
        eprintln!(
            "νR 0.8: a/m {} → {}, b/(mω) {} → {}",
            coarse.added_mass / m,
            fine.added_mass / m,
            coarse.damping / (m * w),
            fine.damping / (m * w)
        );
        assert!((coarse.added_mass - fine.added_mass).abs() < 1e-2 * m);
        assert!((coarse.damping - fine.damping).abs() < 1e-2 * m * w);
    }

    /// The Haskind relation reproduces the diffraction force of the solved
    /// diffraction problem (zero speed, where both are at the wave's own
    /// frequency): a check on the sign chain and on reciprocity.
    #[test]
    fn haskind_matches_the_direct_diffraction_force() {
        for sec in [
            Section::semicircle(1.0, 64),
            Section::rectangle(1.0, 0.6, 64),
        ] {
            for &nu in &[0.3, 0.9, 1.3] {
                let w = (nu * G).sqrt();
                let rad = sec.heave(w, G, RHO).unwrap();
                let hask = rad.diffraction(nu, PI, G, RHO);
                let direct = sec.diffraction_direct(w, G, RHO).unwrap();
                assert!(
                    (hask - direct).abs() < 1e-2 * direct.abs(),
                    "ν {nu}: {hask:?} vs {direct:?}"
                );
            }
        }
    }

    /// The classical semicircle curves (Ursell's multipoles; Frank's and
    /// Vugts' computations) put the heave added-mass coefficient
    /// `a/(ρπR²/2)` near its minimum, about 0.6, around νR = 1, with a
    /// wave-amplitude ratio near 0.8 (this solver: 0.61, 0.79). The bounds
    /// are loose, recalled values rather than a tabulation — a guard against
    /// a sign or factor error; the energy and Haskind checks above are the
    /// sharp ones.
    #[test]
    fn semicircle_near_published_values() {
        let r = 1.0;
        let sol = Section::semicircle(r, 64)
            .heave(omega_for(1.0, r), G, RHO)
            .unwrap();
        let ca = sol.added_mass / (RHO * PI * r * r / 2.0);
        let amp = sol.wave_ratio(G);
        eprintln!("νR 1: Ca {ca}, |ζ/η| {amp}");
        assert!(ca > 0.5 && ca < 1.1, "Ca {ca}");
        assert!(amp > 0.4 && amp < 1.0, "|ζ/η| {amp}");
    }
}

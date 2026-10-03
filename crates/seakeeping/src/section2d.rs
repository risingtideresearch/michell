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
//! A plain source method fails at the **irregular frequencies** of the
//! interior (the sloshing modes of the fluid the section would enclose, the
//! first near `ν ≈ (π/B) coth(πT/B)` for a box of beam `B`, draft `T`), and
//! near them its source strengths carry a spurious interior mode that
//! nearly cancels in the section's own far field but not at other
//! wavenumbers — which is where a Kochin function reads them. So the
//! interior waterplane is closed by a **lid** of sources on which the
//! interior's vertical velocity is held at zero (Ohmatsu's remedy, in its
//! rigid-lid form): the interior problem then has no eigenvalues, the
//! sources are unique at every frequency, and the exterior solution is the
//! plain method's wherever that one is sound (the two converge together on
//! a semicircle). The energy check — damping against the radiated energy —
//! still runs on every solve and bridges a failure by interpolation, as a
//! safeguard.

use crate::green::{far_factor, wave_part};
use crate::linalg::{solve, solve_many};
use hullgeom::quadrature::gauss_legendre;
use hullgeom::C64;
use std::f64::consts::PI;

/// A section's starboard half, as panels from the waterline to the keel.
#[derive(Debug, Clone)]
pub struct Section {
    /// Panel end points `(y, z)`, `z` up, from the waterline to the keel.
    nodes: Vec<[f64; 2]>,
    /// Whether to close the interior waterplane with a lid (see
    /// `Section::lid`); on by default.
    lid: bool,
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
    /// Every source panel: the body's, then the lid's (see `Section::lid`).
    sources: Vec<Panel>,
    panels: Vec<Panel>,
    /// `ψ` at the panel midpoints (starboard; port is equal).
    psi: Vec<C64>,
    /// Source strengths per unit heave velocity (outward volume flux per
    /// unit area, starboard; port is equal).
    sigma: Vec<C64>,
    /// Source strengths of the diffraction problem for a unit-amplitude
    /// incident wave, when solved ([`Section::heave_with_diffraction`]);
    /// empty otherwise.
    sigma_d: Vec<C64>,
}

/// An incident wave for the diffraction sources: wavenumber `k`, heading.
#[derive(Debug, Clone, Copy)]
struct Incident {
    k: f64,
    heading: f64,
    gravity: f64,
}

/// A solved panel system (see `Section::solve_parity`).
struct Solved {
    panels: Vec<Panel>,
    sources: Vec<Panel>,
    s_body: Vec<C64>,
    strengths: Vec<Vec<C64>>,
    nu: f64,
}

impl Solved {
    /// The potential at the body panels' midpoints of the strengths `sigma`.
    fn body_potential(&self, sigma: &[C64]) -> Vec<C64> {
        let size = self.sources.len();
        (0..self.panels.len())
            .map(|i| {
                (0..size).fold(C64::ZERO, |acc, j| {
                    acc + self.s_body[i * size + j] * sigma[j]
                })
            })
            .collect()
    }
}

/// The roll normal `y n_z − z n_y` at a panel's midpoint: the normal
/// velocity of unit roll about the waterline centre, right-handed about +x.
fn roll_normal(p: &Panel) -> f64 {
    p.mid[0] * p.n[1] - p.mid[1] * p.n[0]
}

/// A section's sway and roll radiation solution at one frequency, per unit
/// velocity of each mode; index 0 sway (+y), 1 roll (about the waterline
/// centre, +y side up).
#[derive(Debug, Clone)]
pub struct LateralSolution {
    pub omega: f64,
    /// Added mass `a_ij` per unit length: force in mode `i` from mode `j`.
    pub added_mass: [[f64; 2]; 2],
    /// Damping `b_ij` per unit length.
    pub damping: [[f64; 2]; 2],
    /// Far-field amplitude on the `+y` side per unit velocity of each mode
    /// (`−` of it on the `−y` side): `ψ → ±C e^{νz} e^{iν|y|}`.
    pub far: [C64; 2],
    /// Whether this solution was bridged across a failed energy check.
    pub interpolated: bool,
    density: f64,
    sources: Vec<Panel>,
    panels: Vec<Panel>,
    psi: [Vec<C64>; 2],
    sigma: [Vec<C64>; 2],
    sigma_d: Vec<C64>,
    /// The diffraction potential at the body panels, when solved.
    psi_d: Vec<C64>,
}

impl LateralSolution {
    /// The worst relative disagreement between either mode's damping and the
    /// energy its waves carry, `|b_jj/(ρω|C_j|²) − 1|`.
    pub fn energy_error(&self) -> f64 {
        (0..2)
            .map(|j| {
                let far = self.density * self.omega * self.far[j].abs_sq();
                if far > 0.0 {
                    (self.damping[j][j] / far - 1.0).abs()
                } else if self.damping[j][j].abs() > 0.0 {
                    f64::INFINITY
                } else {
                    0.0
                }
            })
            .fold(0.0, f64::max)
    }

    /// The antisymmetric diffraction force (sway, roll) per unit length and
    /// unit wave amplitude by the Haskind relation on these radiation
    /// potentials:
    /// `f_j = iρgk ∫ ψ_j e^{kz} (n_z sin(k y sin β) + sin β n_y cos(k y sin β)) ds`,
    /// exact within strip theory at the wave's own frequency.
    pub fn diffraction(&self, k: f64, heading: f64, gravity: f64, density: f64) -> [C64; 2] {
        let sb = heading.sin();
        let mut out = [C64::ZERO; 2];
        for (j, o) in out.iter_mut().enumerate() {
            let sum = self
                .panels
                .iter()
                .zip(&self.psi[j])
                .fold(C64::ZERO, |acc, (p, &f)| {
                    let (y, z) = (p.mid[0], p.mid[1]);
                    let arg = k * y * sb;
                    let w = (k * z).exp() * (p.n[1] * arg.sin() + sb * p.n[0] * arg.cos());
                    acc + f.scale(2.0 * w * p.len)
                });
            *o = C64::new(0.0, density * gravity * k) * sum;
        }
        out
    }

    /// The Froude–Krylov force (sway, roll) per unit length and unit wave
    /// amplitude: `−∫ p n_j ds` with `p = ρg e^{kz} e^{iky sin β}`, whose part
    /// odd in `y` is `iρg e^{kz} sin(k y sin β)`.
    pub fn froude_krylov(&self, k: f64, heading: f64, gravity: f64, density: f64) -> [C64; 2] {
        let sb = heading.sin();
        let mut out = [C64::ZERO; 2];
        for (j, o) in out.iter_mut().enumerate() {
            let sum: f64 = self
                .panels
                .iter()
                .map(|p| {
                    let (y, z) = (p.mid[0], p.mid[1]);
                    2.0 * (k * z).exp() * (k * y * sb).sin() * [p.n[0], roll_normal(p)][j] * p.len
                })
                .sum();
            *o = C64::new(0.0, -density * gravity * sum);
        }
        out
    }
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
        Some(Section { nodes, lid: true })
    }

    /// A section from its panel end points `(y, z)`, `z` up, from the
    /// waterline to the keel.
    pub fn from_nodes(nodes: Vec<[f64; 2]>) -> Section {
        Section { nodes, lid: true }
    }

    /// This section without the interior lid: the plain Frank method, with
    /// its irregular frequencies (for comparison).
    pub fn without_lid(mut self) -> Section {
        self.lid = false;
        self
    }

    /// A semicircle of radius `r` (a heaving half-immersed cylinder).
    pub fn semicircle(r: f64, panels: usize) -> Section {
        let nodes = (0..=panels)
            .map(|i| {
                let th = 0.5 * PI * i as f64 / panels as f64;
                [r * th.cos(), -r * th.sin()]
            })
            .collect();
        Section { nodes, lid: true }
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
        Section { nodes, lid: true }
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

    /// The interior waterplane ("lid"), from the waterline to the
    /// centreplane: Ohmatsu's (1975) remedy for the irregular frequencies.
    /// Sources on it, with the interior's vertical velocity held at zero at
    /// its midpoints, leave the interior problem with the hull's potential
    /// on its curve and a rigid lid — no eigenvalues — so the source
    /// representation is unique at every frequency, while the exterior
    /// solution is unchanged. (Holding the lid's *potential* at zero instead
    /// also removes them, but clashes with the hull's potential at the
    /// waterline corner and converges far too slowly.)
    fn lid(&self) -> Vec<Panel> {
        if !self.lid {
            return Vec::new();
        }
        let b = self.half_beam();
        let m = (self.nodes.len() / 3).max(4);
        (0..m)
            .map(|i| {
                let ya = b * 0.5 * (1.0 + (PI * i as f64 / m as f64).cos());
                let yb = b * 0.5 * (1.0 + (PI * (i + 1) as f64 / m as f64).cos());
                Panel {
                    a: [ya, 0.0],
                    b: [yb, 0.0],
                    mid: [0.5 * (ya + yb), 0.0],
                    len: ya - yb,
                    n: [0.0, 1.0],
                }
            })
            .collect()
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
        self.heave_bridged(omega, gravity, density, None)
    }

    /// [`Section::heave`] together with the sources of the diffraction
    /// problem for a unit-amplitude incident wave of wavenumber `k` at
    /// `heading` — its part symmetric about the centreplane, the body
    /// condition `∂φ_D/∂n = −∂φ_I/∂n` with the section's free surface at
    /// the radiation problem's frequency `omega` (strip theory at speed: the
    /// encounter frequency). The far field of both is what
    /// [`HeaveSolution::source_spectrum`] integrates.
    pub fn heave_with_diffraction(
        &self,
        omega: f64,
        k: f64,
        heading: f64,
        gravity: f64,
        density: f64,
    ) -> Option<HeaveSolution> {
        self.heave_bridged(
            omega,
            gravity,
            density,
            Some(Incident {
                k,
                heading,
                gravity,
            }),
        )
    }

    fn heave_bridged(
        &self,
        omega: f64,
        gravity: f64,
        density: f64,
        inc: Option<Incident>,
    ) -> Option<HeaveSolution> {
        let sol = self.heave_at(omega, gravity, density, inc)?;
        if sol.energy_error() <= ENERGY_TOLERANCE {
            return Some(sol);
        }
        for eps in [0.02, 0.04, 0.08] {
            let lo = self.heave_at(omega * (1.0 - eps), gravity, density, inc)?;
            let hi = self.heave_at(omega * (1.0 + eps), gravity, density, inc)?;
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
                    sigma: lo
                        .sigma
                        .iter()
                        .zip(&hi.sigma)
                        .map(|(&a, &b)| mix(a, b))
                        .collect(),
                    sigma_d: lo
                        .sigma_d
                        .iter()
                        .zip(&hi.sigma_d)
                        .map(|(&a, &b)| mix(a, b))
                        .collect(),
                    sources: lo.sources,
                    panels: lo.panels,
                    interpolated: true,
                    density,
                });
            }
        }
        Some(sol)
    }

    /// The panel system of one parity (`1` symmetric, `−1` antisymmetric
    /// about the centreplane), lid included, solved for each of `rhs` (the
    /// normal velocity at each body panel's midpoint): the source panels
    /// (body then lid), the body's potential matrix, and the strengths.
    fn solve_parity(
        &self,
        omega: f64,
        gravity: f64,
        parity: f64,
        rhs: Vec<Vec<C64>>,
    ) -> Option<Solved> {
        let panels = self.panels();
        let lid = self.lid();
        let nu = omega * omega / gravity;
        let n = panels.len();
        let m = lid.len();
        let sources: Vec<Panel> = panels.iter().chain(&lid).copied().collect();
        let size = n + m;
        // Rows: the body's normal velocity, then the lid's zero vertical
        // velocity. Columns: body sources, then lid sources. `s_body` keeps
        // the body's potentials for the pressure.
        let mut d = vec![C64::ZERO; size * size];
        let mut s_body = vec![C64::ZERO; n * size];
        for (i, pi) in panels.iter().enumerate() {
            for (j, pj) in sources.iter().enumerate() {
                let (pot, dn) = pair(pi.mid, pi.n, pj, nu, i == j, parity);
                s_body[i * size + j] = pot;
                d[i * size + j] = dn
                    + if i == j {
                        C64::new(0.5, 0.0)
                    } else {
                        C64::ZERO
                    };
            }
        }
        for (l, pl) in lid.iter().enumerate() {
            for (j, pj) in sources.iter().enumerate() {
                // ∂φ/∂z from below: a lid panel's own layer and its
                // coincident image each give −½ there.
                let (_, dz) = pair(pl.mid, pl.n, pj, nu, j == n + l, parity);
                d[(n + l) * size + j] = dz - if j == n + l { C64::ONE } else { C64::ZERO };
            }
        }
        let rhs = rhs
            .into_iter()
            .map(|mut v| {
                v.resize(size, C64::ZERO);
                v
            })
            .collect();
        let strengths = solve_many(d, rhs)?;
        Some(Solved {
            panels,
            sources,
            s_body,
            strengths,
            nu,
        })
    }

    fn heave_at(
        &self,
        omega: f64,
        gravity: f64,
        density: f64,
        inc: Option<Incident>,
    ) -> Option<HeaveSolution> {
        let panels = self.panels();
        let mut rhs = vec![panels
            .iter()
            .map(|p| C64::new(p.n[1], 0.0))
            .collect::<Vec<_>>()];
        if let Some(inc) = inc {
            // φ_I = −i(g/ω₀) e^{kz} e^{iky sin β}; its normal derivative's part
            // even in y is k φ_I(n_z cos(k y sin β) − sin β n_y sin(k y sin β)),
            // and ∂φ_D/∂n = −∂φ_I/∂n = iω₀ e^{kz}(…).
            let w0 = (inc.gravity * inc.k).sqrt();
            let sb = inc.heading.sin();
            rhs.push(
                panels
                    .iter()
                    .map(|p| {
                        let arg = inc.k * p.mid[0] * sb;
                        let w = (inc.k * p.mid[1]).exp()
                            * (p.n[1] * arg.cos() - sb * p.n[0] * arg.sin());
                        C64::new(0.0, w0 * w)
                    })
                    .collect(),
            );
        }
        let mut solved = self.solve_parity(omega, gravity, 1.0, rhs)?;
        let sigma_d = if solved.strengths.len() > 1 {
            solved.strengths.pop().unwrap()
        } else {
            Vec::new()
        };
        let sigma = solved.strengths.pop().unwrap();
        let psi = solved.body_potential(&sigma);
        let force: C64 = solved
            .panels
            .iter()
            .zip(&psi)
            .fold(C64::ZERO, |acc, (p, &f)| {
                acc + f.scale(2.0 * p.n[1] * p.len)
            });
        let far = far_amplitude(&solved.sources, &sigma, solved.nu);
        Some(HeaveSolution {
            omega,
            added_mass: -density * force.re,
            damping: -omega * density * force.im,
            far,
            sources: solved.sources,
            panels: solved.panels,
            psi,
            sigma,
            sigma_d,
            interpolated: false,
            density,
        })
    }

    /// The **sway and roll** radiation problems at frequency `omega` — the
    /// antisymmetric modes, their mirror sources of opposite sign — and,
    /// with `incident` (`(k, heading)`), the antisymmetric part of the
    /// diffraction problem. Roll is about the section's waterline centre,
    /// right-handed about +x (the `+y` side rises). Irregular frequencies
    /// are removed by the same lid; the energy check still runs on both
    /// modes and bridges a failure.
    pub fn lateral(
        &self,
        omega: f64,
        gravity: f64,
        density: f64,
        incident: Option<(f64, f64)>,
    ) -> Option<LateralSolution> {
        let inc = incident.map(|(k, heading)| Incident {
            k,
            heading,
            gravity,
        });
        let sol = self.lateral_at(omega, gravity, density, inc)?;
        if sol.energy_error() <= ENERGY_TOLERANCE {
            return Some(sol);
        }
        for eps in [0.02, 0.04, 0.08] {
            let lo = self.lateral_at(omega * (1.0 - eps), gravity, density, inc)?;
            let hi = self.lateral_at(omega * (1.0 + eps), gravity, density, inc)?;
            if lo.energy_error() <= ENERGY_TOLERANCE && hi.energy_error() <= ENERGY_TOLERANCE {
                let mix = |a: C64, b: C64| (a + b).scale(0.5);
                let mixv = |a: &[C64], b: &[C64]| {
                    a.iter()
                        .zip(b)
                        .map(|(&x, &y)| mix(x, y))
                        .collect::<Vec<_>>()
                };
                let mut added_mass = [[0.0; 2]; 2];
                let mut damping = [[0.0; 2]; 2];
                for i in 0..2 {
                    for j in 0..2 {
                        added_mass[i][j] = 0.5 * (lo.added_mass[i][j] + hi.added_mass[i][j]);
                        damping[i][j] = 0.5
                            * (lo.damping[i][j] / lo.omega + hi.damping[i][j] / hi.omega)
                            * omega;
                    }
                }
                return Some(LateralSolution {
                    omega,
                    added_mass,
                    damping,
                    far: [mix(lo.far[0], hi.far[0]), mix(lo.far[1], hi.far[1])],
                    psi: [mixv(&lo.psi[0], &hi.psi[0]), mixv(&lo.psi[1], &hi.psi[1])],
                    sigma: [
                        mixv(&lo.sigma[0], &hi.sigma[0]),
                        mixv(&lo.sigma[1], &hi.sigma[1]),
                    ],
                    sigma_d: mixv(&lo.sigma_d, &hi.sigma_d),
                    psi_d: mixv(&lo.psi_d, &hi.psi_d),
                    sources: lo.sources,
                    panels: lo.panels,
                    interpolated: true,
                    density,
                });
            }
        }
        Some(sol)
    }

    fn lateral_at(
        &self,
        omega: f64,
        gravity: f64,
        density: f64,
        inc: Option<Incident>,
    ) -> Option<LateralSolution> {
        let panels = self.panels();
        let mut rhs = vec![
            panels
                .iter()
                .map(|p| C64::new(p.n[0], 0.0))
                .collect::<Vec<_>>(),
            panels
                .iter()
                .map(|p| C64::new(roll_normal(p), 0.0))
                .collect::<Vec<_>>(),
        ];
        if let Some(inc) = inc {
            // The part of −∂φ_I/∂n odd in y:
            // −ω₀ e^{kz}(n_z sin(k y sin β) + sin β n_y cos(k y sin β)).
            let w0 = (inc.gravity * inc.k).sqrt();
            let sb = inc.heading.sin();
            rhs.push(
                panels
                    .iter()
                    .map(|p| {
                        let arg = inc.k * p.mid[0] * sb;
                        let w = (inc.k * p.mid[1]).exp()
                            * (p.n[1] * arg.sin() + sb * p.n[0] * arg.cos());
                        C64::new(-w0 * w, 0.0)
                    })
                    .collect(),
            );
        }
        let mut solved = self.solve_parity(omega, gravity, -1.0, rhs)?;
        let sigma_d = if solved.strengths.len() > 2 {
            solved.strengths.pop().unwrap()
        } else {
            Vec::new()
        };
        let sigma4 = solved.strengths.pop().unwrap();
        let sigma2 = solved.strengths.pop().unwrap();
        let psi = [
            solved.body_potential(&sigma2),
            solved.body_potential(&sigma4),
        ];
        let psi_d = if sigma_d.is_empty() {
            Vec::new()
        } else {
            solved.body_potential(&sigma_d)
        };
        let normals = |p: &Panel| [p.n[0], roll_normal(p)];
        let mut added_mass = [[0.0; 2]; 2];
        let mut damping = [[0.0; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                // Force in mode i from unit velocity in mode j (both sides).
                let f = solved
                    .panels
                    .iter()
                    .zip(&psi[j])
                    .fold(C64::ZERO, |acc, (p, &v)| {
                        acc + v.scale(2.0 * normals(p)[i] * p.len)
                    });
                added_mass[i][j] = -density * f.re;
                damping[i][j] = -omega * density * f.im;
            }
        }
        let far = [
            far_amplitude_odd(&solved.sources, &sigma2, solved.nu),
            far_amplitude_odd(&solved.sources, &sigma4, solved.nu),
        ];
        Some(LateralSolution {
            omega,
            added_mass,
            damping,
            far,
            psi,
            sigma: [sigma2, sigma4],
            sigma_d,
            psi_d,
            sources: solved.sources,
            panels: solved.panels,
            interpolated: false,
            density,
        })
    }

    /// The antisymmetric diffraction force (sway, roll) per unit length and
    /// unit wave amplitude, by solving the diffraction problem at the wave's
    /// own frequency — the reference for [`LateralSolution::diffraction`].
    pub fn lateral_diffraction_direct(
        &self,
        omega: f64,
        heading: f64,
        gravity: f64,
        density: f64,
    ) -> Option<[C64; 2]> {
        let k = omega * omega / gravity;
        let sol = self.lateral_at(
            omega,
            gravity,
            density,
            Some(Incident {
                k,
                heading,
                gravity,
            }),
        )?;
        let mut out = [C64::ZERO; 2];
        for (i, o) in out.iter_mut().enumerate() {
            let integral = sol
                .panels
                .iter()
                .zip(&sol.psi_d)
                .fold(C64::ZERO, |acc, (p, &f)| {
                    acc + f.scale(2.0 * [p.n[0], roll_normal(p)][i] * p.len)
                });
            *o = C64::new(0.0, -omega * density) * integral;
        }
        Some(out)
    }

    /// Heave added mass per unit length at **infinite frequency**, where the
    /// free surface holds `φ = 0` and the source is `ln r − ln r₁`: the limit
    /// the frequency-dependent added mass tends to, and — for a semicircle,
    /// `ρπR²/2` exactly — an analytic check on the panel integrals.
    pub fn added_mass_infinite(&self, density: f64) -> Option<f64> {
        let panels = self.panels();
        let n = panels.len();
        let mirror = |q: [f64; 2]| [-q[0], q[1]];
        let image = |q: [f64; 2]| [q[0], -q[1]];
        let inv2pi = 1.0 / (2.0 * PI);
        let mut s = vec![C64::ZERO; n * n];
        let mut d = vec![C64::ZERO; n * n];
        for (i, pi) in panels.iter().enumerate() {
            for (j, pj) in panels.iter().enumerate() {
                let (mut val, mut grad) = (0.0, [0.0; 2]);
                for (m, (a, b), sign) in [
                    (0, (pj.a, pj.b), 1.0),
                    (1, (mirror(pj.b), mirror(pj.a)), 1.0),
                    (2, (image(pj.a), image(pj.b)), -1.0),
                    (3, (image(mirror(pj.b)), image(mirror(pj.a))), -1.0),
                ] {
                    let (v, g) = log_panel(pi.mid, a, b, m == 0 && i == j);
                    val += sign * v;
                    grad[0] += sign * g[0];
                    grad[1] += sign * g[1];
                }
                s[i * n + j] = C64::new(val * inv2pi, 0.0);
                let dn = (grad[0] * pi.n[0] + grad[1] * pi.n[1]) * inv2pi
                    + if i == j { 0.5 } else { 0.0 };
                d[i * n + j] = C64::new(dn, 0.0);
            }
        }
        let rhs: Vec<C64> = panels.iter().map(|p| C64::new(p.n[1], 0.0)).collect();
        let sigma = solve(d, rhs)?;
        let psi = apply(&s, &sigma, n);
        let force: f64 = panels
            .iter()
            .zip(&psi)
            .map(|(p, f)| f.re * 2.0 * p.n[1] * p.len)
            .sum();
        Some(-density * force)
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
    /// The Froude–Krylov heave force per unit length and unit wave
    /// amplitude, `−∫ p n_z ds` with `p = ρg e^{kz} e^{iky sin β}` (its part
    /// even in `y`, `ρg e^{kz} cos(k y sin β)`, both sides) — the section's
    /// own form of what [`crate::froude_krylov`] gives a whole hull in
    /// closed form.
    pub fn froude_krylov(&self, k: f64, heading: f64, gravity: f64, density: f64) -> C64 {
        let sb = heading.sin();
        let sum: f64 = self
            .panels
            .iter()
            .map(|p| 2.0 * (k * p.mid[1]).exp() * (k * p.mid[0] * sb).cos() * p.n[1] * p.len)
            .sum();
        C64::new(-density * gravity * sum, 0.0)
    }

    /// Drop the diffraction sources (a radiation-only far field).
    pub fn clear_diffraction(&mut self) {
        self.sigma_d.clear();
    }

    /// The section's share of a three-dimensional Kochin function: the
    /// section's sources, heaving at velocity `v` and diffracting an
    /// incident wave of complex amplitude `d`, summed against the far-field
    /// wave of wavenumber `k` travelling at angle `θ` (`sin_theta`),
    ///
    /// ```text
    /// Σ ∫ (v σ + d σ_D) e^{kz} 2 cos(k y sin θ) ds      (both sides)
    /// ```
    ///
    /// per unit length of hull; the longitudinal phase `e^{−ik x cos θ}` is
    /// the caller's. The diffraction part needs the sources of
    /// [`Section::heave_with_diffraction`].
    pub fn source_spectrum(&self, k: f64, sin_theta: f64, v: C64, d: C64) -> C64 {
        const GAUSS2: [f64; 2] = [0.211_324_865_405_187_1, 0.788_675_134_594_812_9];
        let mut total = C64::ZERO;
        for (j, p) in self.sources.iter().enumerate() {
            let strength = self.sigma[j] * v + self.sigma_d.get(j).map_or(C64::ZERO, |&s| s * d);
            let mut w = 0.0;
            for u in GAUSS2 {
                let (y, z) = (
                    p.a[0] + u * (p.b[0] - p.a[0]),
                    p.a[1] + u * (p.b[1] - p.a[1]),
                );
                w += (k * z).exp() * (k * y * sin_theta).cos();
            }
            total = total + strength.scale(w * p.len);
        }
        total
    }

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

/// The far-field amplitude of antisymmetric sources on the `+y` side,
/// `C = −2 Σ_j σ_j ∫_j e^{νζ} sin νη ds`.
fn far_amplitude_odd(panels: &[Panel], sigma: &[C64], nu: f64) -> C64 {
    let (gx, gw) = gauss_legendre(8);
    let mut sum = C64::ZERO;
    for (p, &sg) in panels.iter().zip(sigma) {
        let mut f = 0.0;
        for (&t, &w) in gx.iter().zip(&gw) {
            let u = 0.5 * (1.0 + t);
            let q = [
                p.a[0] + u * (p.b[0] - p.a[0]),
                p.a[1] + u * (p.b[1] - p.a[1]),
            ];
            f += 0.5 * w * p.len * (nu * q[1]).exp() * (nu * q[0]).sin();
        }
        sum = sum + sg.scale(f);
    }
    sum.scale(-2.0)
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
    let mut s = vec![C64::ZERO; n * n];
    let mut d = vec![C64::ZERO; n * n];
    for (i, pi) in panels.iter().enumerate() {
        for (j, pj) in panels.iter().enumerate() {
            let (pot, dn) = pair(pi.mid, pi.n, pj, nu, i == j, 1.0);
            s[i * n + j] = pot;
            d[i * n + j] = dn
                + if i == j {
                    C64::new(0.5, 0.0)
                } else {
                    C64::ZERO
                };
        }
    }
    (s, d)
}

/// The potential and the normal derivative (along `n`) at the field point
/// `p` of a unit source density on panel `pj` and its port mirror, both
/// normalised `1/2π` — without the self term `½` of a panel's own normal
/// derivative. `on_panel` says `p` is `pj`'s own midpoint; a panel on the
/// waterline is its own image, so there the image integral is taken on the
/// panel too.
fn pair(p: [f64; 2], n: [f64; 2], pj: &Panel, nu: f64, on_panel: bool, parity: f64) -> (C64, C64) {
    let (gx, gw) = gauss_legendre(8);
    let inv2pi = 1.0 / (2.0 * PI);
    let mirror = |q: [f64; 2]| [-q[0], q[1]];
    let image = |q: [f64; 2]| [q[0], -q[1]];
    let on_waterline = pj.a[1] == 0.0 && pj.b[1] == 0.0;
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
        let self_line = on_panel && (m == 0 || (m == 2 && on_waterline));
        let (v, g) = log_panel(p, a, b, self_line);
        // The mirror (m = 1, 3) carries the parity's sign.
        let sign = if m % 2 == 1 { parity } else { 1.0 };
        val += sign * v;
        grad[0] += sign * g[0];
        grad[1] += sign * g[1];
        if m >= 2 {
            image_logs += sign * v;
        }
    }
    // The wave part. Its value and horizontal gradient are bounded; its
    // vertical gradient `−2ν Re P` goes as `2ν ln(ν r₁)` where the panel
    // meets the waterline, so that logarithm is subtracted here and added
    // back in closed form from `∫ ln r₁`. Panels are subdivided in
    // proportion to their length over their distance.
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
            for (qy, sign) in [(q[0], 1.0), (-q[0], parity)] {
                let (dy, zs) = (p[0] - qy, p[1] + q[1]);
                let (v, gy, gz) = wave_part(nu, dy, zs);
                let r1 = (dy * dy + zs * zs).sqrt();
                let gz_reg = if r1 > 0.0 {
                    gz - C64::new(2.0 * nu * (nu * r1).ln(), 0.0)
                } else {
                    gz
                };
                wv = wv + v.scale(ds * sign);
                wg[0] = wg[0] + gy.scale(ds * sign);
                wg[1] = wg[1] + gz_reg.scale(ds * sign);
            }
        }
    }
    // ∫ 2ν ln(ν r₁) over the panel and its mirror (the mirror signed).
    wg[1] = wg[1]
        + C64::new(
            2.0 * nu * (image_logs + (1.0 + parity) * pj.len * nu.ln()),
            0.0,
        );
    let pot = (C64::new(val, 0.0) + wv).scale(inv2pi);
    let dn =
        (C64::new(grad[0] * n[0] + grad[1] * n[1], 0.0) + wg[0].scale(n[0]) + wg[1].scale(n[1]))
            .scale(inv2pi);
    (pot, dn)
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
    /// B = 0.3 m, T = 0.1875 m). Without the lid, the plain Frank solve is
    /// polluted there — its energy check fails — and is bridged; with the
    /// lid (the default) there is nothing to bridge, and its sources, whose
    /// far field at other wavenumbers the Kochin function reads, run
    /// smoothly through: the result lies between its neighbours.
    #[test]
    fn an_irregular_frequency_is_removed_by_the_lid() {
        let (b, t) = (0.15, 0.1875);
        let para: Vec<(f64, f64)> = (0..=64)
            .map(|i| {
                let z = t * i as f64 / 64.0;
                (b * (1.0 - (z / t).powi(2)), z)
            })
            .collect();
        let w = |nu: f64| (nu * G).sqrt();
        let bare = Section::from_curve(&para, 24).unwrap().without_lid();
        let raw = bare.heave_at(w(12.0), G, RHO, None).unwrap();
        assert!(
            raw.energy_error() > ENERGY_TOLERANCE,
            "the plain solve should be polluted: {}",
            raw.energy_error()
        );
        assert!(bare.heave(w(12.0), G, RHO).unwrap().interpolated);
        let sec = Section::from_curve(&para, 24).unwrap();
        let at = |nu: f64| sec.heave_with_diffraction(w(nu), 2.0, PI, G, RHO).unwrap();
        let (lo, mid, hi) = (at(11.0), at(12.0), at(13.0));
        assert!(!lo.interpolated && !mid.interpolated && !hi.interpolated);
        assert!(mid.energy_error() < 0.03, "{}", mid.energy_error());
        let between =
            |a: f64, m: f64, c: f64| m > a.min(c) - 0.02 * a.abs() && m < a.max(c) + 0.02 * a.abs();
        let far_d =
            |s: &HeaveSolution, nu: f64| s.source_spectrum(nu, 1.0, C64::ZERO, C64::ONE).abs();
        for (a, m, c) in [
            (lo.added_mass, mid.added_mass, hi.added_mass),
            (
                lo.damping / lo.omega,
                mid.damping / mid.omega,
                hi.damping / hi.omega,
            ),
            (far_d(&lo, 11.0), far_d(&mid, 12.0), far_d(&hi, 13.0)),
        ] {
            assert!(between(a, m, c), "{a} {m} {c}");
        }
    }

    /// A semicircle rolling about its centre moves no water (its normal is
    /// radial): zero roll added mass, damping and coupling.
    #[test]
    fn a_rolling_semicircle_moves_no_water() {
        let r = 1.0;
        let sol = Section::semicircle(r, 48)
            .lateral(omega_for(1.0, r), G, RHO, None)
            .unwrap();
        let m = RHO * PI * r * r / 2.0;
        for (i, j) in [(1, 1), (0, 1), (1, 0)] {
            assert!(
                sol.added_mass[i][j].abs() < 1e-3 * m * r,
                "a{i}{j} {}",
                sol.added_mass[i][j]
            );
            assert!(
                sol.damping[i][j].abs() < 1e-3 * m * r * 3.0,
                "b{i}{j} {}",
                sol.damping[i][j]
            );
        }
        assert!(sol.added_mass[0][0] > 0.0 && sol.damping[0][0] > 0.0);
    }

    /// Sway and roll of a box: the coefficient matrices are symmetric
    /// (reciprocity), each mode's damping is the energy its waves carry,
    /// the cross damping is the far fields' product, and they converge.
    #[test]
    fn lateral_coefficients_are_reciprocal_and_balance_energy() {
        for nu in [0.4, 1.0, 1.8] {
            let w = (nu * G).sqrt();
            let sol = Section::rectangle(1.0, 0.8, 64)
                .lateral(w, G, RHO, None)
                .unwrap();
            let (a, b) = (sol.added_mass, sol.damping);
            assert!(
                (a[0][1] - a[1][0]).abs() < 1e-2 * a[0][0].abs(),
                "ν {nu}: a24 {} a42 {}",
                a[0][1],
                a[1][0]
            );
            assert!(
                (b[0][1] - b[1][0]).abs() < 1e-2 * b[0][0].abs(),
                "ν {nu}: b24 {} b42 {}",
                b[0][1],
                b[1][0]
            );
            for j in 0..2 {
                let far = RHO * w * sol.far[j].abs_sq();
                assert!(
                    (b[j][j] - far).abs() < 2e-2 * far,
                    "ν {nu} mode {j}: b {} vs far {far}",
                    b[j][j]
                );
            }
            let cross = RHO * w * (sol.far[0] * sol.far[1].conj()).re;
            assert!(
                (b[0][1] - cross).abs() < 2e-2 * b[0][0].abs().max(b[1][1].abs()),
                "ν {nu}: b24 {} vs {cross}",
                b[0][1]
            );
            let fine = Section::rectangle(1.0, 0.8, 128)
                .lateral(w, G, RHO, None)
                .unwrap();
            assert!(
                (fine.added_mass[0][0] - a[0][0]).abs() < 3e-2 * a[0][0],
                "ν {nu}: a22 {} → {}",
                a[0][0],
                fine.added_mass[0][0]
            );
        }
    }

    /// The Haskind relation reproduces the solved antisymmetric diffraction
    /// force, in beam and oblique seas.
    #[test]
    fn lateral_haskind_matches_the_direct_diffraction_force() {
        for sec in [
            Section::semicircle(1.0, 64),
            Section::rectangle(1.0, 0.6, 64),
        ] {
            for (nu, heading) in [(0.4, 0.5 * PI), (1.0, 0.5 * PI), (0.8, 2.3)] {
                let w = (nu * G).sqrt();
                let rad = sec.lateral(w, G, RHO, None).unwrap();
                let hask = rad.diffraction(nu, heading, G, RHO);
                let direct = sec.lateral_diffraction_direct(w, heading, G, RHO).unwrap();
                for j in 0..2 {
                    let scale = direct[0].abs().max(direct[1].abs());
                    if direct[j].abs() < 1e-3 * scale {
                        continue;
                    }
                    assert!(
                        (hask[j] - direct[j]).abs() < 2e-2 * direct[j].abs(),
                        "ν {nu} β {heading} mode {j}: {:?} vs {:?}",
                        hask[j],
                        direct[j]
                    );
                }
            }
        }
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

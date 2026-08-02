//! Far-field free-wave spectrum and Kelvin wave-pattern ("wake") evaluation.
//!
//! In thin-ship theory the disturbance far behind a hull is a superposition
//! of free plane waves, one per propagation angle θ ∈ (−π/2, π/2) measured
//! from the track. Deep-water dispersion pins the stationary component at
//! wavenumber k(θ) = ν sec²θ (ν = g/U²), wavevector ν secθ · (1, tanθ), and
//! the complex amplitude density is carried by the same inner integrals
//! Michell's resistance uses:
//!
//! ```text
//! ζ(x, y) = Re ∫_{−π/2}^{π/2} A(θ) exp[i ν secθ (x + y tanθ)] dθ
//! A(θ)    = −(2ν/π) sec³θ · Σ_j conj(F_j(secθ)) exp[−i ν secθ (x_j + y_j tanθ)]
//! ```
//!
//! with F_j = I + iJ per [`crate::michell`] (phases relative to hull j's
//! x-midpoint) and (x_j, y_j) hull j's midpoint in fleet coordinates. The
//! magnitude of A is fixed by the deep-water free-wave resistance identity
//!
//! ```text
//! R_w = ½ π ρ U² ∫_{−π/2}^{π/2} |A(θ)|² cos³θ dθ
//! ```
//!
//! which reproduces Michell's integral term for term; the phase (the −conj)
//! follows Tuck, Scullen & Lazauskas ("Ship-Wave Patterns in the Spirit of
//! Michell", 2001; "Wave Patterns and Minimum Wave Resistance for High-Speed
//! Vessels", 24th Symp. Naval Hydrodynamics, 2002) mapped to this crate's
//! conventions, and reproduces the classical Havelock point-source phase
//! (a source's transverse system trails cos(νx̄ + π/4): bow ⇒ crest first).
//!
//! # Conventions
//!
//! The ship advances toward **+x**: the bow is the high-x end of the hull
//! file and the wake trails toward −x. `y` is transverse (the multihull
//! placement axis) and `ζ` is positive upward. θ > 0 waves propagate with a
//! +y component.
//!
//! # Validity
//!
//! This is the far-field *free-wave* part of the linear solution. It omits
//! the local (non-radiating) disturbance, so it is the full linear answer
//! only aft of the stern; on, abreast of, or ahead of the hull it is not
//! physical. Linear theory also means no breaking: steep bow systems are
//! indicative only.

use crate::conditions::Conditions;
use crate::error::{Error, Result};
use crate::hull::Hull;
use crate::michell::{dipole_weight, InnerIntegral, Placement};
use crate::moments::C64;
use crate::quadrature::gauss_legendre;
use std::f64::consts::{FRAC_PI_2, PI};

/// Free-wave spectrum of a hull or fleet at fixed speed: the complex
/// amplitude density A(θ), the angular resistance density dR_w/dθ, and
/// wave-elevation reconstruction on points, cuts, and grids.
pub struct FreeWaveSpectrum<'h> {
    nu: f64,
    speed: f64,
    rho: f64,
    members: Vec<Member<'h>>,
    /// Longitudinal phase reference: mean hull midpoint.
    x_ref: f64,
    /// Fleet envelope half-extent about `x_ref`, hull lengths included.
    x_half: f64,
    /// Largest |y_j| over the fleet.
    y_span: f64,
    /// Deepest draft in the fleet.
    t_max: f64,
}

struct Member<'h> {
    inner: InnerIntegral<'h>,
    /// Hull midpoint relative to `x_ref`.
    dx: f64,
    /// Transverse centerplane position.
    y: f64,
}

/// One diagonal or pair-interference contribution to `|Σ Aⱼ(θ)|²`.
#[derive(Debug, Clone, Copy)]
pub struct WaveInterferenceContribution {
    /// Left member index in the input fleet.
    pub left: usize,
    /// Right member index (`right >= left`). Equal indices are self terms.
    pub right: usize,
    /// Signed contribution to `|Σ Aⱼ|²` [m²/rad²]. Off-diagonal
    /// interference terms may be negative.
    pub amplitude_squared: f64,
    /// Signed contribution to `dR_w/dθ` [N/rad].
    pub resistance_density: f64,
}

/// Complex per-member wave signature and its pairwise resistance attribution
/// at one propagation angle.
#[derive(Debug, Clone)]
pub struct WaveSignature {
    pub theta: f64,
    /// Complex free-wave amplitudes `Aⱼ(θ)` in input fleet order.
    pub member_amplitudes: Vec<C64>,
    /// `Σ Aⱼ(θ)`.
    pub total_amplitude: C64,
    /// `|Σ Aⱼ(θ)|²`.
    pub total_amplitude_squared: f64,
    /// Upper-triangular self and pair terms, including the factor of two for
    /// off-diagonal terms. These sum to `total_amplitude_squared`.
    pub interference: Vec<WaveInterferenceContribution>,
    /// Total `dR_w/dθ` [N/rad].
    pub total_resistance_density: f64,
}

/// Termination reason for a [`WaveGrid`] spectrum integration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WaveGridOutcome {
    /// A full angular window decayed below the amplitude threshold.
    AmplitudeDecay,
    /// The grid could not resolve shorter waves, which were smoothly tapered.
    ResolutionCap,
    /// The documented λ = 15 spectral cap was reached.
    SpectralCap,
    /// The integration exhausted its θ-node budget.
    EvalCap,
}

impl WaveGridOutcome {
    /// Stable lower-case name for CLI diagnostics.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::AmplitudeDecay => "amplitude_decay",
            Self::ResolutionCap => "resolution_cap",
            Self::SpectralCap => "spectral_cap",
            Self::EvalCap => "eval_cap",
        }
    }
}

/// Wave elevation ζ sampled on a rectangular grid.
#[derive(Debug, Clone)]
pub struct WaveGrid {
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub nx: usize,
    pub ny: usize,
    /// Row-major over y: `zeta[iy * nx + ix]` [m], `iy` from y0 to y1 and
    /// `ix` from x0 to x1 (inclusive endpoints; single row/column when
    /// nx/ny is 1).
    pub zeta: Vec<f64>,
    /// Largest λ = sec θ integrated before truncation.
    pub max_lambda: f64,
    /// Number of θ nodes evaluated (per signed pair).
    pub theta_samples: usize,
    /// Why the angular integration stopped.
    pub outcome: WaveGridOutcome,
    /// True when the θ integral was truncated at the grid's resolution
    /// limit (waves shorter than ~2 pixels smoothly tapered away) rather
    /// than by amplitude decay.
    pub resolution_limited: bool,
}

impl WaveGrid {
    /// x coordinate of column `ix`.
    pub fn x(&self, ix: usize) -> f64 {
        grid_coord(self.x0, self.x1, self.nx, ix)
    }

    /// y coordinate of row `iy`.
    pub fn y(&self, iy: usize) -> f64 {
        grid_coord(self.y0, self.y1, self.ny, iy)
    }

    /// ζ at column `ix`, row `iy` [m].
    pub fn get(&self, ix: usize, iy: usize) -> f64 {
        self.zeta[iy * self.nx + ix]
    }
}

fn grid_coord(a: f64, b: f64, n: usize, i: usize) -> f64 {
    if n <= 1 {
        a
    } else {
        a + (b - a) * i as f64 / (n - 1) as f64
    }
}

impl<'h> FreeWaveSpectrum<'h> {
    /// Build the spectrum of a fleet (same member format as
    /// [`crate::multihull_wave_resistance`]).
    pub fn new(members: &[(&'h Hull, Placement)], cond: &Conditions) -> Result<Self> {
        cond.validate()?;
        if members.is_empty() {
            return Err(Error::InvalidConditions(
                "at least one hull is required".into(),
            ));
        }
        if members
            .iter()
            .any(|(_, p)| !(p.x.is_finite() && p.y.is_finite()))
        {
            return Err(Error::InvalidConditions(
                "hull placements must be finite".into(),
            ));
        }
        let nu = cond.gravity / (cond.speed * cond.speed);
        let n = members.len() as f64;
        let x_ref = members.iter().map(|(h, p)| h.x_center() + p.x).sum::<f64>() / n;
        Ok(FreeWaveSpectrum {
            nu,
            speed: cond.speed,
            rho: cond.fluid.density,
            x_half: members
                .iter()
                .map(|(h, p)| (h.x_center() + p.x - x_ref).abs() + h.x_half_extent())
                .fold(0.0, f64::max),
            y_span: members.iter().map(|(_, p)| p.y.abs()).fold(0.0, f64::max),
            t_max: members.iter().map(|(h, _)| h.draft()).fold(0.0, f64::max),
            members: members
                .iter()
                .map(|(h, p)| Member {
                    inner: InnerIntegral::new(h, nu),
                    dx: h.x_center() + p.x - x_ref,
                    y: p.y,
                })
                .collect(),
            x_ref,
        })
    }

    /// Fundamental (transverse-wave) wavenumber ν = g/U² [1/m].
    pub fn wavenumber(&self) -> f64 {
        self.nu
    }

    /// Wavelength of the transverse (θ = 0) waves, 2π U²/g [m].
    pub fn transverse_wavelength(&self) -> f64 {
        2.0 * PI / self.nu
    }

    /// Complex free-wave amplitude density A(θ) [m/rad], θ ∈ (−π/2, π/2).
    ///
    /// Longitudinal phases are referenced to the fleet's mean hull midpoint;
    /// |A| does not depend on that choice. Returns zero outside the physical
    /// range (and where the amplitude underflows).
    pub fn amplitude(&mut self, theta: f64) -> C64 {
        if !(theta.is_finite() && theta.abs() < FRAC_PI_2) {
            return C64::ZERO;
        }
        let sec = 1.0 / theta.cos();
        if sec > 1e8 {
            return C64::ZERO;
        }
        let (plus, minus) = self.amp_pair(sec, self.nu * sec * theta.tan().abs());
        if theta >= 0.0 {
            plus
        } else {
            minus
        }
    }

    /// dR_w/dθ [N/rad]: the angular density of wave resistance,
    /// ½ π ρ U² |A(θ)|² cos³θ. Integrating over (−π/2, π/2) recovers the
    /// Michell wave resistance.
    pub fn resistance_density(&mut self, theta: f64) -> f64 {
        let a = self.amplitude(theta);
        let c = theta.cos();
        0.5 * PI * self.rho * self.speed * self.speed * a.abs_sq() * c * c * c
    }

    /// Per-member complex amplitudes and pairwise interference attribution at
    /// one signed propagation angle `θ`.
    ///
    /// Member amplitudes are in the same fleet order passed to [`Self::new`]
    /// and sum exactly to [`WaveSignature::total_amplitude`]. The upper-
    /// triangular interference array expands `|Σ Aⱼ|²`: diagonal entries
    /// are `|Aⱼ|²`, while off-diagonal entries are
    /// `2 Re(Aⱼ conj(Aⱼ))` and can be negative at favourable-interference
    /// angles. Outside `|θ| < π/2`, every amplitude and contribution is zero.
    pub fn signature(&mut self, theta: f64) -> WaveSignature {
        let member_amplitudes = self.member_amplitudes(theta);
        let total_amplitude = member_amplitudes
            .iter()
            .copied()
            .fold(C64::ZERO, |sum, amplitude| sum + amplitude);
        let density_scale = if theta.is_finite() {
            let cosine = theta.cos();
            0.5 * PI * self.rho * self.speed * self.speed * cosine * cosine * cosine
        } else {
            0.0
        };
        let mut interference =
            Vec::with_capacity(member_amplitudes.len() * (member_amplitudes.len() + 1) / 2);
        for left in 0..member_amplitudes.len() {
            for right in left..member_amplitudes.len() {
                let product = member_amplitudes[left] * conjugate(member_amplitudes[right]);
                let amplitude_squared = if left == right {
                    product.re
                } else {
                    2.0 * product.re
                };
                interference.push(WaveInterferenceContribution {
                    left,
                    right,
                    amplitude_squared,
                    resistance_density: density_scale * amplitude_squared,
                });
            }
        }
        let total_amplitude_squared = total_amplitude.abs_sq();
        WaveSignature {
            theta,
            member_amplitudes,
            total_amplitude,
            total_amplitude_squared,
            interference,
            total_resistance_density: density_scale * total_amplitude_squared,
        }
    }

    /// ζ at a single point [m]. See [`Self::elevation_grid`].
    pub fn elevation_at(&mut self, x: f64, y: f64) -> Result<f64> {
        Ok(self.elevation_grid(x, x, y, y, 1, 1)?.zeta[0])
    }

    /// Wave elevation on an inclusive rectangular grid (`nx` × `ny` points
    /// spanning [x0, x1] × [y0, y1]; a cut is `ny = 1`).
    ///
    /// The θ integral marches oscillation-rate-sized Gauss–Legendre panels
    /// and truncates when the amplitude has decayed, or — on grids too
    /// coarse for the shortest diverging waves — at the grid's resolution
    /// limit, tapering smoothly to avoid aliasing artefacts (reported via
    /// [`WaveGrid::resolution_limited`]).
    ///
    /// For wall-sided hulls |A(θ)| decays only algebraically toward θ =
    /// ±π/2, so the diverging system carries ever-shorter ripple near the
    /// track; components beyond λ = secθ = 15 (wavelength below the
    /// transverse wave's /225) are always tapered away — in reality
    /// viscosity destroys them (cf. the eddy-viscosity factor of Tuck,
    /// Scullen & Lazauskas 2002).
    pub fn elevation_grid(
        &mut self,
        x0: f64,
        x1: f64,
        y0: f64,
        y1: f64,
        nx: usize,
        ny: usize,
    ) -> Result<WaveGrid> {
        const MAX_NODES: usize = 400_000;
        /// Amplitude-truncation: a full quiet window relative to the peak.
        const QUIET_REL: f64 = 1e-4;
        const QUIET_WINDOW_PHASE: f64 = 6.0 * PI;

        for v in [x0, x1, y0, y1] {
            if !v.is_finite() {
                return Err(Error::InvalidConditions(
                    "grid extents must be finite".into(),
                ));
            }
        }
        if nx == 0 || ny == 0 {
            return Err(Error::InvalidConditions(
                "grid must have at least one point per axis".into(),
            ));
        }
        if (nx > 1 && x1 <= x0) || (ny > 1 && y1 <= y0) {
            return Err(Error::InvalidConditions(
                "grid extents must be increasing".into(),
            ));
        }
        let hx = if nx > 1 {
            (x1 - x0) / (nx - 1) as f64
        } else {
            0.0
        };
        let hy = if ny > 1 {
            (y1 - y0) / (ny - 1) as f64
        } else {
            0.0
        };
        let nu = self.nu;

        // Resolution cutoff: the largest θ whose wave the grid can render
        // (kx·hx and ky·hy both under π), found by bisection; the last 30%
        // (in sec²θ) before the cutoff is cosine-tapered.
        let resolvable = |theta: f64| -> bool {
            let sec = 1.0 / theta.cos();
            nu * sec * hx < PI && nu * sec * theta.tan() * hy < PI
        };
        /// Spectral cap: components with λ = secθ beyond this are dropped.
        const MAX_SEC: f64 = 15.0;
        let theta_hi = FRAC_PI_2 - 1e-9;
        if !resolvable(0.0) {
            return Err(Error::InvalidConditions(format!(
                "grid spacing {hx:.3} m cannot resolve the transverse wavelength \
                 {:.3} m; use at least ~4 points per wavelength",
                2.0 * PI / nu
            )));
        }
        let theta_pix = if resolvable(theta_hi) {
            theta_hi
        } else {
            let (mut lo, mut hi) = (0.0f64, theta_hi);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if resolvable(mid) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        };
        let theta_end = theta_pix.min((1.0 / MAX_SEC).acos());
        let grid_capped = theta_end < (1.0 / MAX_SEC).acos();
        let sec_end = 1.0 / theta_end.cos();
        let taper_s0 = 1.0f64.max(0.7 * sec_end * sec_end);
        let taper = |sec: f64| -> f64 {
            let s = sec * sec;
            if s <= taper_s0 {
                1.0
            } else {
                0.5 * (1.0 + (PI * (s - taper_s0) / (sec_end * sec_end - taper_s0)).cos())
            }
        };

        // Phase-rate bound over the grid for panel sizing.
        let x_tot = self.x_half + (x0 - self.x_ref).abs().max((x1 - self.x_ref).abs());
        let y_tot = self.y_span + y0.abs().max(y1.abs());
        let t_max = self.t_max;
        let rate = |theta: f64| -> f64 {
            let sec = 1.0 / theta.cos();
            let tan = theta.tan();
            nu * sec * tan * (x_tot + 2.0 * t_max * sec)
                + nu * y_tot * sec * (sec * sec + tan * tan)
                + 2.0
        };
        // Near θ = 0 the longitudinal phase is quadratic in θ; cap the first
        // panels to one period of it (cf. the resistance integrator).
        let cap = (2.0 * PI / (2.0 * nu * x_tot).sqrt().max(1.0)).min(0.12);
        const FRAC: f64 = 0.6;

        let (gx, gw) = gauss_legendre(16);
        let mut zeta = vec![0.0f64; nx * ny];
        let mut p_re = vec![0.0f64; nx];
        let mut p_im = vec![0.0f64; nx];

        let mut theta = 0.0f64;
        let mut nodes = 0usize;
        let mut amp_peak = 0.0f64;
        let mut window_phase = 0.0f64;
        let mut window_peak = 0.0f64;
        let mut resolution_limited = grid_capped;
        let mut outcome = if grid_capped {
            WaveGridOutcome::ResolutionCap
        } else {
            WaveGridOutcome::SpectralCap
        };
        while theta < theta_end - 1e-12 {
            let local_rate = rate(theta);
            let dt = (FRAC * 2.0 * PI / local_rate)
                .min(FRAC * cap)
                .min(theta_end - theta)
                .max(1e-15);
            let half = dt / 2.0;
            let mid = theta + half;
            let mut panel_peak = 0.0f64;
            for (i, &xi) in gx.iter().enumerate() {
                let th = mid + half * xi;
                let sec = 1.0 / th.cos();
                let ky = nu * sec * th.tan();
                let (mut ap, mut am) = self.amp_pair(sec, ky);
                let mag = ap.abs().max(am.abs());
                panel_peak = panel_peak.max(mag);
                if mag == 0.0 {
                    continue;
                }
                let w = gw[i] * half * taper(sec);
                ap = ap.scale(w);
                am = am.scale(w);

                // Column phasors e^{i kx (x - x_ref)} by recurrence.
                let kx = nu * sec;
                let step = C64::cis(kx * hx);
                let mut cur = C64::cis(kx * (x0 - self.x_ref));
                for ix in 0..nx {
                    p_re[ix] = cur.re;
                    p_im[ix] = cur.im;
                    cur = cur * step;
                }
                for iy in 0..ny {
                    let y = grid_coord(y0, y1, ny, iy);
                    let ph = C64::cis(ky * y);
                    // c = A(+θ) e^{i ky y} + A(−θ) e^{−i ky y}
                    let c = ap * ph + am * C64::new(ph.re, -ph.im);
                    let row = &mut zeta[iy * nx..(iy + 1) * nx];
                    for ix in 0..nx {
                        row[ix] += c.re * p_re[ix] - c.im * p_im[ix];
                    }
                }
            }
            nodes += 16;
            theta += dt;
            amp_peak = amp_peak.max(panel_peak);

            // Amplitude truncation past λ = 2: stop once a full window of
            // accumulated phase stayed far below the spectrum's peak.
            if 1.0 / theta.cos() > 2.0 {
                window_phase += local_rate * dt;
                window_peak = window_peak.max(panel_peak);
                if window_phase >= QUIET_WINDOW_PHASE {
                    if window_peak <= QUIET_REL * amp_peak {
                        resolution_limited = false;
                        outcome = WaveGridOutcome::AmplitudeDecay;
                        break;
                    }
                    window_phase = 0.0;
                    window_peak = 0.0;
                }
            }
            if nodes >= MAX_NODES {
                outcome = WaveGridOutcome::EvalCap;
                break;
            }
        }

        Ok(WaveGrid {
            x0,
            x1,
            y0,
            y1,
            nx,
            ny,
            zeta,
            max_lambda: 1.0 / theta.cos().max(1e-300),
            theta_samples: nodes,
            outcome,
            resolution_limited,
        })
    }

    /// A(+θ) and A(−θ) for one |θ|, sharing the inner-integral evaluations
    /// (F depends only on sec θ; only the transverse placement phase flips).
    fn amp_pair(&mut self, sec: f64, ky_abs: f64) -> (C64, C64) {
        let kx = self.nu * sec;
        let mut plus = C64::ZERO;
        let mut minus = C64::ZERO;
        for m in self.members.iter_mut() {
            let (source, camber) = m.inner.eval_pair(sec);
            let weighted_camber = camber.map_or(C64::ZERO, |value| value.scale(dipole_weight(sec)));
            let system_plus = source - weighted_camber;
            let system_minus = source + weighted_camber;
            if system_plus == C64::ZERO && system_minus == C64::ZERO {
                continue;
            }
            plus = plus + conjugate(system_plus) * C64::cis(-(kx * m.dx + ky_abs * m.y));
            minus = minus + conjugate(system_minus) * C64::cis(-(kx * m.dx - ky_abs * m.y));
        }
        if plus == C64::ZERO && minus == C64::ZERO {
            return (C64::ZERO, C64::ZERO);
        }
        let scale = -(2.0 * self.nu / PI) * sec * sec * sec;
        (plus.scale(scale), minus.scale(scale))
    }

    fn member_amplitudes(&mut self, theta: f64) -> Vec<C64> {
        if !(theta.is_finite() && theta.abs() < FRAC_PI_2) {
            return vec![C64::ZERO; self.members.len()];
        }
        let sec = 1.0 / theta.cos();
        if sec > 1e8 {
            return vec![C64::ZERO; self.members.len()];
        }
        let kx = self.nu * sec;
        let ky = self.nu * sec * theta.tan();
        let scale = -(2.0 * self.nu / PI) * sec * sec * sec;
        self.members
            .iter_mut()
            .map(|member| {
                let (source, camber) = member.inner.eval_pair(sec);
                let weighted_camber =
                    camber.map_or(C64::ZERO, |value| value.scale(dipole_weight(sec)));
                let system = if theta >= 0.0 {
                    source - weighted_camber
                } else {
                    source + weighted_camber
                };
                (conjugate(system) * C64::cis(-(kx * member.dx + ky * member.y))).scale(scale)
            })
            .collect()
    }
}

#[inline]
fn conjugate(value: C64) -> C64 {
    C64::new(value.re, -value.im)
}

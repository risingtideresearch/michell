//! # michell-seakeeping
//!
//! Linear **strip-theory seakeeping** for sectional hulls: heave and pitch in
//! regular waves, in the frequency domain, following Salvesen, Tuck &
//! Faltinsen (1970). It is built on [`michell_geometry`] alone — the
//! thin-ship wave kernel (`michell`) plays no part in it; the two share the
//! hull geometry and the equilibrium attitude it is linearised about.
//!
//! ## Conventions
//!
//! - Hull frame as in [`michell_geometry`]: `x` forward (the bow is the
//!   high-`x` end), sections measure `z` **down** from the waterline.
//! - Heave `η₃` is positive **up**; pitch `η₅` is positive **bow up**, about a
//!   pivot `x_ref` (usually the longitudinal centre of gravity).
//! - Waves travel at heading `β` to the ship's course: `β = 0` following
//!   seas, `β = π` head seas. The incident elevation, in the frame moving
//!   with the ship at speed `U`, is `ζ = Re[ζ_a e^{i k (x − x_ref) cos β}
//!   e^{−i ω_e t}]` with the encounter frequency `ω_e = ω − k U cos β`
//!   ([`waves`]): complex amplitudes here are per unit wave amplitude, with
//!   their phase relative to a crest at `x_ref`.
//!
//! ## What is here
//!
//! - [`restoring`]: hydrostatic restoring coefficients `C₃₃, C₃₅, C₅₅`
//!   from the waterplane;
//! - [`froude_krylov`]: the incident-wave heave force and pitch moment, in
//!   closed form per x-span from the sectional hull's depth integrals
//!   ([`michell_geometry::SectionalHull::x_transform`]);
//! - [`green`], [`section2d`]: the 2-D pulsating source and Frank's
//!   close-fit method on each station's section curve — added mass,
//!   damping, radiated waves, and the diffraction force by Haskind;
//! - [`strip`]: Salvesen–Tuck–Faltinsen assembly at forward speed (transom
//!   terms included) and the heave–pitch response in regular waves;
//! - [`strip::added_resistance`]: mean added resistance by Gerritsma &
//!   Beukelman's radiated energy (no short-wave correction);
//! - [`sea`]: Bretschneider and JONSWAP spectra, significant motions and
//!   accelerations, mean added resistance in a sea.
//!
//! Still to come: a short-wave added-resistance correction; oblique-sea Froude–Krylov with
//! its transverse variation.
//!
//! Checks: the source's principal value against quadrature, the section
//! damping against the energy its waves carry, the exact infinite-frequency
//! added mass of a semicircle (and the approach to it), the Haskind diffraction
//! force against the solved diffraction problem, long-wave limits of the
//! full response. The external target is Journée's (1992) Wigley data —
//! heave and pitch RAOs and added resistance at Fn 0.2–0.4 — not yet
//! compared (the report was not reachable while this was written).

pub mod froude_krylov;
pub mod green;
pub mod linalg;
pub mod restoring;
pub mod sea;
pub mod section2d;
pub mod strip;
pub mod waves;

pub use froude_krylov::{froude_krylov, WaveLoad};
pub use restoring::{restoring, Restoring};

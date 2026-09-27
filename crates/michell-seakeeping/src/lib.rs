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
//! ## What is here, and what comes next
//!
//! Done — the parts that need only the geometry:
//!
//! - [`restoring`]: hydrostatic restoring coefficients `C₃₃, C₃₅, C₅₅`
//!   from the waterplane;
//! - [`froude_krylov`]: the incident-wave (Froude–Krylov) heave force and
//!   pitch moment, in closed form per x-span from the sectional hull's
//!   depth integrals `Z(x; k)` ([`michell_geometry::SectionalHull::x_transform`]).
//!
//! To come:
//!
//! 1. 2-D section hydrodynamics: added mass `a₃₃(x, ω_e)`, damping
//!    `b₃₃(x, ω_e)` and the diffraction force per station — Lewis forms
//!    first, then a Frank close-fit source method on the section curves
//!    ([`michell_geometry::SectionalHull::curves`]);
//! 2. assembly into the forward-speed coefficients `A_jk, B_jk`
//!    (Salvesen–Tuck–Faltinsen), the mass matrix (needs a pitch radius of
//!    gyration on the load case) and the 2×2 heave–pitch RAOs;
//! 3. irregular seas: spectra, encounter-frequency mapping, significant
//!    motions and accelerations, and added resistance (Gerritsma–Beukelman
//!    with a short-wave correction).
//!
//! Validation targets: the Wigley hulls I–III of Journée (1992) — heave and
//! pitch RAOs and added resistance at Fn 0.2, 0.3 and 0.4.

pub mod froude_krylov;
pub mod restoring;
pub mod waves;

pub use froude_krylov::{froude_krylov, WaveLoad};
pub use restoring::{restoring, Restoring};

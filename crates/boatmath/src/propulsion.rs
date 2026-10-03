//! A case's hulls at an attitude, ready to say what wake a propeller sits
//! in and what thrust deduction it costs (see `thinship::propulsion`).

use crate::params::CaseParams;
use crate::{setup, LoftRequest};
use hullgeom::{Conditions, Placement, SectionalHull};
use thinship::nearfield::NearFieldOptions;
use thinship::propulsion::{Disc, Interaction};
use thinship::TransomClosure;

/// The hulls at an attitude and speed, their singularities placed.
pub struct Placed {
    pub interaction: Interaction,
    pub mounts: Mounts,
}

/// The platform's attitude, for carrying a propeller fixed to a hull to
/// where the attitude takes it.
pub struct Mounts {
    sinkage: f64,
    trim: f64,
    pivot: f64,
}

fn mounts_of(attitude: (f64, f64), pivot: f64) -> Mounts {
    Mounts {
        sinkage: attitude.0,
        trim: attitude.1,
        pivot,
    }
}

/// The case at an attitude, for placing propellers (no cut, no flow).
pub fn mounts(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    c: &CaseParams,
    attitude: (f64, f64),
) -> Result<Mounts, String> {
    let s = setup(name, bytes, cut, c)?;
    Ok(mounts_of(
        attitude,
        s.platform(attitude.0, attitude.1).pivot_x,
    ))
}

pub fn placed(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    c: &CaseParams,
    attitude: (f64, f64),
    speed: f64,
    closure: TransomClosure,
) -> Result<Placed, String> {
    let s = setup(name, bytes, cut, c)?;
    let platform = s.platform(attitude.0, attitude.1);
    let owned = s.situate(&platform)?;
    let members: Vec<(&SectionalHull, Placement)> = owned.iter().map(|(h, p)| (h, *p)).collect();
    let cond = Conditions::seawater(speed);
    let opts = NearFieldOptions {
        closure,
        ..NearFieldOptions::default()
    };
    let interaction = Interaction::new(&members, &cond, &opts).map_err(|e| e.to_string())?;
    Ok(Placed {
        interaction,
        mounts: mounts_of(attitude, platform.pivot_x),
    })
}

impl Mounts {
    /// A disc of `radius` whose centre is at `(x, y, z)` in the hull's design
    /// frame (z up from the design waterline), moved with the hull to the
    /// attitude.
    pub fn disc_at(&self, x: f64, y: f64, z: f64, radius: f64) -> Disc {
        Disc {
            x,
            y,
            depth: -z * self.trim.cos() + self.sinkage - (x - self.pivot) * self.trim.sin(),
            radius,
            hub: 0.2,
        }
    }

    /// The discs of a case's mount, one per drive, each hub at least the
    /// pod it sits on (with a margin, so the blades clear it).
    pub fn mount_discs(&self, thrusts: &[crate::mount::Thrust], radius: f64) -> Vec<Disc> {
        thrusts
            .iter()
            .map(|t| Disc {
                hub: (1.05 * t.hub / radius).max(0.2),
                ..self.disc_at(t.x, t.y, t.z, radius)
            })
            .collect()
    }
}

//! A case's hulls at an attitude, ready to say what wake a propeller sits
//! in and what thrust deduction it costs (see `michell::propulsion`).

use crate::params::CaseParams;
use crate::{setup, LoftRequest};
use michell::nearfield::NearFieldOptions;
use michell::propulsion::{Disc, Interaction};
use michell::TransomClosure;
use michell_geometry::{Conditions, Placement, SectionalHull};

/// Where each hull's propellers sit, in its own design frame: forward of
/// its aft end, out from its centreplane (two propellers per hull go
/// either side), and the shaft's depth below the design waterline.
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize, PartialEq)]
pub struct Position {
    pub forward_of_aft: f64,
    pub outboard: f64,
    pub depth: f64,
}

/// The hulls at an attitude and speed, their singularities placed.
pub struct Placed {
    pub interaction: Interaction,
    /// Each hull's aft end and centreplane at the attitude [m, fleet x, y].
    pub hulls: Vec<(f64, f64)>,
    sinkage: f64,
    trim: f64,
    pivot: f64,
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
        hulls: members
            .iter()
            .map(|(h, p)| (h.x_range().0 + p.x, p.y))
            .collect(),
        sinkage: attitude.0,
        trim: attitude.1,
        pivot: platform.pivot_x,
    })
}

impl Placed {
    /// The discs of `per_hull` propellers of `radius` on each hull, moved
    /// with it: its sinkage takes them down, and a bow-up trim takes the
    /// ones aft of the pivot deeper.
    pub fn discs(&self, at: &Position, per_hull: usize, radius: f64) -> Vec<Disc> {
        let sides: &[f64] = if per_hull >= 2 { &[1.0, -1.0] } else { &[1.0] };
        let mut out = Vec::new();
        for &(x_aft, y) in &self.hulls {
            let x = x_aft + at.forward_of_aft;
            let depth =
                at.depth * self.trim.cos() + self.sinkage - (x - self.pivot) * self.trim.sin();
            for &side in sides {
                out.push(Disc {
                    x,
                    y: y + side * at.outboard,
                    depth,
                    radius,
                    hub: 0.2,
                });
            }
        }
        out
    }
}

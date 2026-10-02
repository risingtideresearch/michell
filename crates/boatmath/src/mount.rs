//! A propulsion drive on a hull: its leg and pod as geometry, and where its
//! propeller sits. The parts are members of the platform alongside the hull
//! — appended to the hull's JSON geometry as hulls of their own, each with
//! its `role` and viscous form factor — so the solver cuts, floats and
//! integrates them as it does the hull: their wave resistance, their near
//! field and their part in the attitude are the thin-ship model's.
//!
//! Frame: the hull's design frame, x forward, y to port, z up with the
//! design waterline at z = 0, metres.

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    /// A leg down through the hull bottom to a pod (or a gear housing).
    Saildrive,
    /// A leg down from a bracket on the transom, piercing the water.
    Outboard,
    /// A pod hung close under the hull on a short leg (strut).
    Pod,
}

/// A vertical leg or strut: a symmetric NACA 4-digit section.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Leg {
    pub chord: f64,
    pub thickness: f64,
}

/// A pod: a body of revolution on the shaft line.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pod {
    pub length: f64,
    pub diameter: f64,
    /// How far its nose is ahead of the leg's mid-chord.
    pub nose_ahead: f64,
}

/// A drive on each hull of a case.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Mount {
    pub kind: Kind,
    /// The stock drive it's taken from, if one.
    #[serde(default)]
    pub stock: Option<String>,
    /// The leg's mid-chord forward of the hull's aft end (an outboard's is
    /// astern of it: negative) [m].
    pub x: f64,
    /// Out from the hull's centreplane [m].
    #[serde(default)]
    pub y: f64,
    /// The shaft centreline below the hull's bottom at the leg (saildrive,
    /// pod), or below the transom's bottom (outboard) [m].
    pub shaft_depth: f64,
    pub leg: Leg,
    #[serde(default)]
    pub pod: Option<Pod>,
    /// The propeller ahead of the pod rather than astern of it.
    #[serde(default)]
    pub tractor: bool,
    /// The propeller's diameter, if the drive comes with one [m].
    #[serde(default)]
    pub prop_diameter: Option<f64>,
    /// The shaft's angle to the hull's baseline, bow up [deg].
    #[serde(default)]
    pub shaft_angle_deg: f64,
}

/// Where the propeller sits and points, in the hull's design frame.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Thrust {
    /// The propeller plane's centre.
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// The shaft's angle to the baseline, bow up [rad].
    pub angle: f64,
}

/// Hoerner's form factor `k` (as in `(1 + k) C_F`) of a foil of thickness
/// ratio `t/c`.
pub fn foil_form_factor(t_over_c: f64) -> f64 {
    2.0 * t_over_c + 60.0 * t_over_c.powi(4)
}

/// Hoerner's form factor of a body of revolution of diameter `d` and
/// length `l`.
pub fn body_form_factor(d: f64, l: f64) -> f64 {
    let r = d / l;
    1.5 * r.powf(1.5) + 7.0 * r.powi(3)
}

/// The NACA 4-digit half-thickness at `s` along the chord, closed trailing
/// edge, for thickness ratio `t`.
fn naca(s: f64, t: f64) -> f64 {
    5.0 * t
        * (0.2969 * s.sqrt() - 0.1260 * s - 0.3516 * s * s + 0.2843 * s.powi(3)
            - 0.1036 * s.powi(4))
}

/// A hull's keel depth below the design waterline at `x`, from its cut
/// (`boatmath::loft` of the bare hull): the station nearest `x`.
pub fn keel_depth_at(cut: &Value, hull: usize, x: f64) -> Option<f64> {
    let h = &cut["hulls"][hull];
    let px = h["placement"]["x"].as_f64().unwrap_or(0.0);
    let keel: Vec<(f64, f64)> = h["keel"]
        .as_array()?
        .iter()
        .filter_map(|p| Some((p.get(0)?.as_f64()? + px, p.get(1)?.as_f64()?)))
        .collect();
    // Interpolated between the stations either side.
    let i = keel.iter().position(|k| k.0 >= x)?;
    if i == 0 {
        return Some(keel[0].1);
    }
    let ((xa, da), (xb, db)) = (keel[i - 1], keel[i]);
    let t = if xb > xa { (x - xa) / (xb - xa) } else { 0.0 };
    Some(da + t * (db - da))
}

/// The aft end of a hull and its centreplane, from its cut.
fn aft_end(cut: &Value, hull: usize) -> Option<(f64, f64)> {
    let h = &cut["hulls"][hull];
    let px = h["placement"]["x"].as_f64().unwrap_or(0.0);
    let x = h["keel"]
        .as_array()?
        .iter()
        .filter_map(|k| k.get(0)?.as_f64())
        .fold(f64::INFINITY, f64::min);
    x.is_finite()
        .then(|| (x + px, h["placement"]["y"].as_f64().unwrap_or(0.0)))
}

/// The leg as two patches (starboard and port), its section swept from
/// `z_top` down to `z_bottom`, mid-chord at `x_mid`, on the plane `y`.
fn leg_patches(
    leg: &Leg,
    x_mid: f64,
    y: f64,
    z_top: f64,
    z_bottom: f64,
) -> Result<Vec<Value>, String> {
    let (nu, nv) = (33, 13);
    let t = leg.thickness / leg.chord;
    let x_le = x_mid + 0.5 * leg.chord;
    // The foot closes: over the last fifth of the span the section thins to
    // nothing, rounded, so a station through the leg is a closed section
    // (its foot sits inside the pod, if there is one).
    let foot = |v: f64| -> f64 {
        if v < 0.8 {
            1.0
        } else {
            (1.0 - ((v - 0.8) / 0.2).powi(2)).max(0.0).sqrt()
        }
    };
    let mut out = Vec::new();
    for side in [1.0, -1.0] {
        let mut pts = Vec::with_capacity(nu * nv);
        for iu in 0..nu {
            // Cosine spacing: fine at the nose and the tail.
            let s = 0.5 * (1.0 - (std::f64::consts::PI * iu as f64 / (nu - 1) as f64).cos());
            for iv in 0..nv {
                let v = iv as f64 / (nv - 1) as f64;
                let z = z_top + (z_bottom - z_top) * v;
                // Along the chord from the leading edge (forward) aft.
                pts.push([
                    x_le - s * leg.chord,
                    y + side * leg.chord * naca(s, t) * foot(v),
                    z,
                ]);
            }
        }
        let s =
            michell_geometry::iges::interpolate_grid(nu, nv, &pts).map_err(|e| e.to_string())?;
        out.push(crate::native::patch_value(&s));
    }
    Ok(out)
}

/// The pod as two patches, a body of revolution from its nose at `x_nose`
/// aft, its axis at `(y, z)`: a half-cosine nose and tail on a parallel
/// middle.
fn pod_patches(pod: &Pod, x_nose: f64, y: f64, z: f64) -> Result<Vec<Value>, String> {
    let (nu, nv) = (25, 9);
    let r_max = 0.5 * pod.diameter;
    // Radius along the length: a fair nose over the first 30 %, parallel,
    // and a longer run aft over the last 45 %.
    let radius = |s: f64| -> f64 {
        if s < 0.3 {
            r_max * (1.0 - (1.0 - s / 0.3).powi(2)).sqrt()
        } else if s < 0.55 {
            r_max
        } else {
            let u = (s - 0.55) / 0.45;
            r_max * (1.0 - u * u).max(0.0).sqrt()
        }
    };
    let mut out = Vec::new();
    for side in [1.0, -1.0] {
        let mut pts = Vec::with_capacity(nu * nv);
        for iu in 0..nu {
            let s = 0.5 * (1.0 - (std::f64::consts::PI * iu as f64 / (nu - 1) as f64).cos());
            let r = radius(s);
            for iv in 0..nv {
                // Round from the top (θ = 0) down the side to the bottom.
                let th = std::f64::consts::PI * iv as f64 / (nv - 1) as f64;
                pts.push([
                    x_nose - s * pod.length,
                    y + side * r * th.sin(),
                    z + r * th.cos(),
                ]);
            }
        }
        let s =
            michell_geometry::iges::interpolate_grid(nu, nv, &pts).map_err(|e| e.to_string())?;
        out.push(crate::native::patch_value(&s));
    }
    Ok(out)
}

/// A hull's geometry with `m` on each of its hulls: the parts appended as
/// members (`role` leg or pod, with their form factors), and where each
/// propeller sits. `cut` is the bare hull's (`boatmath::loft`).
pub fn mounted(geometry: &Value, cut: &Value, m: &Mount) -> Result<(Value, Vec<Thrust>), String> {
    if geometry["kind"] != "nurbs" {
        return Err("a mount needs a hull of B-spline patches (from IGES), not a mesh".into());
    }
    let n_hulls = geometry["hulls"].as_array().map_or(0, Vec::len);
    let mut g = geometry.clone();
    let hulls = g["hulls"].as_array_mut().ok_or("geometry without hulls")?;
    for h in hulls.iter_mut() {
        if h.get("role").is_none() {
            h["role"] = json!("hull");
        }
    }
    let leg_k = foil_form_factor(m.leg.thickness / m.leg.chord);
    let mut thrusts = Vec::new();
    for i in 0..n_hulls {
        let (x_aft, yc) = aft_end(cut, i).ok_or("the hull's cut has no stations")?;
        let x_mid = x_aft + m.x;
        let y = yc + m.y;
        let (z_shaft, z_top) = match m.kind {
            Kind::Outboard => {
                let transom = keel_depth_at(cut, i, x_aft).unwrap_or(0.0);
                // The leg rises clear of the water astern of the transom.
                (-(transom + m.shaft_depth), 0.15 + 0.5 * m.leg.chord)
            }
            Kind::Saildrive | Kind::Pod => {
                let keel = keel_depth_at(cut, i, x_mid).ok_or("no keel depth at the drive")?;
                // A little into the hull, so the leg meets it.
                (-(keel + m.shaft_depth), -keel + 0.02)
            }
        };
        if z_shaft >= z_top {
            return Err(format!(
                "the drive's shaft ({z_shaft:.3} m) is above its leg's top ({z_top:.3} m)"
            ));
        }
        hulls.push(json!({
            "role": "leg",
            "form_factor": leg_k,
            "patches": leg_patches(&m.leg, x_mid, y, z_top, z_shaft)?,
        }));
        let x_le = x_mid + 0.5 * m.leg.chord;
        let prop_x = match m.pod {
            Some(p) => {
                let x_nose = x_mid + p.nose_ahead;
                hulls.push(json!({
                    "role": "pod",
                    "form_factor": body_form_factor(p.diameter, p.length),
                    "patches": pod_patches(&p, x_nose, y, z_shaft)?,
                }));
                if m.tractor {
                    x_nose + 0.05 * p.length
                } else {
                    x_nose - p.length
                }
            }
            // A saildrive's gear housing: the propeller a little aft of the
            // leg's trailing edge.
            None => x_le - 1.6 * m.leg.chord,
        };
        thrusts.push(Thrust {
            x: prop_x,
            y,
            z: z_shaft,
            angle: m.shaft_angle_deg.to_radians(),
        });
    }
    Ok((g, thrusts))
}

/// Each member's role and form factor, from a geometry's JSON (`hull` and
/// none for a bare hull's).
pub fn roles(geometry: &Value) -> Vec<(String, Option<f64>)> {
    geometry["hulls"]
        .as_array()
        .map(|hs| {
            hs.iter()
                .map(|h| {
                    (
                        h["role"].as_str().unwrap_or("hull").to_string(),
                        h["form_factor"].as_f64(),
                    )
                })
                .collect()
        })
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::LoftRequest;

    fn wigley() -> Value {
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell_geometry::iges::write(&surfaces, "wigley").unwrap();
        crate::native::from_file("w.igs", text.into_bytes(), &LoftRequest::default()).unwrap()
    }

    fn cut(g: &Value) -> Value {
        crate::loft(
            "g.json",
            serde_json::to_vec(g).unwrap(),
            &LoftRequest::default(),
        )
        .unwrap()
    }

    fn saildrive() -> Mount {
        Mount {
            kind: Kind::Saildrive,
            stock: None,
            x: 2.5,
            y: 0.0,
            shaft_depth: 0.35,
            leg: Leg {
                chord: 0.25,
                thickness: 0.05,
            },
            pod: Some(Pod {
                length: 0.6,
                diameter: 0.14,
                nose_ahead: 0.2,
            }),
            tractor: false,
            prop_diameter: None,
            shaft_angle_deg: 0.0,
        }
    }

    /// The pod's volume, by its radius profile.
    fn pod_volume(p: &Pod) -> f64 {
        let n = 20000;
        let r_max = 0.5 * p.diameter;
        (0..n)
            .map(|i| {
                let s = (i as f64 + 0.5) / n as f64;
                let r = if s < 0.3 {
                    r_max * (1.0 - (1.0 - s / 0.3).powi(2)).sqrt()
                } else if s < 0.55 {
                    r_max
                } else {
                    let u = (s - 0.55) / 0.45;
                    r_max * (1.0 - u * u).max(0.0).sqrt()
                };
                std::f64::consts::PI * r * r * p.length / n as f64
            })
            .sum()
    }

    /// The parts are cut as members of their own, below the hull, with the
    /// volumes their shapes have: a pod wholly under water, and a leg from
    /// the hull's bottom down to it.
    #[test]
    fn a_saildrive_cuts_as_members() {
        let g = wigley();
        let bare = cut(&g);
        let m = saildrive();
        let (mounted, thrusts) = mounted(&g, &bare, &m).unwrap();
        let roles: Vec<String> = roles(&mounted).into_iter().map(|r| r.0).collect();
        assert_eq!(roles, ["hull", "leg", "pod"]);
        let c = cut(&mounted);
        let hulls = c["hulls"].as_array().unwrap();
        assert_eq!(hulls.len(), 3, "{}", c["notes"]);
        let vol = |i: usize| hulls[i]["displaced_volume"].as_f64().unwrap();
        let draft = |i: usize| hulls[i]["draft"].as_f64().unwrap();
        let (hull_v, leg_v, pod_v) = (vol(0), vol(1), vol(2));
        let exact_hull = 4.0 / 9.0 * 10.0 * 0.625;
        assert!((hull_v - exact_hull).abs() < 1e-6 * exact_hull, "{hull_v}");

        let p = m.pod.unwrap();
        let pv = pod_volume(&p);
        eprintln!(
            "pod: cut {pod_v:.6} m3, exact {pv:.6} m3 ({:+.2}%)",
            100.0 * (pod_v / pv - 1.0)
        );
        assert!((pod_v - pv).abs() < 0.02 * pv);

        // The leg from 0.02 m inside the hull's bottom (0.625 m down) to the
        // shaft 0.35 m below it.
        let span = 0.35 + 0.02;
        let lv = 0.6851
            * m.leg.thickness
            * m.leg.chord
            * span
            * (0.8 + 0.2 * std::f64::consts::FRAC_PI_4);
        eprintln!(
            "leg: cut {leg_v:.6} m3, exact {lv:.6} m3 ({:+.2}%)",
            100.0 * (leg_v / lv - 1.0)
        );
        assert!((leg_v - lv).abs() < 0.03 * lv);
        assert!(
            (draft(2) - (0.625 + 0.35 + 0.07)).abs() < 0.01,
            "pod's draft {}",
            draft(2)
        );

        // The propeller at the pod's tail, on the shaft line.
        let t = thrusts[0];
        assert!((t.z + 0.625 + 0.35).abs() < 1e-9);
        assert!((t.x - (-5.0 + 2.5 + 0.2 - 0.6)).abs() < 1e-6, "{t:?}");
    }
}

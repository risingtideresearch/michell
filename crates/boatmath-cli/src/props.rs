//! Propellers and motors as records: a `prop` is the best B-series
//! propeller for an operating point, with its rpm curve; a `drive` is one
//! motor driving it, at its best reduction and point on that curve; a
//! `motor` is a database entry. See the `propeller` crate for the models.

use crate::records::{expect, study_record};
use crate::store::{id_of, short, Store};
use boatmath::params::StudyParams;
use boatmath::SOLVER_VERSION;
use propeller::bseries::{sweep, Config, CurvePoint, Inputs, SweepOptions};
use propeller::motor::{scatter, Database, Dot, Motor, Options, RankBy, Rating};
use serde_json::{json, Value};

/// A calm-water result's operating point: its speed [m/s] and the thrust
/// [N] its resistance asks for, `R_t / (1 − t)`, and its hulls.
pub fn operating_point(r: &Value, thrust_deduction: f64) -> Result<(f64, f64), String> {
    let id = expect(r, "result")?;
    if r["kind"] != "calm" {
        return Err(format!("result {}: not a calm-water result", short(id)));
    }
    let speed = r["speed"].as_f64().ok_or("a result without a speed")?;
    let rt = r["forces"]["rt"].as_f64().ok_or("a result without R_t")?;
    Ok((speed, rt / (1.0 - thrust_deduction)))
}

/// The number of hulls a result's platform floats on: a shaft each.
pub fn hull_count(store: &Store, r: &Value) -> usize {
    let case = store
        .get("study", r["study"].as_str().unwrap_or(""))
        .ok()
        .flatten()
        .and_then(|s| store.get("case", s["case"].as_str()?).ok().flatten());
    match case {
        Some(c) if !c["params"]["span"].is_null() => 2,
        Some(c) => store
            .get("hull", c["hull"].as_str().unwrap_or(""))
            .ok()
            .flatten()
            .and_then(|h| h["summary"]["hulls"].as_array().map(Vec::len))
            .unwrap_or(1),
        None => 1,
    }
}

/// The result of the same study at Froude number `froude`: the second
/// operating point, which must already be in the store.
pub fn result_at_froude(store: &Store, r: &Value, froude: f64) -> Result<Value, String> {
    let study = store.need("study", r["study"].as_str().unwrap_or(""))?;
    let mut p: StudyParams = serde_json::from_value(study["params"].clone())
        .map_err(|e| format!("study params: {e}"))?;
    p.froude = froude;
    let other = study_record(study["case"].as_str().unwrap_or(""), p.canonical()?);
    let id = other["id"].as_str().unwrap_or("");
    store
        .get("result", id)?
        .filter(|r| r["solver_version"] == SOLVER_VERSION)
        .ok_or_else(|| {
            format!(
                "no result at Fn {froude} for case {}: run that study first \
                 (boatmath study --froude {froude} … | boatmath run)",
                short(study["case"].as_str().unwrap_or(""))
            )
        })
}

/// The best propeller for `inputs`, and its curve, as a record (reused from
/// the store when it's there).
pub fn prop(
    store: &Store,
    inputs: &Inputs,
    result: Option<&str>,
    power_cap: f64,
) -> Result<Value, String> {
    let id = id_of(&json!({ "inputs": inputs, "result": result, "power_cap": power_cap }));
    if let Some(old) = store.get("prop", &id)? {
        if old["solver_version"] == SOLVER_VERSION {
            return Ok(old);
        }
    }
    let c = Config::new(inputs);
    let s = sweep(
        &c,
        &SweepOptions {
            power_cap,
            ..SweepOptions::default()
        },
    );
    let qpc = s.best.as_ref().map(|b| b.eta0 * c.eta_h);
    let r = json!({
        "type": "prop",
        "id": id,
        "result": result,
        "inputs": inputs,
        // Per shaft, as the propeller sees it; and the whole boat's.
        "shaft": { "V_A": c.v_a, "T": c.t_target, "top_T": c.top.map(|t| t.t_target) },
        "boat": { "R_T": c.r_t, "P_E": c.p_e, "eta_H": c.eta_h, "QPC": qpc },
        "feasible": s.feasible,
        "best": s.best,
        "best_unconstrained": s.bestUnconstrained,
        "top_feasible": s.topFeasible,
        "window": { "rpm_lo": s.rpmLo, "rpm_hi": s.rpmHi },
        "feasible_rpm": { "lo": s.feasLo, "hi": s.feasHi },
        "curve": s.points,
        "solver_version": SOLVER_VERSION,
    });
    store.put(&r)?;
    Ok(r)
}

// ------------------------------------------------------------- motors

/// Which motors to consider, as the web page's filters.
#[derive(Default, Clone)]
pub struct Filter {
    pub vendors: Vec<String>,
    pub mapped_only: bool,
    pub single_only: bool,
    pub max_mass: Option<f64>,
    pub max_od: Option<f64>,
    pub max_price: Option<f64>,
}

impl Filter {
    /// A size or price limit is hard: a motor that doesn't publish the
    /// number can't meet it.
    pub fn admits(&self, m: &Motor) -> bool {
        (self.vendors.is_empty()
            || self
                .vendors
                .iter()
                .any(|v| v.eq_ignore_ascii_case(&m.vendor)))
            && (!self.mapped_only || matches!(m.tier.as_deref(), Some("A") | Some("A-")))
            && (!self.single_only || m.stack.is_none_or(|s| s == 0.0))
            && self
                .max_mass
                .is_none_or(|l| m.mass_kg.is_some_and(|v| v <= l))
            && self.max_od.is_none_or(|l| m.od_mm.is_some_and(|v| v <= l))
            && self
                .max_price
                .is_none_or(|l| m.price_usd.is_some_and(|v| v <= l))
    }
}

/// A motor's database entry as a record.
pub fn motor_record(m: &Motor) -> Value {
    let mut r = m.record.clone();
    r["type"] = json!("motor");
    r
}

pub struct MatchOptions {
    pub filter: Filter,
    pub rank_by: RankBy,
    pub all_windings: bool,
    /// A fixed reduction (1: direct drive), else the best for each motor.
    pub ratio: Option<f64>,
    pub ratio_max: f64,
    pub gear_eta: f64,
    pub controller_loss: bool,
    pub peak: bool,
    /// Require the prop's second operating point, through the same gearbox.
    pub check_top: bool,
}

/// Every admitted motor at its best point on `prop`'s curve, ranked, as
/// `drive` records; and the ids of those that can't drive it.
pub fn drives(
    store: &Store,
    db: &Database,
    prop: &Value,
    o: &MatchOptions,
) -> Result<(Vec<Value>, Vec<String>), String> {
    let prop_id = expect(prop, "prop")?;
    let curve: Vec<CurvePoint> = serde_json::from_value(prop["curve"].clone())
        .map_err(|e| format!("prop {}: curve: {e}", short(prop_id)))?;
    let has_top = !prop["inputs"]["top"].is_null();
    let mo = Options {
        rating: if o.peak { Rating::Peak } else { Rating::Cont },
        gear_eta: o.gear_eta,
        controller_eta: if o.controller_loss {
            db.controller_eta
        } else {
            1.0
        },
        ratio: o.ratio,
        ratio_max: o.ratio_max,
        check_top: o.check_top && has_top,
        ..Options::default()
    };
    let sc = scatter(
        db,
        &curve,
        &mo,
        |m| o.filter.admits(m),
        o.rank_by,
        o.all_windings,
    );
    let settings = json!({
        "filter": {
            "vendors": o.filter.vendors, "mapped_only": o.filter.mapped_only,
            "single_only": o.filter.single_only, "max_mass": o.filter.max_mass,
            "max_od": o.filter.max_od, "max_price": o.filter.max_price,
        },
        "ratio": o.ratio, "ratio_max": o.ratio_max, "gear_eta": o.gear_eta,
        "controller_eta": mo.controller_eta, "rating": if o.peak { "peak" } else { "cont" },
        "check_top": mo.check_top,
    });
    let mut out = Vec::new();
    for (rank, d) in sc.dots.iter().enumerate() {
        let r = drive_record(prop_id, d, rank + 1, &settings);
        store.put(&r)?;
        out.push(r);
    }
    Ok((out, sc.unreachable.iter().map(|m| m.id.clone()).collect()))
}

fn drive_record(prop_id: &str, d: &Dot, rank: usize, settings: &Value) -> Value {
    let p = d.point.design.as_ref();
    let run = &d.run;
    json!({
        "type": "drive",
        "id": id_of(&json!({ "prop": prop_id, "motor": d.motor.id, "settings": settings })),
        "prop": prop_id,
        "motor": d.motor.id,
        "rank": rank,
        "settings": settings,
        // The propeller it runs: the curve's point at this shaft speed.
        "rpm": d.point.rpm,
        "P_shaft": p.map(|p| p.P_shaft),
        "propeller": p.map(|p| json!({ "D": p.D, "PD": p.PD, "EAR": p.EAR, "Z": p.Z, "eta0": p.eta0 })),
        // The motor, through its reduction.
        "P_elec": run.P_elec,
        "eta": run.eta,
        "etaWithController": run.etaWithController,
        "etaDrive": run.etaDrive,
        "ratio": run.ratio,
        "gearEta": run.gearEta,
        "motorRpm": run.rpm,
        "motorQ": run.Q,
        "V_bus": run.V_bus,
        "I_dc": run.I_dc,
        "I_arms": run.I_arms,
        "load": d.load,
        "extrapolated": run.extrapolated,
        "top": run.top,
        "tier": d.motor.tier,
    })
}

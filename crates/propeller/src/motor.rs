//! Motors: efficiency and limits at an operating point, the best reduction
//! for a propeller's demand, and each motor's best point along a
//! propeller's curve. Ported from propopt's `web/motorcore.js`, against the
//! motor database vendored in `data/motors.js` (see `data/README.md`).
//!
//! A record's efficiency comes in four kinds:
//!
//! - `map`: the manufacturer's efficiency map, digitised and interpolated,
//!   shifted to this winding's copper loss; off the map, a loss model
//!   fitted to it;
//! - `loss`: `P_loss = c_cu Q² + k_h ω + k_e ω² + k_0`, fitted to the
//!   datasheet;
//! - `pmsm`: the same loss model, with current as the binding limit;
//! - `bldc`: the classic brushless model, `I = Q/K_t + I_0`,
//!   `V = rpm/K_v + I R`, against the bus voltage and current limit.
//!
//! All four meet the same envelope (speed, torque, power, bus voltage, DC
//! current), continuous ratings at the design point and peak at a second.

use crate::bseries::CurvePoint;
use serde::Deserialize;
use serde_json::Value;
use std::collections::HashMap;

/// One stage of reduction.
pub const GEAR_ETA_DEFAULT: f64 = 0.97;
/// The rated point sits exactly on the bounds and published figures are
/// rounded: a hair of slack, as the original.
const SLACK: f64 = 1.0001;

#[derive(Clone, Copy, Debug, Deserialize)]
pub struct Loss {
    pub c_cu: f64,
    pub k_h: f64,
    pub k_e: f64,
    #[serde(default)]
    pub k_0: f64,
}

impl Loss {
    fn at(&self, q: f64, w: f64) -> f64 {
        self.c_cu * q * q + self.k_h * w + self.k_e * w * w + self.k_0
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum Eff {
    Map {
        map: String,
        ccu: f64,
        ccu_map: f64,
        #[serde(default)]
        fallback: Option<Loss>,
    },
    Loss(Loss),
    Pmsm {
        kt: f64,
        r_ll: f64,
        i_cont: f64,
        i_peak: f64,
        #[serde(flatten)]
        loss: Loss,
    },
    Bldc {
        kt: f64,
        kv: f64,
        r: f64,
        i0: f64,
        i_max: f64,
        v_bus: f64,
    },
    #[serde(other)]
    Other,
}

/// The fields of a motor record the model reads. The record itself is
/// kept whole in [`Motor::record`].
#[derive(Clone, Debug, Deserialize)]
pub struct Motor {
    pub id: String,
    pub vendor: String,
    #[serde(default)]
    pub family: Option<String>,
    #[serde(default)]
    pub tier: Option<String>,
    #[serde(default)]
    pub stack: Option<f64>,
    pub rpm_max: Option<f64>,
    pub tq_peak_nm: Option<f64>,
    pub tq_cont_nm: Option<f64>,
    pub p_peak_kw: Option<f64>,
    pub p_cont_kw: Option<f64>,
    #[serde(default)]
    pub v_dc: Option<f64>,
    #[serde(default)]
    pub kt: Option<f64>,
    #[serde(default)]
    pub kv: Option<f64>,
    #[serde(default)]
    pub kv_nom: Option<f64>,
    #[serde(default)]
    pub kv_peak: Option<f64>,
    #[serde(default)]
    pub i_cont: Option<f64>,
    #[serde(default)]
    pub i_peak: Option<f64>,
    #[serde(default)]
    pub i_dc_cont: Option<f64>,
    #[serde(default)]
    pub i_dc_peak: Option<f64>,
    #[serde(default)]
    pub eff_includes_controller: bool,
    #[serde(default)]
    pub mass_kg: Option<f64>,
    #[serde(default)]
    pub od_mm: Option<f64>,
    #[serde(default)]
    pub price_usd: Option<f64>,
    #[serde(default)]
    pub eff: Option<Eff>,
    #[serde(skip)]
    pub record: Value,
}

/// A digitised efficiency map: `rpm` and `tq` as `[lo, hi, n]`, `eta` in
/// units of 1e-4 on the `n_rpm × n_tq` grid, `inside` marking the cells
/// within the published envelope.
#[derive(Clone, Debug, Deserialize)]
pub struct EffMap {
    pub rpm: [f64; 3],
    pub tq: [f64; 3],
    pub eta: Vec<f64>,
    pub inside: Vec<u8>,
}

pub struct Database {
    pub motors: Vec<Motor>,
    pub maps: HashMap<String, EffMap>,
    /// The inverter efficiency put back on every motor-only record.
    pub controller_eta: f64,
    /// The whole database object, as published.
    pub raw: Value,
}

const VENDORED: &str = include_str!("../data/motors.js");

impl Database {
    /// The vendored database.
    pub fn vendored() -> &'static Database {
        static DB: std::sync::OnceLock<Database> = std::sync::OnceLock::new();
        DB.get_or_init(|| Database::from_js(VENDORED).expect("the vendored motor database"))
    }

    /// A database from propopt's `motors.js` (`const MOTOR_DB = {…};`) or
    /// from the bare JSON object.
    pub fn from_js(text: &str) -> Result<Database, String> {
        let start = text.find('{').ok_or("motor database: no object")?;
        let mut de = serde_json::Deserializer::from_str(&text[start..]);
        let raw = Value::deserialize(&mut de).map_err(|e| format!("motor database: {e}"))?;
        Database::from_value(raw)
    }

    pub fn from_value(raw: Value) -> Result<Database, String> {
        let motors = raw["motors"]
            .as_array()
            .ok_or("motor database: no motors")?
            .iter()
            .map(|m| {
                let mut motor: Motor = serde_json::from_value(m.clone())
                    .map_err(|e| format!("motor {}: {e}", m["id"]))?;
                motor.record = m.clone();
                Ok(motor)
            })
            .collect::<Result<Vec<_>, String>>()?;
        let maps = match raw.get("maps") {
            Some(m) => serde_json::from_value(m.clone()).map_err(|e| format!("motor maps: {e}"))?,
            None => HashMap::new(),
        };
        Ok(Database {
            controller_eta: raw["controller_eta"].as_f64().unwrap_or(1.0),
            motors,
            maps,
            raw,
        })
    }

    pub fn motor(&self, id: &str) -> Option<&Motor> {
        self.motors.iter().find(|m| m.id == id)
    }
}

// ------------------------------------------------------------- one point

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rating {
    Cont,
    Peak,
}

/// A second operating point at the same gearbox: the propeller's shaft
/// speed and torque there.
#[derive(Clone, Copy, Debug)]
pub struct Second {
    pub rpm: f64,
    pub q: f64,
}

#[derive(Clone, Debug)]
pub struct Options {
    pub rating: Rating,
    pub gear_eta: f64,
    /// Inverter efficiency, put back on motor-only records.
    pub controller_eta: f64,
    pub second: Option<Second>,
    pub second_rating: Rating,
    pub ratio_min: f64,
    pub ratio_max: f64,
    /// A fixed reduction (1 for direct drive) instead of the best one.
    pub ratio: Option<f64>,
    pub coarse: bool,
    /// Require each point's second condition, through the same gearbox.
    pub check_top: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            rating: Rating::Cont,
            gear_eta: GEAR_ETA_DEFAULT,
            controller_eta: 1.0,
            second: None,
            second_rating: Rating::Peak,
            ratio_min: 1.0,
            ratio_max: 25.0,
            ratio: None,
            coarse: false,
            check_top: false,
        }
    }
}

/// A motor at an operating point, or why it can't be.
#[derive(Clone, Debug)]
pub enum At {
    Ok(Running),
    Limit(String),
}

impl At {
    fn limit(s: &str) -> At {
        At::Limit(s.to_string())
    }

    pub fn ok(self) -> Option<Running> {
        match self {
            At::Ok(r) => Some(r),
            At::Limit(_) => None,
        }
    }
}

#[derive(Clone, Debug, serde::Serialize)]
#[allow(non_snake_case)]
pub struct Running {
    pub eta: f64,
    pub etaWithController: f64,
    pub extrapolated: bool,
    /// At the motor.
    pub rpm: f64,
    pub Q: f64,
    pub P_shaft: f64,
    pub P_elec: f64,
    pub I_arms: Option<f64>,
    pub V_bus: Option<f64>,
    pub I_dc: Option<f64>,
    /// Through a gearbox: its ratio and efficiency, shaft-to-battery, and
    /// the second point.
    pub ratio: f64,
    pub gearEta: f64,
    pub etaDrive: f64,
    pub top: Option<TopRunning>,
}

#[derive(Clone, Debug, serde::Serialize)]
#[allow(non_snake_case)]
pub struct TopRunning {
    pub rpm: f64,
    pub Q: f64,
    pub eta: f64,
    pub P_elec: f64,
    pub V_bus: Option<f64>,
    pub I_dc: Option<f64>,
    pub I_arms: Option<f64>,
    pub shaftRpm: f64,
}

/// The DC bus voltage an operating point needs, or `None` when the record
/// doesn't publish enough to say.
pub fn bus_voltage(m: &Motor, rpm: f64, i_arms: Option<f64>) -> Option<f64> {
    match &m.eff {
        Some(Eff::Pmsm { kt, r_ll, .. }) => {
            let kb_hot = kt / 0.01654;
            let v_ll = kb_hot * rpm / 1000.0 + i_arms.unwrap_or(0.0) * r_ll;
            return Some(std::f64::consts::SQRT_2 * v_ll);
        }
        Some(Eff::Bldc { kv, r, .. }) => return Some(rpm / kv + i_arms.unwrap_or(0.0) * r),
        _ => {}
    }
    let kv0 = m.kv.filter(|&k| k != 0.0)?;
    let mut kv = kv0;
    if let (Some(kv_nom), Some(i_cont), Some(i)) =
        (m.kv_nom, m.i_cont.filter(|&c| c != 0.0), i_arms)
    {
        if i <= i_cont {
            kv = kv0 + (kv_nom - kv0) * (i / i_cont);
        } else if let (Some(kv_peak), Some(i_peak)) = (m.kv_peak, m.i_peak.filter(|&p| p > i_cont))
        {
            let f = ((i - i_cont) / (i_peak - i_cont)).min(1.0);
            kv = kv_nom + (kv_peak - kv_nom) * f;
        } else {
            kv = kv_nom;
        }
    }
    (kv > 0.0).then(|| rpm / kv)
}

fn sample_map(map: &EffMap, rpm: f64, q: f64) -> Option<(f64, bool)> {
    let [r0, r1, nr] = map.rpm;
    let [q0, q1, nq] = map.tq;
    if rpm < r0 || rpm > r1 || q < q0 || q > q1 {
        return None;
    }
    let (nr, nq) = (nr as usize, nq as usize);
    let fr = (rpm - r0) / (r1 - r0) * (nr - 1) as f64;
    let fq = (q - q0) / (q1 - q0) * (nq - 1) as f64;
    let i0 = (fr.floor().max(0.0) as usize).min(nr - 2);
    let j0 = (fq.floor().max(0.0) as usize).min(nq - 2);
    let (tr, tq) = (fr - i0 as f64, fq - j0 as f64);
    let at = |i: usize, j: usize| map.eta[i * nq + j] / 10000.0;
    let inside = |i: usize, j: usize| map.inside[i * nq + j] == 1;
    let e = at(i0, j0) * (1.0 - tr) * (1.0 - tq)
        + at(i0 + 1, j0) * tr * (1.0 - tq)
        + at(i0, j0 + 1) * (1.0 - tr) * tq
        + at(i0 + 1, j0 + 1) * tr * tq;
    let solid =
        inside(i0, j0) && inside(i0 + 1, j0) && inside(i0, j0 + 1) && inside(i0 + 1, j0 + 1);
    Some((e, !solid))
}

/// A motor at a shaft speed and torque at the motor.
pub fn motor_point(db: &Database, m: &Motor, rpm: f64, q: f64, o: &Options) -> Option<At> {
    let rating = o.rating;
    if rpm < 0.0 || q <= 0.0 {
        return None;
    }
    // A bound the record doesn't publish doesn't bind, as in the original.
    let over = |v: f64, lim: Option<f64>| lim.is_some_and(|l| v > l * SLACK);
    if over(rpm, m.rpm_max) {
        return Some(At::limit("speed"));
    }
    let tq_lim = if rating == Rating::Peak {
        m.tq_peak_nm
    } else {
        m.tq_cont_nm
    };
    if over(q, tq_lim) {
        return Some(At::limit("torque"));
    }
    let w = 2.0 * std::f64::consts::PI * rpm / 60.0;
    let p_shaft = q * w;
    let p_lim = if rating == Rating::Peak {
        m.p_peak_kw
    } else {
        m.p_cont_kw
    };
    if over(p_shaft, p_lim.map(|p| 1000.0 * p)) {
        return Some(At::limit("power"));
    }
    let n = m.stack.filter(|&s| s != 0.0).unwrap_or(1.0);
    let qm = q / n;
    let p_one = qm * w;
    let eta;
    // Inside the map (or any other model): taken as measured.
    let extrapolated = false;
    let Some(e) = &m.eff else {
        return Some(At::limit("no-efficiency-data"));
    };
    match e {
        Eff::Map {
            map,
            ccu,
            ccu_map,
            fallback,
        } => {
            let s = db.maps.get(map).and_then(|mp| sample_map(mp, rpm, qm));
            match s {
                Some((se, false)) => {
                    let p_loss_one = p_one * (1.0 / se - 1.0);
                    let p_loss = n * (p_loss_one + (ccu - ccu_map) * qm * qm);
                    eta = p_shaft / (p_shaft + p_loss.max(1e-9));
                }
                _ => {
                    let Some(f) = fallback else {
                        return Some(At::limit("outside-map"));
                    };
                    let p_loss = n * f.at(qm, w);
                    let et = p_shaft / (p_shaft + p_loss.max(1e-9));
                    if !(et > 0.0) || !(et < 1.0) {
                        return Some(At::limit("outside-map"));
                    }
                    return Some(finish(db, m, rpm, q, qm, p_shaft, et, true, o));
                }
            }
        }
        Eff::Bldc {
            kt,
            kv,
            r,
            i0,
            i_max,
            v_bus,
        } => {
            let i = qm / kt + i0;
            if i > i_max * SLACK {
                return Some(At::limit("current"));
            }
            let v_req = rpm / kv + i * r;
            if v_req > v_bus * SLACK {
                return Some(At::limit("voltage"));
            }
            eta = p_one / (v_req * i);
        }
        Eff::Pmsm {
            kt,
            i_cont,
            i_peak,
            loss,
            ..
        } => {
            let i = qm / kt;
            let lim = if rating == Rating::Peak {
                *i_peak
            } else {
                *i_cont
            };
            if i > lim * SLACK {
                return Some(At::limit("current"));
            }
            let p_loss = n * loss.at(qm, w);
            eta = p_shaft / (p_shaft + p_loss.max(1e-9));
        }
        Eff::Loss(loss) => {
            let p_loss = n * loss.at(qm, w);
            eta = p_shaft / (p_shaft + p_loss.max(1e-9));
        }
        Eff::Other => return Some(At::limit("no-efficiency-data")),
    }
    if !(eta > 0.0) || !(eta < 1.0) {
        return Some(At::limit("model"));
    }
    Some(finish(db, m, rpm, q, qm, p_shaft, eta, extrapolated, o))
}

/// The electrical operating point, and the bus-voltage and DC-current limits.
#[allow(clippy::too_many_arguments)]
fn finish(
    _db: &Database,
    m: &Motor,
    rpm: f64,
    q: f64,
    qm: f64,
    p_shaft: f64,
    eta: f64,
    extrapolated: bool,
    o: &Options,
) -> At {
    let i_arms = m.kt.filter(|&k| k != 0.0).map(|k| qm / k);
    let v_bus = bus_voltage(m, rpm, i_arms);
    if let (Some(v), Some(lim)) = (v_bus, m.v_dc.filter(|&v| v != 0.0)) {
        if v > lim * SLACK {
            return At::limit("voltage");
        }
    }
    let ce = if o.controller_eta != 0.0 {
        o.controller_eta
    } else {
        1.0
    };
    let eta_drv = if m.eff_includes_controller {
        eta
    } else {
        eta * ce
    };
    let p_elec = p_shaft / eta_drv;
    let v_for_i = v_bus.filter(|&v| v > 0.0).or(m.v_dc.filter(|&v| v != 0.0));
    let i_dc = v_for_i.map(|v| p_elec / v);
    if let Some(i) = i_dc {
        let lim = if o.rating == Rating::Peak {
            m.i_dc_peak
        } else {
            m.i_dc_cont
        };
        if lim.is_some_and(|l| l != 0.0 && i > l * SLACK) {
            return At::limit("dc-current");
        }
    }
    At::Ok(Running {
        eta,
        etaWithController: eta_drv,
        extrapolated,
        rpm,
        Q: q,
        P_shaft: p_shaft,
        P_elec: p_elec,
        I_arms: i_arms,
        V_bus: v_bus,
        I_dc: i_dc,
        ratio: 1.0,
        gearEta: 1.0,
        etaDrive: eta_drv,
        top: None,
    })
}

// ------------------------------------------------------------- gearing

/// A propeller's demand at `rpm` and `q`, through a reduction of `ratio`.
pub fn through_gear(
    db: &Database,
    m: &Motor,
    rpm: f64,
    q: f64,
    ratio: f64,
    o: &Options,
) -> Option<At> {
    let g = if ratio <= 0.0 { 1.0 } else { ratio };
    let eff_gear = if g == 1.0 { 1.0 } else { o.gear_eta };
    let r = match motor_point(db, m, rpm * g, q / (g * eff_gear), o)? {
        At::Ok(r) => r,
        lim => return Some(lim),
    };
    let mut top = None;
    if let Some(s) = o.second {
        let o2 = Options {
            rating: o.second_rating,
            second: None,
            ..o.clone()
        };
        match motor_point(db, m, s.rpm * g, s.q / (g * eff_gear), &o2) {
            Some(At::Ok(t)) => {
                top = Some(TopRunning {
                    rpm: s.rpm * g,
                    Q: s.q / (g * eff_gear),
                    eta: t.eta,
                    P_elec: t.P_elec,
                    V_bus: t.V_bus,
                    I_dc: t.I_dc,
                    I_arms: t.I_arms,
                    shaftRpm: s.rpm,
                })
            }
            Some(At::Limit(l)) => return Some(At::Limit(format!("top-speed {l}"))),
            None => return Some(At::limit("top-speed unreachable")),
        }
    }
    let p_elec = r.P_elec;
    Some(At::Ok(Running {
        ratio: g,
        gearEta: eff_gear,
        top,
        etaDrive: (q * 2.0 * std::f64::consts::PI * rpm / 60.0) / p_elec,
        ..r
    }))
}

/// The best reduction for one demand: least electrical power.
pub fn best_gear(db: &Database, m: &Motor, rpm: f64, q: f64, o: &Options) -> Option<Running> {
    let (lo, hi) = (o.ratio_min, o.ratio_max);
    let nn = if o.coarse { 26 } else { 80 };
    let mut best: Option<Running> = None;
    for i in 0..nn {
        let g = lo * (hi / lo).powf(i as f64 / (nn - 1) as f64);
        if let Some(r) = through_gear(db, m, rpm, q, g, o).and_then(At::ok) {
            if best.as_ref().is_none_or(|b| r.P_elec < b.P_elec) {
                best = Some(r);
            }
        }
    }
    let mut best = best?;
    if o.coarse {
        return Some(best);
    }
    const PHI: f64 = 0.6180339887;
    let (mut a, mut b) = (lo.max(best.ratio / 1.25), hi.min(best.ratio * 1.25));
    let f = |g: f64| through_gear(db, m, rpm, q, g, o).and_then(At::ok);
    let mut x1 = b - PHI * (b - a);
    let mut x2 = a + PHI * (b - a);
    let mut f1 = f(x1);
    let mut f2 = f(x2);
    let mut i = 0;
    while i < 30 && b - a > 1e-4 {
        let p1 = f1.as_ref().map_or(f64::INFINITY, |r| r.P_elec);
        let p2 = f2.as_ref().map_or(f64::INFINITY, |r| r.P_elec);
        if p1 < p2 {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - PHI * (b - a);
            f1 = f(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + PHI * (b - a);
            f2 = f(x2);
        }
        i += 1;
    }
    for c in [f1, f2].into_iter().flatten() {
        if c.P_elec < best.P_elec {
            best = c;
        }
    }
    Some(best)
}

// ------------------------------------------------------------- the scatter

/// One motor at its best point on a propeller's curve.
#[derive(Clone, Debug)]
pub struct Dot<'a> {
    pub motor: &'a Motor,
    /// The curve's point (the propeller) it runs at.
    pub point: CurvePoint,
    pub run: Running,
    /// How much of the motor's continuous torque the job uses.
    pub load: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RankBy {
    Power,
    Mass,
    Price,
}

pub struct Scatter<'a> {
    /// One per family (or every winding), ranked.
    pub dots: Vec<Dot<'a>>,
    /// Every winding, ranked.
    pub all: Vec<Dot<'a>>,
    pub unreachable: Vec<&'a Motor>,
}

/// Each included motor's best point along `points`: least electrical power,
/// its reduction chosen for that point (or fixed by `o.ratio`).
pub fn scatter<'a>(
    db: &'a Database,
    points: &[CurvePoint],
    o: &Options,
    include: impl Fn(&Motor) -> bool,
    rank_by: RankBy,
    all_windings: bool,
) -> Scatter<'a> {
    let feasible: Vec<&CurvePoint> = points
        .iter()
        .filter(|p| p.ok && (!o.check_top || p.design.as_ref().is_some_and(|d| d.topOk)))
        .collect();
    let mut out = Scatter {
        dots: Vec::new(),
        all: Vec::new(),
        unreachable: Vec::new(),
    };
    if feasible.is_empty() {
        return out;
    }
    let stride = (feasible.len() / 48).max(1);
    let coarse = Options {
        coarse: true,
        ..o.clone()
    };
    let solve = |m: &Motor, p: &CurvePoint, oo: &Options| -> Option<Running> {
        let d = p.design.as_ref()?;
        let mut oo = oo.clone();
        if o.check_top {
            if let Some(t) = &d.top {
                oo.second = Some(Second { rpm: t.rpm, q: t.Q });
            }
        }
        match o.ratio {
            Some(g) => through_gear(db, m, d.rpm, d.Q, g, &oo).and_then(At::ok),
            None => best_gear(db, m, d.rpm, d.Q, &oo),
        }
    };
    let mut dots: Vec<Dot> = Vec::new();
    for m in &db.motors {
        if !include(m) {
            continue;
        }
        let mut bi: Option<usize> = None;
        let mut bp: Option<f64> = None;
        let mut i = 0;
        while i < feasible.len() {
            if let Some(r) = solve(m, feasible[i], &coarse) {
                if bp.is_none_or(|b| r.P_elec < b) {
                    bp = Some(r.P_elec);
                    bi = Some(i);
                }
            }
            i += stride;
        }
        let Some(bi) = bi else {
            out.unreachable.push(m);
            continue;
        };
        let mut best: Option<(Running, usize)> = None;
        let lo = bi.saturating_sub(stride);
        let hi = (bi + stride).min(feasible.len() - 1);
        for (k, p) in feasible.iter().enumerate().take(hi + 1).skip(lo) {
            if let Some(r) = solve(m, p, o) {
                if best.as_ref().is_none_or(|b| r.P_elec < b.0.P_elec) {
                    best = Some((r, k));
                }
            }
        }
        let Some((run, k)) = best else {
            out.unreachable.push(m);
            continue;
        };
        dots.push(Dot {
            motor: m,
            point: feasible[k].clone(),
            load: m.tq_cont_nm.filter(|&t| t != 0.0).map(|t| run.Q / t),
            run,
        });
    }
    let key = |d: &Dot| match rank_by {
        RankBy::Power => None,
        RankBy::Mass => Some(d.motor.mass_kg),
        RankBy::Price => Some(d.motor.price_usd),
    };
    let rank = |a: &Dot, b: &Dot| -> std::cmp::Ordering {
        let by_power = a.run.P_elec.total_cmp(&b.run.P_elec);
        match (key(a), key(b)) {
            (None, None) => by_power,
            (Some(av), Some(bv)) => match (av, bv) {
                (None, None) => by_power,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (Some(_), None) => std::cmp::Ordering::Less,
                (Some(x), Some(y)) => x.total_cmp(&y).then(by_power),
            },
            _ => by_power,
        }
    };
    dots.sort_by(rank);
    let shown = if all_windings {
        dots.clone()
    } else {
        let mut seen: Vec<Option<&str>> = Vec::new();
        let mut keep = Vec::new();
        for d in &dots {
            let fam = d.motor.family.as_deref();
            if !seen.contains(&fam) {
                seen.push(fam);
                keep.push(d.clone());
            }
        }
        keep.sort_by(rank);
        keep
    };
    out.dots = shown;
    out.all = dots;
    out
}

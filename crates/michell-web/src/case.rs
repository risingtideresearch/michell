//! A case as the queue stores it: what to compute on a saved hull, in a
//! canonical form, so that asking for the same case twice finds the first.
//!
//! The hull's import settings (waterline, stations, units, ...) belong to the
//! hull, not the case: a case is only the speed, closure, attitude, load and
//! layout. [`CaseParams::canonical`] fills in every default and checks every
//! range, so two requests that mean the same computation serialize — and
//! hash — alike.

use crate::{FlowRequest, LoftRequest, MassBy};
use michell::TransomClosure;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

/// The transom closure, as stored.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Closure {
    /// The free-fall hollow; `coeff` defaults to the solver's.
    Ballistic {
        #[serde(default)]
        coeff: Option<f64>,
    },
    /// A hollow of fixed length [m].
    Fixed {
        length: f64,
    },
    Off,
}

impl Default for Closure {
    fn default() -> Self {
        Closure::Ballistic { coeff: None }
    }
}

impl Closure {
    fn transom(self) -> TransomClosure {
        match self {
            Closure::Ballistic { coeff: Some(c) } => TransomClosure::Ballistic { coeff: c },
            Closure::Ballistic { coeff: None } => TransomClosure::default(),
            Closure::Fixed { length } => TransomClosure::Fixed { length },
            Closure::Off => TransomClosure::None,
        }
    }
}

/// What a case computes on its hull. `spans` makes it a fast catamaran span
/// sweep (one result per span, all at the middle span's attitude); otherwise
/// it is one flow, optionally with the hull doubled at `span`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CaseParams {
    /// Length Froude number on the longest hull.
    pub froude: f64,
    #[serde(default)]
    pub closure: Closure,
    /// Free-surface grid columns.
    #[serde(default = "default_grid")]
    pub grid: usize,
    /// Float at the dynamic equilibrium rather than the design attitude.
    #[serde(default = "default_dynamic")]
    pub dynamic: bool,
    /// Load [kg]; default the design displacement.
    #[serde(default)]
    pub mass: Option<f64>,
    /// LCG [m, fleet x]; default the LCB.
    #[serde(default)]
    pub lcg: Option<f64>,
    #[serde(default)]
    pub mass_by: MassBy,
    /// Catamaran centre span [m].
    #[serde(default)]
    pub span: Option<f64>,
    /// A span sweep's centre spans [m], ascending.
    #[serde(default)]
    pub spans: Option<Vec<f64>>,
}

fn default_grid() -> usize {
    640
}

fn default_dynamic() -> bool {
    true
}

impl CaseParams {
    /// `"flow"` or `"span_sweep"`.
    pub fn kind(&self) -> &'static str {
        if self.spans.is_some() {
            "span_sweep"
        } else {
            "flow"
        }
    }

    /// Checked, with every default filled in (so equal computations compare
    /// equal): the ballistic coefficient resolved, the grid clamped, the
    /// sweep's spans sorted and deduplicated.
    pub fn canonical(mut self) -> Result<CaseParams, String> {
        if !(self.froude > 0.0 && self.froude < 5.0) {
            return Err(format!("froude {}: expected 0 < Fn < 5", self.froude));
        }
        self.closure = match self.closure {
            Closure::Ballistic { coeff } => Closure::Ballistic {
                coeff: Some(match coeff {
                    Some(c) if c.is_finite() => c.max(0.0),
                    Some(c) => return Err(format!("closure coeff {c}: not a number")),
                    None => match TransomClosure::default() {
                        TransomClosure::Ballistic { coeff } => coeff,
                        _ => unreachable!("the default closure is ballistic"),
                    },
                }),
            },
            Closure::Fixed { length } if length.is_finite() => Closure::Fixed {
                length: length.max(0.0),
            },
            Closure::Fixed { length } => {
                return Err(format!("closure length {length}: not a number"))
            }
            Closure::Off => Closure::Off,
        };
        self.grid = self.grid.clamp(40, 1200);
        if self.mass.is_some_and(|m| !(m > 0.0 && m.is_finite())) {
            return Err("mass must be positive".into());
        }
        if self.lcg.is_some_and(|x| !x.is_finite()) {
            return Err("lcg: not a number".into());
        }
        if self.mass.is_none() && self.mass_by != MassBy::Sinking {
            // Scaling to the design displacement is no scaling at all.
            self.mass_by = MassBy::Sinking;
        }
        if self.span.is_some_and(|s| !(s > 0.0 && s.is_finite())) {
            return Err("span: expected a positive centre span".into());
        }
        if let Some(spans) = &mut self.spans {
            if self.span.is_some() {
                return Err("give span or spans, not both".into());
            }
            if spans.is_empty() || spans.iter().any(|s| !(*s > 0.0 && s.is_finite())) {
                return Err("spans: expected positive centre spans".into());
            }
            spans.sort_by(f64::total_cmp);
            spans.dedup();
        }
        Ok(self)
    }

    /// Hex SHA-256 of the canonical JSON: the case's identity on its hull.
    pub fn hash(&self) -> String {
        hex(&Sha256::digest(self.to_json().as_bytes()))
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(self).expect("params serialize")
    }

    /// The solver request, on the hull's cut, starting the equilibrium from
    /// `warm` when given.
    pub fn flow_request(&self, cut: LoftRequest, warm: Option<(f64, f64)>) -> FlowRequest {
        FlowRequest {
            cut,
            froude: self.froude,
            closure: self.closure.transom(),
            grid: self.grid,
            dynamic: self.dynamic,
            mass: self.mass,
            lcg: self.lcg,
            mass_by: self.mass_by,
            span: self.span,
            warm,
            hold: None,
        }
    }

    /// The same case at another speed: what a warm start may borrow from.
    pub fn same_but_speed(&self, other: &CaseParams) -> bool {
        CaseParams {
            froude: other.froude,
            grid: other.grid,
            ..self.clone()
        } == *other
    }
}

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> CaseParams {
        serde_json::from_str::<CaseParams>(s)
            .unwrap()
            .canonical()
            .unwrap()
    }

    /// Defaults spelt out and left out are the same case.
    #[test]
    fn equal_cases_hash_alike() {
        let a = parse(r#"{"froude": 0.3}"#);
        let b = parse(
            r#"{"froude": 0.3, "grid": 640, "dynamic": true, "mass_by": "sinking",
                "closure": {"type": "ballistic"}}"#,
        );
        assert_eq!(a, b);
        assert_eq!(a.hash(), b.hash());
        assert_ne!(a.hash(), parse(r#"{"froude": 0.31}"#).hash());
        let s = parse(r#"{"froude": 0.3, "spans": [3, 2, 3]}"#);
        assert_eq!(s.spans, Some(vec![2.0, 3.0]));
        assert_eq!(s.kind(), "span_sweep");
    }

    #[test]
    fn bad_cases_are_refused() {
        for s in [
            r#"{"froude": 0}"#,
            r#"{"froude": 0.3, "mass": -1}"#,
            r#"{"froude": 0.3, "span": 2, "spans": [2]}"#,
            r#"{"froude": 0.3, "spans": []}"#,
        ] {
            let p: CaseParams = serde_json::from_str(s).unwrap();
            assert!(p.canonical().is_err(), "{s}");
        }
        assert!(serde_json::from_str::<CaseParams>(r#"{"froude": 0.3, "speed": 2}"#).is_err());
    }
}

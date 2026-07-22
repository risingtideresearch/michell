//! Round-trip the example manifests from the README/repo through
//! load → model → save → load and check the model survives. Also checks that a
//! freshly serialized manifest re-parses (so `to_string` output stays loadable).

use michell_editor::model::*;
use michell_editor::{jsonio, validate};

// The ama study shipped in the repo root.
const AMA: &str = r#"{
  "name": "ama catamaran: weight x spacing",
  "fluid": "seawater",
  "hulls": [
    { "id": "ama_s", "file": "ama.hull", "pose": { "dy": 3.1580 } },
    { "id": "ama_p", "file": "ama.hull", "pose": { "dy": 0.6580 } }
  ],
  "sweep": [
    { "target": "speed", "unit": "knots", "range": [6, 10], "step": 2 },
    { "target": "weight", "range": [700, 1100], "step": 200 },
    { "target": ["ama_s", "ama_p"], "param": "spread", "range": [0, 0.75], "step": 0.375 }
  ],
  "output": { "format": "csv", "file": "ama-study.csv" }
}"#;

// The README placement study (lcg scalar, per-hull trim, options).
const PLACEMENT: &str = r#"{
  "name": "ama placement study",
  "fluid": "seawater",
  "hulls": [
    { "id": "vaka",  "file": "boat-center.hull" },
    { "id": "ama_s", "file": "boat-starboard.hull" },
    { "id": "ama_p", "file": "boat-port.hull", "pose": { "trim": 0.5 } }
  ],
  "sweep": [
    { "target": "speed", "unit": "knots", "range": [4, 10], "step": 0.5 },
    { "target": "weight", "range": [1800, 2600], "step": 200 },
    { "target": "lcg", "value": -5.8 },
    { "target": ["ama_s", "ama_p"], "param": "spread", "range": [1.5, 2.5] },
    { "target": "ama_s", "param": "trim", "values": [-2, 0, 2] }
  ],
  "output": { "format": "csv", "file": "study.csv" },
  "options": { "rel_tol": 1e-5, "form_factor": 0.05 }
}"#;

// The README GZ-curve example (heel + vcg).
const GZ: &str = r#"{
  "name": "gz",
  "fluid": "seawater",
  "hulls": [ { "id": "vaka", "file": "c.hull" } ],
  "sweep": [
    { "target": "speed", "unit": "knots", "value": 8 },
    { "target": "weight", "value": 2200 },
    { "target": "vcg", "value": 1.1 },
    { "target": "heel", "range": [-15, 15], "step": 1 }
  ]
}"#;

fn reparse(m: &Manifest) -> Manifest {
    let text = jsonio::to_string(m).expect("serialize");
    jsonio::from_str(&text).expect("reparse serialized output")
}

#[test]
fn ama_loads_and_survives_roundtrip() {
    let m = jsonio::from_str(AMA).unwrap();
    assert_eq!(m.name, "ama catamaran: weight x spacing");
    assert_eq!(m.fluid, Fluid::Seawater);
    assert_eq!(m.hulls.len(), 2);
    assert!(m.hulls[0].pose.enabled);
    assert_eq!(m.hulls[0].pose.dy, 3.158);
    assert_eq!(m.axes.len(), 3);
    assert_eq!(m.axes[0].kind, AxisKind::Speed);
    assert_eq!(m.axes[0].unit, SpeedUnit::Knots);
    assert_eq!(m.axes[2].kind, AxisKind::Pose);
    assert_eq!(m.axes[2].param, PoseParam::Spread);
    assert_eq!(m.axes[2].targets, vec!["ama_s".to_string(), "ama_p".into()]);
    assert!(validate::validate(&m)
        .iter()
        .all(|i| i.level != validate::Level::Error));

    let m2 = reparse(&m);
    assert_eq!(m2.hulls.len(), 2);
    assert_eq!(m2.axes.len(), 3);
    assert_eq!(m2.axes[2].targets.len(), 2);
    assert_eq!(m2.output.file, "ama-study.csv");
}

#[test]
fn placement_scalar_and_list_axes_survive() {
    let m = jsonio::from_str(PLACEMENT).unwrap();
    let lcg = m.axes.iter().find(|a| a.kind == AxisKind::Lcg).unwrap();
    assert_eq!(lcg.values.mode, ValueMode::Scalar);
    assert_eq!(lcg.values.scalar, "-5.8");
    let trim = m
        .axes
        .iter()
        .find(|a| a.kind == AxisKind::Pose && a.param == PoseParam::Trim)
        .unwrap();
    assert_eq!(trim.values.mode, ValueMode::List);
    assert!(m.options.rel_tol.enabled);
    assert!(m.options.form_factor.enabled);
    assert!(!m.options.gravity.enabled);

    let m2 = reparse(&m);
    let lcg2 = m2.axes.iter().find(|a| a.kind == AxisKind::Lcg).unwrap();
    assert_eq!(lcg2.values.scalar, "-5.8");
    assert!(m2.options.rel_tol.enabled);
}

#[test]
fn gz_heel_requires_vcg_and_validates() {
    let m = jsonio::from_str(GZ).unwrap();
    assert!(m.axes.iter().any(|a| a.kind == AxisKind::Heel));
    assert!(m.axes.iter().any(|a| a.kind == AxisKind::Vcg));
    let issues = validate::validate(&m);
    assert!(
        issues.iter().all(|i| i.level != validate::Level::Error),
        "unexpected: {:?}",
        issues.iter().map(|i| &i.msg).collect::<Vec<_>>()
    );
}

#[test]
fn validator_flags_missing_speed_and_lonely_lcg() {
    let mut m = Manifest::default();
    m.hulls.push(HullSpec {
        id: "a".into(),
        file: "a.hull".into(),
        pose: Pose::default(),
    });
    m.axes = vec![AxisSpec::new(AxisKind::Lcg)]; // no speed, no weight
    let issues = validate::validate(&m);
    let msgs: Vec<&str> = issues.iter().map(|i| i.msg.as_str()).collect();
    assert!(
        msgs.iter().any(|s| s.contains("needs a speed axis")),
        "{msgs:?}"
    );
    assert!(
        msgs.iter()
            .any(|s| s.contains("lcg axis requires a weight")),
        "{msgs:?}"
    );
}

#[test]
fn whole_numbers_serialize_without_trailing_zero() {
    let m = jsonio::from_str(AMA).unwrap();
    let text = jsonio::to_string(&m).unwrap();
    // weight range 700..1100 should read as integers, not 700.0.
    assert!(text.contains("700"));
    assert!(!text.contains("700.0"));
}

//! Round-trip the example manifests from the README through
//! load → model → save → load and check the model survives. Also checks that a
//! freshly serialized manifest re-parses (so `to_string` output stays loadable).

use michell_editor::model::*;
use michell_editor::{jsonio, validate};

// The README placement study: per-hull loads, a point load, scale + mass axes.
const PLACEMENT: &str = r#"{
  "name": "ama placement study",
  "fluid": "seawater",
  "hulls": [
    { "id": "vaka",  "file": "boat-center.hull",
      "load": { "mass": 1800, "lcg": -5.8, "vcg": 1.1 },
      "points": [
        { "id": "battery", "mass": 200, "dx": 1.0, "dz": 0.6 },
        { "id": "crew",    "mass": 160, "dx": -2.0, "dz": -0.4 }
      ] },
    { "id": "ama_s", "file": "boat-starboard.hull", "load": { "mass": 120, "vcg": 0.4 } },
    { "id": "ama_p", "file": "boat-port.hull", "pose": { "trim": 0.5 }, "load": { "mass": 120, "vcg": 0.4 } }
  ],
  "sweep": [
    { "target": "speed", "unit": "knots", "range": [4, 10], "step": 0.5 },
    { "target": "vaka", "param": "mass", "range": [0, 800], "step": 200 },
    { "target": "battery", "param": "dz", "values": [0.0, 0.6, 1.2] },
    { "target": ["ama_s", "ama_p"], "param": "spread", "range": [1.5, 2.5] },
    { "target": "ama_s", "param": "trim", "values": [-2, 0, 2] },
    { "target": "vaka", "param": "scale", "values": [0.9, 1.0, 1.1] }
  ],
  "output": { "format": "csv", "file": "study.csv" },
  "options": { "rel_tol": 1e-5, "form_factor": 0.05 }
}"#;

// The README heel-rollup example: symmetric loads, options.heel, no CG axis.
const HEEL: &str = r#"{
  "name": "gz",
  "fluid": "seawater",
  "options": { "heel": { "resistance_angles": [5, 10], "gz_step": 2.5, "gz_max": 90 } },
  "hulls": [
    { "id": "port", "file": "boat-port.hull", "load": { "mass": 1100, "vcg": 1.1 } },
    { "id": "stbd", "file": "boat-starboard.hull", "load": { "mass": 1100, "vcg": 1.1 } }
  ],
  "sweep": [
    { "target": "speed", "unit": "knots", "value": 8 }
  ]
}"#;

fn reparse(m: &Manifest) -> Manifest {
    let text = jsonio::to_string(m).expect("serialize");
    jsonio::from_str(&text).expect("reparse serialized output")
}

fn no_errors(m: &Manifest) -> bool {
    validate::validate(m)
        .iter()
        .all(|i| i.level != validate::Level::Error)
}

#[test]
fn placement_loads_points_and_axes_survive() {
    let m = jsonio::from_str(PLACEMENT).unwrap();
    assert_eq!(m.hulls.len(), 3);

    // vaka: load with lcg set, two point loads.
    let vaka = &m.hulls[0];
    assert!(vaka.load.enabled);
    assert_eq!(vaka.load.mass, 1800.0);
    assert!(vaka.load.lcg_set);
    assert_eq!(vaka.load.lcg, -5.8);
    assert_eq!(vaka.points.len(), 2);
    assert_eq!(vaka.points[0].id, "battery");
    assert_eq!(vaka.points[0].dz, 0.6);

    // ama_s load without lcg.
    assert!(m.hulls[1].load.enabled);
    assert!(!m.hulls[1].load.lcg_set);

    // Axes: a hull mass axis, a point dz axis, a coupled spread, a scale axis.
    let mass = m
        .axes
        .iter()
        .find(|a| a.kind == AxisKind::Hull && a.hull_param == HullParam::Mass)
        .unwrap();
    assert_eq!(mass.targets, vec!["vaka".to_string()]);
    let batt = m.axes.iter().find(|a| a.kind == AxisKind::Point).unwrap();
    assert_eq!(batt.point_param, PointParam::Dz);
    assert_eq!(batt.targets, vec!["battery".to_string()]);
    let spread = m
        .axes
        .iter()
        .find(|a| a.hull_param == HullParam::Spread)
        .unwrap();
    assert_eq!(spread.targets.len(), 2);
    assert!(m
        .axes
        .iter()
        .any(|a| a.hull_param == HullParam::Scale && a.kind == AxisKind::Hull));

    assert!(
        no_errors(&m),
        "{:?}",
        validate::validate(&m)
            .iter()
            .map(|i| &i.msg)
            .collect::<Vec<_>>()
    );

    let m2 = reparse(&m);
    assert_eq!(m2.hulls[0].points.len(), 2);
    assert_eq!(m2.axes.len(), m.axes.len());
    let batt2 = m2.axes.iter().find(|a| a.kind == AxisKind::Point).unwrap();
    assert_eq!(batt2.targets, vec!["battery".to_string()]);
    assert!(no_errors(&m2));
}

#[test]
fn heel_rollup_options_survive() {
    let m = jsonio::from_str(HEEL).unwrap();
    assert!(m.options.heel.enabled);
    assert_eq!(m.options.heel.resistance_angles, "5, 10");
    assert_eq!(m.options.heel.gz_step, "2.5");
    // Both hulls carry mass → float mode → heel metrics valid, no errors.
    assert!(
        no_errors(&m),
        "{:?}",
        validate::validate(&m)
            .iter()
            .map(|i| &i.msg)
            .collect::<Vec<_>>()
    );

    let m2 = reparse(&m);
    assert!(m2.options.heel.enabled);
    assert_eq!(m2.options.heel.resistance_angles, "5, 10");
}

#[test]
fn scale_defaults_to_one_and_omitted() {
    // A hull with a pose but no scale should serialize without a scale key.
    let mut m = Manifest::default();
    m.hulls.push(HullSpec {
        id: "h".into(),
        file: "h.hull".into(),
        pose: Pose {
            enabled: true,
            dz: 0.1,
            ..Default::default()
        },
        ..Default::default()
    });
    let text = jsonio::to_string(&m).unwrap();
    assert!(
        !text.contains("scale"),
        "default scale should be omitted:\n{text}"
    );
    assert_eq!(m.hulls[0].pose.scale, 1.0);
}

#[test]
fn validator_flags_cg_axis_without_mass() {
    let mut m = Manifest::default();
    m.hulls.push(HullSpec {
        id: "a".into(),
        file: "a.hull".into(),
        ..Default::default()
    });
    // A vcg axis but no mass anywhere.
    let mut speed = AxisSpec::new(AxisKind::Speed);
    speed.values.mode = ValueMode::Scalar;
    speed.values.scalar = "5".into();
    let mut vcg = AxisSpec::new(AxisKind::Hull);
    vcg.hull_param = HullParam::Vcg;
    vcg.targets = vec!["a".into()];
    m.axes = vec![speed, vcg];
    let msgs: Vec<String> = validate::validate(&m).into_iter().map(|i| i.msg).collect();
    assert!(
        msgs.iter()
            .any(|s| s.contains("needs the fleet to carry mass")),
        "{msgs:?}"
    );
}

#[test]
fn validator_flags_missing_speed() {
    let mut m = Manifest::default();
    m.hulls.push(HullSpec {
        id: "a".into(),
        file: "a.hull".into(),
        ..Default::default()
    });
    m.axes = vec![]; // no speed
    let msgs: Vec<String> = validate::validate(&m).into_iter().map(|i| i.msg).collect();
    assert!(
        msgs.iter().any(|s| s.contains("needs a speed axis")),
        "{msgs:?}"
    );
}

#[test]
fn whole_numbers_serialize_without_trailing_zero() {
    let m = jsonio::from_str(PLACEMENT).unwrap();
    let text = jsonio::to_string(&m).unwrap();
    assert!(text.contains("1800"));
    assert!(!text.contains("1800.0"));
}

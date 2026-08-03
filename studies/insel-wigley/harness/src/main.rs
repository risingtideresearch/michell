use michell::{
    hulls, multihull_wave_resistance_with, Conditions, Fluid, Hull, Placement, WaveOptions,
    WaveResistance,
};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

const LENGTH: f64 = 1.8;
const BEAM: f64 = 0.18;
const DRAFT: f64 = 0.1125;
const DEMIHULL_WETTED_SURFACE: f64 = 0.482;
const DENSITY: f64 = 1000.0;
const KINEMATIC_VISCOSITY: f64 = 1.141e-6;
const GRAVITY: f64 = 9.80665;
const REL_TOL: f64 = 1e-6;
const MAX_REFINEMENTS: usize = 8;
const FN_SCALE: f64 = 10_000_000.0;

const FIXED_GRID: u8 = 1;
const FREE_GRID: u8 = 2;
const FINE_GRID: u8 = 4;

#[derive(Clone, Copy)]
struct Configuration {
    name: &'static str,
    file_suffix: &'static str,
    separation_over_length: Option<f64>,
}

const CONFIGURATIONS: [Configuration; 5] = [
    Configuration {
        name: "monohull",
        file_suffix: "monohull",
        separation_over_length: None,
    },
    Configuration {
        name: "s_l_0_2",
        file_suffix: "s_l_0_2",
        separation_over_length: Some(0.2),
    },
    Configuration {
        name: "s_l_0_3",
        file_suffix: "s_l_0_3",
        separation_over_length: Some(0.3),
    },
    Configuration {
        name: "s_l_0_4",
        file_suffix: "s_l_0_4",
        separation_over_length: Some(0.4),
    },
    Configuration {
        name: "s_l_0_5",
        file_suffix: "s_l_0_5",
        separation_over_length: Some(0.5),
    },
];

fn fn_key(value: f64) -> i64 {
    (value * FN_SCALE).round() as i64
}

fn fn_value(key: i64) -> f64 {
    key as f64 / FN_SCALE
}

fn read_cwp_grid(path: &Path, source: u8, grid: &mut BTreeMap<i64, u8>) {
    let file = File::open(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let mut header = None;
    for line in BufReader::new(file).lines() {
        let line = line.unwrap();
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if header.is_none() {
            header = Some(line.split(',').map(str::to_owned).collect::<Vec<String>>());
            continue;
        }
        let fields = line.split(',').collect::<Vec<_>>();
        let names = header.as_ref().unwrap();
        let observable = fields[names.iter().position(|name| name == "observable").unwrap()];
        if observable != "cwp" {
            continue;
        }
        let fn_ = fields[names.iter().position(|name| name == "fn").unwrap()]
            .parse::<f64>()
            .unwrap();
        *grid.entry(fn_key(fn_)).or_default() |= source;
    }
}

fn grid_for(data: &Path, configuration: Configuration) -> BTreeMap<i64, u8> {
    let mut grid = BTreeMap::new();
    for (attitude, source) in [("fixed", FIXED_GRID), ("free", FREE_GRID)] {
        read_cwp_grid(
            &data.join(format!("{attitude}_{}.csv", configuration.file_suffix)),
            source,
            &mut grid,
        );
    }
    for index in 0..=60 {
        let fn_ = 0.25 + 0.005 * index as f64;
        *grid.entry(fn_key(fn_)).or_default() |= FINE_GRID;
    }
    grid
}

fn grid_source(mask: u8) -> String {
    let mut names = Vec::new();
    if mask & FIXED_GRID != 0 {
        names.push("measured_fixed");
    }
    if mask & FREE_GRID != 0 {
        names.push("measured_free");
    }
    if mask & FINE_GRID != 0 {
        names.push("fine");
    }
    names.join("+")
}

fn conditions(fn_: f64) -> Conditions {
    Conditions {
        speed: fn_ * (GRAVITY * LENGTH).sqrt(),
        fluid: Fluid {
            density: DENSITY,
            kinematic_viscosity: KINEMATIC_VISCOSITY,
        },
        gravity: GRAVITY,
    }
}

fn require_converged(label: &str, fn_: f64, result: WaveResistance) -> WaveResistance {
    assert!(
        result.outcome.is_converged(),
        "{label} Fn={fn_:.7}: outcome={}, est_rel_error={:.3e}, evaluations={}",
        result.outcome.as_str(),
        result.est_rel_error,
        result.inner_evaluations,
    );
    assert!(result.resistance.is_finite() && result.resistance >= 0.0);
    assert!(result.est_rel_error.is_finite() && result.est_rel_error >= 0.0);
    result
}

fn solo_result(
    hull: &Hull,
    fn_: f64,
    options: &WaveOptions,
    cache: &mut BTreeMap<i64, WaveResistance>,
) -> WaveResistance {
    let key = fn_key(fn_);
    if let Some(result) = cache.get(&key) {
        return *result;
    }
    let result = require_converged(
        "standalone multihull-path solve",
        fn_,
        multihull_wave_resistance_with(&[(hull, Placement::default())], &conditions(fn_), options)
            .unwrap(),
    );
    cache.insert(key, result);
    result
}

fn output_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../data/predictions/predictions.csv")
}

fn main() {
    let study = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let digitized = study.join("data/digitized");
    let output = output_path();
    fs::create_dir_all(output.parent().unwrap()).unwrap();

    let hull = hulls::wigley(LENGTH, BEAM, DRAFT).expect("valid exact Wigley geometry");
    let options = WaveOptions {
        rel_tol: REL_TOL,
        max_refinements: MAX_REFINEMENTS,
    };
    let mut solo_cache = BTreeMap::new();
    let mut writer = BufWriter::new(File::create(&output).unwrap());
    writeln!(
        writer,
        "# geometry: C2 Wigley L=1.8m B=0.18m T=0.1125m; exact biquadratic library constructor"
    )
    .unwrap();
    writeln!(
        writer,
        "# fluid: rho=1000kg/m3 nu=1.141e-6m2/s g=9.80665m/s2; model scale"
    )
    .unwrap();
    writeln!(writer, "# solver: multihull_wave_resistance_with rel_tol=1e-6 max_refinements=8; all outcomes required converged").unwrap();
    writeln!(
        writer,
        "configuration,grid_source,fn,speed_m_s,member_count,separation_over_length,separation_m,wetted_surface_m2,rw_n,cw,interference,solo_rw_n,method,outcome,est_rel_error,evaluations,max_lambda,solo_method,solo_outcome,solo_est_rel_error,solo_evaluations,solo_max_lambda,rel_tol"
    )
    .unwrap();

    let mut rows = 0usize;
    for configuration in CONFIGURATIONS {
        for (key, source) in grid_for(&digitized, configuration) {
            let fn_ = fn_value(key);
            let condition = conditions(fn_);
            let solo = solo_result(&hull, fn_, &options, &mut solo_cache);
            let (member_count, separation, area, result, interference) =
                if let Some(separation_over_length) = configuration.separation_over_length {
                    let separation = separation_over_length * LENGTH;
                    let members = [
                        (
                            &hull,
                            Placement {
                                x: 0.0,
                                y: -separation / 2.0,
                            },
                        ),
                        (
                            &hull,
                            Placement {
                                x: 0.0,
                                y: separation / 2.0,
                            },
                        ),
                    ];
                    let result = require_converged(
                        configuration.name,
                        fn_,
                        multihull_wave_resistance_with(&members, &condition, &options).unwrap(),
                    );
                    (
                        2,
                        separation,
                        2.0 * DEMIHULL_WETTED_SURFACE,
                        result,
                        result.resistance / (2.0 * solo.resistance),
                    )
                } else {
                    (1, 0.0, DEMIHULL_WETTED_SURFACE, solo, 1.0)
                };
            assert!(interference.is_finite() && interference >= 0.0);
            let dynamic_pressure_area = 0.5 * DENSITY * condition.speed * condition.speed * area;
            let cw = result.resistance / dynamic_pressure_area;
            writeln!(
                writer,
                "{},{},{:.7},{:.12},{},{},{:.12},{:.12},{:.12e},{:.12e},{:.12e},{:.12e},{},{},{:.12e},{},{:.12e},{},{},{:.12e},{},{:.12e},{:.1e}",
                configuration.name,
                grid_source(source),
                fn_,
                condition.speed,
                member_count,
                configuration
                    .separation_over_length
                    .map(|value| format!("{value:.1}"))
                    .unwrap_or_default(),
                separation,
                area,
                result.resistance,
                cw,
                interference,
                solo.resistance,
                result.method.as_str(),
                result.outcome.as_str(),
                result.est_rel_error,
                result.inner_evaluations,
                result.max_lambda,
                solo.method.as_str(),
                solo.outcome.as_str(),
                solo.est_rel_error,
                solo.inner_evaluations,
                solo.max_lambda,
                REL_TOL,
            )
            .unwrap();
            rows += 1;
        }
    }
    writer.flush().unwrap();
    eprintln!(
        "wrote {rows} converged rows to {} ({} unique standalone solves)",
        output.display(),
        solo_cache.len()
    );
}

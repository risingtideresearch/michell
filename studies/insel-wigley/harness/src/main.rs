use insel_wigley_harness::canal::{
    converged_resistance, converged_resistance_variant, Canal, Flow, HistoricalVariant,
    ModalOptions, WigleyHull,
};
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
const THEORY_GRID: u8 = 8;

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

fn theory_grid() -> BTreeMap<i64, u8> {
    (0..=150)
        .map(|index| (fn_key(0.20 + 0.005 * index as f64), THEORY_GRID))
        .collect()
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
    if mask & THEORY_GRID != 0 {
        names.push("theory");
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

fn output_path(theory: bool) -> PathBuf {
    let filename = if theory {
        "theory_predictions.csv"
    } else {
        "predictions.csv"
    };
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../data/predictions")
        .join(filename)
}

fn run_canal() {
    let output =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../data/predictions/canal_predictions.csv");
    fs::create_dir_all(output.parent().unwrap()).unwrap();
    let hull = WigleyHull {
        length: LENGTH,
        beam: BEAM,
        draft: DRAFT,
    };
    let canal = Canal {
        width: 3.7,
        depth: 1.85,
    };
    let options = ModalOptions::default();
    let mut writer = BufWriter::new(File::create(&output).unwrap());
    writeln!(
        writer,
        "# geometry: C2 Wigley L=1.8m B=0.18m T=0.1125m; continuous analytic source amplitude"
    )
    .unwrap();
    writeln!(
        writer,
        "# fluid: rho=1000kg/m3 g=9.80665m/s2; model scale; canal W=3.7m H=1.85m"
    )
    .unwrap();
    writeln!(writer, "# solver: Insel equations 4.25-4.50; modes doubled from 32 through 1048576; resistance_rel_tol=5e-6 interference_abs_tol=2e-6").unwrap();
    writeln!(writer, "configuration,fn,speed_m_s,member_count,separation_over_length,separation_m,wetted_surface_m2,rw_n,cw,interference,solo_rw_n,modes,resistance_rel_change,interference_abs_change,method,outcome").unwrap();

    let mut rows = 0usize;
    for configuration in CONFIGURATIONS {
        let separation = configuration.separation_over_length.unwrap_or(0.0) * LENGTH;
        for index in 0..=150 {
            let fn_ = 0.20 + 0.005 * index as f64;
            let speed = fn_ * (GRAVITY * LENGTH).sqrt();
            let result = converged_resistance(
                hull,
                Flow {
                    speed,
                    density: DENSITY,
                    gravity: GRAVITY,
                },
                canal,
                separation,
                options,
            )
            .unwrap_or_else(|error| panic!("{} Fn={fn_:.3}: {error}", configuration.name));
            let (members, area, resistance, interference) =
                if configuration.separation_over_length.is_some() {
                    (
                        2,
                        2.0 * DEMIHULL_WETTED_SURFACE,
                        result.catamaran_resistance,
                        result.interference,
                    )
                } else {
                    (1, DEMIHULL_WETTED_SURFACE, result.monohull_resistance, 1.0)
                };
            let cw = resistance / (0.5 * DENSITY * speed.powi(2) * area);
            writeln!(
                writer,
                "{},{fn_:.7},{speed:.12},{members},{},{separation:.12},{area:.12},{resistance:.12e},{cw:.12e},{interference:.12e},{:.12e},{},{:.12e},{:.12e},insel_finite_canal_modal_continuous_wigley,converged",
                configuration.name,
                configuration
                    .separation_over_length
                    .map(|value| format!("{value:.1}"))
                    .unwrap_or_default(),
                result.monohull_resistance,
                result.modes,
                result.resistance_rel_change,
                result.interference_abs_change,
            )
            .unwrap();
            rows += 1;
        }
    }
    writer.flush().unwrap();
    eprintln!("wrote {rows} converged canal rows to {}", output.display());
}

fn run_separation_grid() {
    let output = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../data/predictions/separation_grid_predictions.csv");
    fs::create_dir_all(output.parent().unwrap()).unwrap();
    let hull = hulls::wigley(LENGTH, BEAM, DRAFT).expect("valid exact Wigley geometry");
    let options = WaveOptions {
        rel_tol: REL_TOL,
        max_refinements: 10,
    };
    let mut solo_cache = BTreeMap::new();
    let mut writer = BufWriter::new(File::create(&output).unwrap());
    writeln!(writer, "# criteria: CRITERIA-SEPARATION.md").unwrap();
    writeln!(
        writer,
        "# geometry: exact C2 Wigley L=1.8m B=0.18m T=0.1125m; centreline separation"
    )
    .unwrap();
    writeln!(
        writer,
        "# grid: S/L=0.080:0.005:0.550; Fn=0.150:0.005:1.000; rel_tol=1e-6; max_refinements=10"
    )
    .unwrap();
    writeln!(writer, "separation_over_length,separation_m,fn,speed_m_s,interference,pair_rw_n,solo_rw_n,pair_method,pair_outcome,pair_est_rel_error,pair_evaluations,pair_max_lambda,solo_method,solo_outcome,solo_est_rel_error,solo_evaluations,solo_max_lambda,rel_tol").unwrap();
    let mut rows = 0usize;
    for separation_index in 0..=94 {
        let separation_over_length = 0.080 + 0.005 * separation_index as f64;
        let separation = separation_over_length * LENGTH;
        for fn_index in 0..=170 {
            let fn_ = 0.150 + 0.005 * fn_index as f64;
            let condition = conditions(fn_);
            let solo = solo_result(&hull, fn_, &options, &mut solo_cache);
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
            let pair = require_converged(
                "separation grid pair",
                fn_,
                multihull_wave_resistance_with(&members, &condition, &options).unwrap(),
            );
            let interference = pair.resistance / (2.0 * solo.resistance);
            assert!(interference.is_finite() && interference >= 0.0);
            writeln!(
                writer,
                "{separation_over_length:.3},{separation:.12},{fn_:.7},{:.12},{interference:.12e},{:.12e},{:.12e},{},{},{:.12e},{},{:.12e},{},{},{:.12e},{},{:.12e},{:.1e}",
                condition.speed,
                pair.resistance,
                solo.resistance,
                pair.method.as_str(),
                pair.outcome.as_str(),
                pair.est_rel_error,
                pair.inner_evaluations,
                pair.max_lambda,
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
        "wrote {rows} converged separation-grid rows to {} ({} unique standalone solves)",
        output.display(),
        solo_cache.len()
    );
}

fn run_historical_variants() {
    let output = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../data/predictions/historical_variant_predictions.csv");
    fs::create_dir_all(output.parent().unwrap()).unwrap();
    let hull = WigleyHull {
        length: LENGTH,
        beam: BEAM,
        draft: DRAFT,
    };
    let canal = Canal {
        width: 3.7,
        depth: 1.85,
    };
    let options = ModalOptions::default();
    let variants = [
        HistoricalVariant::Baseline,
        HistoricalVariant::LiteralEquation429,
        HistoricalVariant::DoubledCrossTerm,
        HistoricalVariant::HalfNonzeroModeMultiplicity,
        HistoricalVariant::DoubleNonzeroModeMultiplicity,
        HistoricalVariant::CosineInsteadOfCosineSquared,
    ];
    let mut writer = BufWriter::new(File::create(&output).unwrap());
    writeln!(writer, "# criteria: time-boxed historical-variant probe requested after CRITERIA-SEPARATION.md was refuted").unwrap();
    writeln!(writer, "# geometry: exact C2 Wigley L=1.8m B=0.18m T=0.1125m; canal W=3.7m H=1.85m; face-value centreline separations").unwrap();
    writeln!(writer, "# grid: S/L=0.2,0.3,0.4,0.5; Fn=0.20:0.005:0.95; every variant must pass unchanged CRITERIA-THEORY.md gates on all four panels").unwrap();
    writeln!(writer, "variant,configuration,fn,speed_m_s,separation_over_length,separation_m,interference,catamaran_rw_n,monohull_rw_n,modes,resistance_rel_change,interference_abs_change,outcome").unwrap();

    let mut rows = 0usize;
    for variant in variants {
        for configuration in &CONFIGURATIONS[1..] {
            let separation_over_length = configuration.separation_over_length.unwrap();
            let separation = separation_over_length * LENGTH;
            for fn_index in 0..=150 {
                let fn_ = 0.20 + 0.005 * fn_index as f64;
                let speed = fn_ * (GRAVITY * LENGTH).sqrt();
                let result = converged_resistance_variant(
                    hull,
                    Flow {
                        speed,
                        density: DENSITY,
                        gravity: GRAVITY,
                    },
                    canal,
                    separation,
                    options,
                    variant,
                )
                .unwrap_or_else(|error| {
                    panic!(
                        "{} {} Fn={fn_:.3}: {error}",
                        variant.as_str(),
                        configuration.name
                    )
                });
                writeln!(
                    writer,
                    "{},{},{fn_:.7},{speed:.12},{separation_over_length:.1},{separation:.12},{:.12e},{:.12e},{:.12e},{},{:.12e},{:.12e},converged",
                    variant.as_str(),
                    configuration.name,
                    result.interference,
                    result.catamaran_resistance,
                    result.monohull_resistance,
                    result.modes,
                    result.resistance_rel_change,
                    result.interference_abs_change,
                )
                .unwrap();
                rows += 1;
            }
        }
    }
    writer.flush().unwrap();
    eprintln!(
        "wrote {rows} converged historical-variant rows to {}",
        output.display()
    );
}

fn main() {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    if arguments.as_slice() == ["--canal"] {
        run_canal();
        return;
    }
    if arguments.as_slice() == ["--separation-grid"] {
        run_separation_grid();
        return;
    }
    if arguments.as_slice() == ["--historical-variants"] {
        run_historical_variants();
        return;
    }
    let theory = match arguments.as_slice() {
        [] => false,
        [flag] if flag == "--theory" => true,
        _ => panic!("usage: insel-wigley-harness [--theory | --canal | --separation-grid | --historical-variants]"),
    };
    let study = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let digitized = study.join("data/digitized");
    let output = output_path(theory);
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
    let configurations = if theory {
        &CONFIGURATIONS[1..]
    } else {
        &CONFIGURATIONS[..]
    };
    for &configuration in configurations {
        let grid = if theory {
            theory_grid()
        } else {
            grid_for(&digitized, configuration)
        };
        for (key, source) in grid {
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

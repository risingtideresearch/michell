//! Physics checks for the free-wave spectrum and wake reconstruction.

use michell::iges::{self, HullPose, Platform, SectionalOptions};
use michell::{sectional, Conditions, FreeWaveSpectrum, Placement, SectionalHull, TransomClosure};

/// The Wigley hull cut into sections from its exact CAD surfaces.
fn wigley(length: f64, beam: f64, draft: f64) -> SectionalHull {
    let surfaces = iges::wigley_surfaces(length, beam, draft).unwrap();
    let src = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
    src.situate_sectional(
        0,
        0.0,
        &HullPose::default(),
        &Platform::default(),
        &SectionalOptions::default(),
    )
    .unwrap()
    .expect("wet")
    .hull
}

fn spectrum<'h>(
    members: &[(&'h SectionalHull, Placement)],
    cond: &Conditions,
) -> FreeWaveSpectrum<'h> {
    FreeWaveSpectrum::new_sectional(members, cond, TransomClosure::default()).unwrap()
}

/// Integrate dR/dθ over (−π/2, π/2) by fine trapezoid.
fn resistance_from_spectrum(spec: &mut FreeWaveSpectrum, n: usize) -> f64 {
    let lim = 89.9f64.to_radians();
    let h = 2.0 * lim / n as f64;
    let mut sum = 0.0;
    for i in 0..=n {
        let theta = -lim + i as f64 * h;
        let w = if i == 0 || i == n { 0.5 } else { 1.0 };
        sum += w * spec.resistance_density(theta);
    }
    sum * h
}

#[test]
fn spectrum_reproduces_michell_resistance_monohull() {
    let hull = wigley(10.0, 1.0, 0.625);
    let cond = Conditions::seawater(3.0);
    let members = [(&hull, Placement::default())];
    let rw = sectional::wave_resistance(&hull, &cond, &Default::default())
        .unwrap()
        .resistance;
    let mut spec = spectrum(&members, &cond);
    let rw_spec = resistance_from_spectrum(&mut spec, 4_000);
    assert!(
        (rw_spec - rw).abs() <= 1e-3 * rw,
        "spectrum {rw_spec} vs michell {rw} (rel {:.1e})",
        (rw_spec - rw).abs() / rw
    );
}

#[test]
fn spectrum_reproduces_michell_resistance_staggered_catamaran() {
    let hull = wigley(8.0, 0.8, 0.5);
    let cond = Conditions::seawater(2.5);
    let members = [
        (&hull, Placement { x: 0.0, y: 1.4 }),
        (&hull, Placement { x: 1.7, y: -1.4 }),
    ];
    let rw = sectional::multihull_wave_resistance(&members, &cond, &Default::default())
        .unwrap()
        .resistance;
    let mut spec = spectrum(&members, &cond);
    let rw_spec = resistance_from_spectrum(&mut spec, 8_000);
    assert!(
        (rw_spec - rw).abs() <= 2e-3 * rw,
        "spectrum {rw_spec} vs michell {rw} (rel {:.1e})",
        (rw_spec - rw).abs() / rw
    );
}

#[test]
fn transverse_wavelength_far_astern() {
    let hull = wigley(10.0, 1.0, 0.625);
    let cond = Conditions::seawater(3.0);
    let members = [(&hull, Placement::default())];
    let mut spec = spectrum(&members, &cond);
    let want = spec.transverse_wavelength(); // 2π U²/g ≈ 5.766 m

    // Longitudinal cut on the track, 10–45 m astern of the stern.
    let (x0, x1, n) = (-50.0, -15.0, 3501);
    let grid = spec.elevation_grid(x0, x1, 0.0, 0.0, n, 1).unwrap();
    assert!(!grid.resolution_limited);
    let peak = grid.zeta.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    assert!(peak > 1e-3, "wake unexpectedly flat: peak {peak} m");

    // Mean spacing of zero crossings = half a wavelength.
    let mut crossings = Vec::new();
    for i in 1..n {
        let (a, b) = (grid.zeta[i - 1], grid.zeta[i]);
        if a == 0.0 || a.signum() == b.signum() {
            continue;
        }
        let x = grid.x(i - 1) + (grid.x(i) - grid.x(i - 1)) * a / (a - b);
        crossings.push(x);
    }
    assert!(
        crossings.len() >= 8,
        "too few crossings: {}",
        crossings.len()
    );
    let spacing = (crossings.last().unwrap() - crossings[0]) / (crossings.len() - 1) as f64;
    let got = 2.0 * spacing;
    assert!(
        (got - want).abs() <= 0.02 * want,
        "wavelength {got} vs 2πU²/g = {want}"
    );
}

#[test]
fn pattern_translates_with_placement() {
    let hull = wigley(10.0, 1.0, 0.625);
    let cond = Conditions::seawater(3.0);
    let base_members = [(&hull, Placement::default())];
    let shifted_members = [(&hull, Placement { x: 7.0, y: 0.0 })];
    let mut base = spectrum(&base_members, &cond);
    let mut shifted = spectrum(&shifted_members, &cond);
    for (x, y) in [(-20.0, 0.0), (-27.3, 0.0), (-33.0, 0.0), (-25.0, 4.0)] {
        let a = base.elevation_at(x, y).unwrap();
        let b = shifted.elevation_at(x + 7.0, y).unwrap();
        assert!(
            (a - b).abs() <= 1e-9,
            "at ({x}, {y}): base {a} vs shifted {b}"
        );
    }
}

#[test]
fn monohull_pattern_is_y_symmetric() {
    let hull = wigley(10.0, 1.0, 0.625);
    let cond = Conditions::seawater(3.0);
    let members = [(&hull, Placement::default())];
    let mut spec = spectrum(&members, &cond);
    let (nx, ny) = (40, 25);
    let grid = spec
        .elevation_grid(-30.0, -10.0, -12.0, 12.0, nx, ny)
        .unwrap();
    let peak = grid.zeta.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    for iy in 0..ny / 2 {
        for ix in 0..nx {
            let a = grid.get(ix, iy);
            let b = grid.get(ix, ny - 1 - iy);
            assert!(
                (a - b).abs() <= 1e-12 * peak.max(1e-12),
                "asymmetry at ix={ix} iy={iy}: {a} vs {b}"
            );
        }
    }
}

#[test]
fn fleet_pattern_superposes() {
    let hull = wigley(10.0, 1.0, 0.625);
    let cond = Conditions::seawater(3.0);
    let pa = Placement { x: 1.2, y: 2.5 };
    let pb = Placement { x: 0.0, y: -2.5 };
    let (x0, x1, y0, y1, nx, ny) = (-35.0, -8.0, -10.0, 10.0, 60, 41);

    let cat_members = [(&hull, pa), (&hull, pb)];
    let mut cat = spectrum(&cat_members, &cond);
    let g = cat.elevation_grid(x0, x1, y0, y1, nx, ny).unwrap();

    let solo_a_members = [(&hull, pa)];
    let solo_b_members = [(&hull, pb)];
    let mut sa = spectrum(&solo_a_members, &cond);
    let mut sb = spectrum(&solo_b_members, &cond);
    let ga = sa.elevation_grid(x0, x1, y0, y1, nx, ny).unwrap();
    let gb = sb.elevation_grid(x0, x1, y0, y1, nx, ny).unwrap();

    let peak = g.zeta.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    for i in 0..g.zeta.len() {
        let sum = ga.zeta[i] + gb.zeta[i];
        assert!(
            (g.zeta[i] - sum).abs() <= 1e-4 * peak,
            "at {i}: fleet {} vs sum {sum}",
            g.zeta[i]
        );
    }
}

#[test]
fn too_coarse_grid_is_rejected() {
    let hull = wigley(10.0, 1.0, 0.625);
    let cond = Conditions::seawater(3.0);
    let members = [(&hull, Placement::default())];
    let mut spec = spectrum(&members, &cond);
    // λ_t ≈ 5.77 m; 3 m pixels cannot carry it.
    assert!(spec.elevation_grid(-40.0, -10.0, 0.0, 0.0, 11, 1).is_err());
}

#[test]
fn amplitude_is_zero_outside_domain() {
    let hull = wigley(10.0, 1.0, 0.625);
    let cond = Conditions::seawater(3.0);
    let members = [(&hull, Placement::default())];
    let mut spec = spectrum(&members, &cond);
    for theta in [std::f64::consts::FRAC_PI_2, 2.0, -2.0, f64::NAN] {
        let a = spec.amplitude(theta);
        assert!(a.re == 0.0 && a.im == 0.0, "A({theta}) = {a:?}");
    }
}

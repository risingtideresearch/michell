//! End-to-end validation against closed-form results for the Wigley hull,
//! for which Michell's inner integrals have an analytic expression.

use michell::{hulls, BSplineSurface, Conditions, Hull, Placement, WaveOptions};

const G: f64 = michell::STANDARD_GRAVITY;

fn rel_err(got: f64, want: f64) -> f64 {
    (got - want).abs() / want.abs().max(f64::MIN_POSITIVE)
}

/// Analytic Michell inner integrals for the Wigley hull
/// f = (B/2)(1-(2x/L)²)(1-(z/T)²), x ∈ [-L/2, L/2], z ∈ [0, T]:
/// I = 0 (fore-aft symmetry about the midpoint), and
/// J(λ) = -(4B/L²) · X(k) · Z(κ), k = νλ, κ = νλ²,
/// X = ∫ x sin(kx) dx = 2 (sin(ka)/k² - a cos(ka)/k), a = L/2,
/// Z = ∫ (1-(z/T)²) e^{-κz} dz = (1-E)/κ - (2 - E(κ²T²+2κT+2))/(κ³T²).
fn wigley_j(l: f64, b: f64, t: f64, nu: f64, lambda: f64) -> f64 {
    let a = l / 2.0;
    let k = nu * lambda;
    let kappa = nu * lambda * lambda;
    let x_int = 2.0 * ((k * a).sin() / (k * k) - a * (k * a).cos() / k);
    let e = (-kappa * t).exp();
    let z_int = (1.0 - e) / kappa
        - (2.0 - e * (kappa * kappa * t * t + 2.0 * kappa * t + 2.0)) / (kappa.powi(3) * t * t);
    -(4.0 * b / (l * l)) * x_int * z_int
}

#[test]
fn wigley_inner_integrals_match_analytic() {
    let (l, b, t) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(l, b, t).unwrap();
    for fn_ in [0.2, 0.3, 0.45] {
        let u = fn_ * (G * l).sqrt();
        let cond = Conditions::freshwater(u);
        let nu = G / (u * u);
        for lambda in [1.0, 1.05, 1.3, 2.0, 3.7, 8.0, 20.0] {
            let (i, j) = michell::inner_integrals(&hull, &cond, lambda).unwrap();
            let want = wigley_j(l, b, t, nu, lambda);
            let scale = want.abs().max(1e-12);
            assert!(
                i.abs() < 1e-10 * scale.max(1.0),
                "Fn={fn_} λ={lambda}: I = {i}, expected 0"
            );
            assert!(
                (j - want).abs() < 1e-9 * scale,
                "Fn={fn_} λ={lambda}: J = {j}, want {want}"
            );
        }
    }
}

/// Reference outer integral by brute-force Simpson in θ using the analytic J,
/// with a multihull interference factor applied to the integrand.
#[allow(clippy::too_many_arguments)]
fn reference_wave_resistance_factored(
    l: f64,
    b: f64,
    t: f64,
    u: f64,
    rho: f64,
    lambda_max: f64,
    n: usize,
    factor: impl Fn(f64) -> f64,
) -> f64 {
    let nu = G / (u * u);
    let theta_max = (1.0 / lambda_max).acos();
    let dt = theta_max / n as f64;
    let f = |theta: f64| -> f64 {
        let sec = 1.0 / theta.cos();
        let j = wigley_j(l, b, t, nu, sec);
        j * j * factor(theta) * sec * sec * sec
    };
    let mut s = f(0.0) + f(theta_max);
    for i in 1..n {
        let w = if i % 2 == 1 { 4.0 } else { 2.0 };
        s += w * f(i as f64 * dt);
    }
    let integral = s * dt / 3.0;
    4.0 * rho * G * G / (std::f64::consts::PI * u * u) * integral
}

fn reference_wave_resistance(l: f64, b: f64, t: f64, u: f64, rho: f64) -> f64 {
    reference_wave_resistance_factored(l, b, t, u, rho, 200.0, 2_000_000, |_| 1.0)
}

/// Doctors & Beck, "Numerical Aspects of the Neumann-Kelvin Problem",
/// Journal of Ship Research 31(1), 1987, Table 1, reports `10^3 Cw = 1.2486`
/// for a Wigley hull at Fn=0.35 with B/L=0.1 and T/L=0.0625.
/// DOI: https://doi.org/10.5957/jsr.1987.31.1.1
///
/// The value and proportions are primary-source verified. This reproduction
/// uses the library's existing wetted-area coefficient convention; it does
/// not claim that the source text independently establishes an identical
/// normalization convention.
#[test]
fn wigley_matches_published_thin_ship_value() {
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(length, beam, draft).unwrap();
    let speed = 0.35 * (G * length).sqrt();
    let cond = Conditions::freshwater(speed);
    let wave = michell::wave_resistance_with(
        &hull,
        &cond,
        &WaveOptions {
            rel_tol: 1e-8,
            max_refinements: 7,
        },
    )
    .unwrap();
    let cw = wave.resistance / (0.5 * cond.fluid.density * speed * speed * hull.wetted_surface());
    let published_cw = 1.2486e-3;

    // This 0.5% reproduction tolerance covers convention and physical-
    // constant differences; it is not presented as the paper's uncertainty.
    assert!(
        rel_err(cw, published_cw) < 5e-3,
        "Cw={cw:.10e}, published={published_cw:.10e}"
    );
}

/// Geometrically similar hulls at equal length Froude number have invariant
/// wave-resistance coefficient and wave resistance proportional to L^3.
#[test]
fn froude_similarity_scales_wave_resistance_cubically() {
    let base = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let fn_ = 0.32;
    let base_speed = fn_ * (G * base.length()).sqrt();
    let base_cond = Conditions::freshwater(base_speed);
    let base_wave = michell::wave_resistance(&base, &base_cond)
        .unwrap()
        .resistance;

    for scale in [0.25_f64, 0.5, 2.0, 4.0] {
        let hull = hulls::wigley(10.0 * scale, scale, 0.625 * scale).unwrap();
        let speed = fn_ * (G * hull.length()).sqrt();
        let cond = Conditions::freshwater(speed);
        let wave = michell::wave_resistance(&hull, &cond).unwrap().resistance;
        assert!(
            rel_err(wave, base_wave * scale.powi(3)) < 2e-8,
            "scale={scale}: Rw={wave}, expected {}",
            base_wave * scale.powi(3)
        );
    }
}

/// Requested accuracy is checked against an analytic Wigley inner amplitude
/// integrated by a separate dense Simpson outer quadrature. This checks the
/// achieved bound for these cases; it does not claim every tolerance setting
/// forces an additional refinement.
#[test]
fn requested_accuracy_matches_analytic_reference() {
    let (l, b, t) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(l, b, t).unwrap();
    for fn_ in [0.20, 0.35, 0.50] {
        let u = fn_ * (G * l).sqrt();
        let cond = Conditions::freshwater(u);
        let reference = reference_wave_resistance(l, b, t, u, cond.fluid.density);
        let mut previous_actual = f64::INFINITY;
        for rel_tol in [1e-3, 1e-5] {
            let result = michell::wave_resistance_with(
                &hull,
                &cond,
                &WaveOptions {
                    rel_tol,
                    max_refinements: 7,
                },
            )
            .unwrap();
            let actual = rel_err(result.resistance, reference);
            assert!(
                actual <= previous_actual * 1.05,
                "Fn={fn_} tol={rel_tol}: actual error grew to {actual:.3e}"
            );
            assert!(
                actual < 5.0 * rel_tol,
                "Fn={fn_} tol={rel_tol}: actual error {actual:.3e}"
            );
            previous_actual = actual;
        }
    }
}

/// Internal consistency check for the production outer-quadrature error
/// estimator. Both calculations share the same moments and quadrature code,
/// so this does not serve as an independent physics oracle.
#[test]
fn reported_outer_error_is_internally_consistent() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    for fn_ in [0.08, 0.20, 0.35, 0.50] {
        let speed = fn_ * (G * hull.length()).sqrt();
        let cond = Conditions::freshwater(speed);
        let reference = michell::wave_resistance_with(
            &hull,
            &cond,
            &WaveOptions {
                // Six forced refinements remain below the per-pass safety cap
                // at Fn=0.08; finer passes can truncate the reference earlier.
                rel_tol: 1e-14,
                max_refinements: 6,
            },
        )
        .unwrap();
        let result = michell::wave_resistance_with(
            &hull,
            &cond,
            &WaveOptions {
                rel_tol: 1e-5,
                max_refinements: 7,
            },
        )
        .unwrap();
        assert!(reference.max_lambda >= 0.95 * result.max_lambda);
        let actual = rel_err(result.resistance, reference.resistance);
        assert!(
            actual <= 10.0 * result.est_rel_error.max(1e-10),
            "Fn={fn_}: estimated {:.3e}, internal difference {actual:.3e}",
            result.est_rel_error
        );
    }
}

#[test]
fn wigley_wave_resistance_matches_reference() {
    let (l, b, t) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(l, b, t).unwrap();
    let opts = WaveOptions {
        rel_tol: 1e-7,
        max_refinements: 6,
    };
    for fn_ in [0.25, 0.35] {
        let u = fn_ * (G * l).sqrt();
        let cond = Conditions::freshwater(u);
        let got = michell::wave_resistance_with(&hull, &cond, &opts).unwrap();
        let want = reference_wave_resistance(l, b, t, u, cond.fluid.density);
        let rel = (got.resistance - want).abs() / want;
        assert!(
            rel < 2e-5,
            "Fn={fn_}: got {} N, reference {} N (rel {rel:.2e}, est {:.2e})",
            got.resistance,
            want,
            got.est_rel_error
        );
        assert!(got.resistance > 0.0);
    }
}

#[test]
fn coincident_pair_quadruples_wave_resistance() {
    // Two coincident hulls: amplitudes add, |2F|² = 4|F|².
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let cond = Conditions::seawater(3.0);
    let single = michell::wave_resistance(&hull, &cond).unwrap().resistance;
    let pair = [(&hull, Placement::default()), (&hull, Placement::default())];
    let both = michell::multihull_wave_resistance(&pair, &cond)
        .unwrap()
        .resistance;
    assert!(
        (both - 4.0 * single).abs() < 1e-9 * both,
        "pair {both} vs 4x single {}",
        4.0 * single
    );
}

#[test]
fn catamaran_matches_analytic_interference() {
    // Two Wigley demihulls at y = ±s/2: interference factor
    // 4 cos²(½ ν s sec θ tan θ).
    let (l, b, t) = (10.0, 1.0, 0.625);
    let s = 2.0f64;
    let hull = hulls::wigley(l, b, t).unwrap();
    let fn_ = 0.35;
    let u = fn_ * (G * l).sqrt();
    let cond = Conditions::freshwater(u);
    let nu = G / (u * u);
    let members = [
        (&hull, Placement { x: 0.0, y: s / 2.0 }),
        (
            &hull,
            Placement {
                x: 0.0,
                y: -s / 2.0,
            },
        ),
    ];
    let opts = WaveOptions {
        rel_tol: 1e-7,
        max_refinements: 6,
    };
    let got = michell::multihull_wave_resistance_with(&members, &cond, &opts)
        .unwrap()
        .resistance;
    let want = reference_wave_resistance_factored(
        l,
        b,
        t,
        u,
        cond.fluid.density,
        100.0,
        4_000_000,
        |theta| {
            let sec = 1.0 / theta.cos();
            let alpha = 0.5 * nu * s * sec * theta.tan();
            4.0 * alpha.cos().powi(2)
        },
    );
    let rel = (got - want).abs() / want;
    assert!(
        rel < 2e-4,
        "catamaran: got {got} N, reference {want} N (rel {rel:.2e})"
    );
}

#[test]
fn tandem_matches_analytic_interference() {
    // Two Wigley hulls in line, staggered by d: interference factor
    // 4 cos²(½ ν d sec θ) from the longitudinal phase.
    let (l, b, t) = (10.0, 1.0, 0.625);
    let d = 6.0f64;
    let hull = hulls::wigley(l, b, t).unwrap();
    let fn_ = 0.35;
    let u = fn_ * (G * l).sqrt();
    let cond = Conditions::freshwater(u);
    let nu = G / (u * u);
    let members = [
        (&hull, Placement { x: d / 2.0, y: 0.0 }),
        (
            &hull,
            Placement {
                x: -d / 2.0,
                y: 0.0,
            },
        ),
    ];
    let opts = WaveOptions {
        rel_tol: 1e-7,
        max_refinements: 6,
    };
    let got = michell::multihull_wave_resistance_with(&members, &cond, &opts)
        .unwrap()
        .resistance;
    let want = reference_wave_resistance_factored(
        l,
        b,
        t,
        u,
        cond.fluid.density,
        100.0,
        4_000_000,
        |theta| {
            let alpha = 0.5 * nu * d / theta.cos();
            4.0 * alpha.cos().powi(2)
        },
    );
    let rel = (got - want).abs() / want;
    assert!(
        rel < 2e-4,
        "tandem: got {got} N, reference {want} N (rel {rel:.2e})"
    );
}

#[test]
fn multihull_breakdown_is_consistent() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let cond = Conditions::seawater(3.0);
    let members = [
        (&hull, Placement { x: 0.0, y: 1.5 }),
        (&hull, Placement { x: 0.0, y: -1.5 }),
    ];
    let r = michell::multihull_resistance(&members, &cond).unwrap();
    assert_eq!(r.viscous.len(), 2);
    let solo = michell::wave_resistance(&hull, &cond).unwrap().resistance;
    assert!((r.solo_wave_total - 2.0 * solo).abs() < 1e-9 * r.solo_wave_total);
    assert!((r.interference - r.wave.resistance / r.solo_wave_total).abs() < 1e-12);
    assert!((r.total - r.wave.resistance - r.viscous_total).abs() < 1e-9 * r.total);
    assert!((r.effective_power - r.total * 3.0).abs() < 1e-9 * r.effective_power);
    assert!((r.wetted_surface - 2.0 * hull.wetted_surface()).abs() < 1e-9 * r.wetted_surface);
}

#[test]
fn wigley_geometry_integrals() {
    let (l, b, t) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(l, b, t).unwrap();
    // Volume: ∇ = 2 (B/2)(2L/3)(2T/3) = 4 B L T / 9.
    let vol_exact = 4.0 * b * l * t / 9.0;
    assert!(
        (hull.displaced_volume() - vol_exact).abs() < 1e-12 * vol_exact,
        "volume {} vs {vol_exact}",
        hull.displaced_volume()
    );
    // Wetted surface: 2D Simpson reference with analytic derivatives.
    let n = 800usize;
    let (dx, dz) = (l / n as f64, t / n as f64);
    let integrand = |x: f64, z: f64| -> f64 {
        let fx = -4.0 * b / (l * l) * x * (1.0 - (z / t).powi(2));
        let fz = -(b / 2.0) * (1.0 - (2.0 * x / l).powi(2)) * 2.0 * z / (t * t);
        (1.0 + fx * fx + fz * fz).sqrt()
    };
    let mut s = 0.0;
    for i in 0..=n {
        let wx = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let x = -l / 2.0 + i as f64 * dx;
        for j in 0..=n {
            let wz = if j == 0 || j == n {
                1.0
            } else if j % 2 == 1 {
                4.0
            } else {
                2.0
            };
            s += wx * wz * integrand(x, j as f64 * dz);
        }
    }
    let wetted_ref = 2.0 * s * dx * dz / 9.0;
    let rel = (hull.wetted_surface() - wetted_ref).abs() / wetted_ref;
    assert!(
        rel < 1e-8,
        "wetted {} vs {wetted_ref} (rel {rel:.2e})",
        hull.wetted_surface()
    );
    assert!((hull.length() - l).abs() < 1e-12);
    assert!((hull.draft() - t).abs() < 1e-12);
    // Hydrostatics closed forms for the Wigley hull (symmetric about x = 0):
    // LCB = LCF = 0, A_w = 2BL/3, second moment = BL³/30.
    assert!(hull.lcb_x().abs() < 1e-10, "lcb {}", hull.lcb_x());
    let aw_exact = 2.0 * b * l / 3.0;
    assert!(
        (hull.waterplane_area() - aw_exact).abs() < 1e-12 * aw_exact,
        "Aw {} vs {aw_exact}",
        hull.waterplane_area()
    );
    assert!(hull.lcf_x().abs() < 1e-10);
    let ixx_exact = b * l.powi(3) / 30.0;
    assert!(
        (hull.waterplane_second_moment() - ixx_exact).abs() < 1e-12 * ixx_exact,
        "Ixx {} vs {ixx_exact}",
        hull.waterplane_second_moment()
    );
}

#[test]
fn multispan_inner_integrals_match_brute_force() {
    // Multi-span cubic × quadratic hull; compare the closed-form inner
    // integrals against direct 2D Gauss-Legendre quadrature of
    // ∬ fx e^{-κz} (cos, sin)(k(x - x_mid)) dx dz at low λ (mild oscillation,
    // so a dense per-span rule is a trustworthy reference).
    let knots_x = vec![0.0, 0.0, 0.0, 0.0, 2.5, 5.0, 7.5, 10.0, 10.0, 10.0, 10.0];
    let knots_z = vec![0.0, 0.0, 0.0, 0.6, 1.2, 1.2, 1.2];
    let (nx, nz) = (7usize, 4usize);
    let mut control = vec![0.0; nx * nz];
    for i in 0..nx {
        for j in 0..nz {
            // Positive, closing at the ends in x; smooth-ish in z.
            let gi = [0.0, 0.35, 0.8, 1.0, 0.8, 0.35, 0.0][i];
            control[i * nz + j] = gi * (1.5 - 0.31 * j as f64);
        }
    }
    let surf = BSplineSurface::new(3, 2, knots_x, knots_z, control).unwrap();
    let hull = Hull::new(surf).unwrap();

    let u = 6.0;
    let cond = Conditions::seawater(u);
    let nu = G / (u * u);
    let x_mid = 5.0;

    // Dense Gauss-Legendre per span (nodes/weights on [-1,1], order 30).
    let order = 30usize;
    let (gn, gw) = gauss_legendre_ref(order);
    let xs = [0.0, 2.5, 5.0, 7.5, 10.0];
    let zs = [0.0, 0.6, 1.2];

    for lambda in [1.0, 1.2, 1.6] {
        let k = nu * lambda;
        let kappa = nu * lambda * lambda;
        let (mut i_ref, mut j_ref) = (0.0f64, 0.0f64);
        for w in xs.windows(2) {
            let (x0, x1) = (w[0], w[1]);
            for v in zs.windows(2) {
                let (z0, z1) = (v[0], v[1]);
                let jac = (x1 - x0) / 2.0 * ((z1 - z0) / 2.0);
                for (a, &na) in gn.iter().enumerate() {
                    let x = x0 + (x1 - x0) * (na + 1.0) / 2.0;
                    for (b, &nb) in gn.iter().enumerate() {
                        let z = z0 + (z1 - z0) * (nb + 1.0) / 2.0;
                        let fx = hull.surface().eval_deriv(x, z, 1, 0);
                        let common = gw[a] * gw[b] * jac * fx * (-kappa * z).exp();
                        i_ref += common * (k * (x - x_mid)).cos();
                        j_ref += common * (k * (x - x_mid)).sin();
                    }
                }
            }
        }
        let (i_got, j_got) = michell::inner_integrals(&hull, &cond, lambda).unwrap();
        let scale = (i_ref * i_ref + j_ref * j_ref).sqrt().max(1e-12);
        assert!(
            (i_got - i_ref).abs() < 1e-10 * scale && (j_got - j_ref).abs() < 1e-10 * scale,
            "λ={lambda}: got ({i_got}, {j_got}), want ({i_ref}, {j_ref})"
        );
    }
}

/// Local Gauss-Legendre reference (independent of the crate's internal one).
fn gauss_legendre_ref(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..200 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        nodes[i] = -x;
        weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    (nodes, weights)
}

#[test]
fn wedge_prism_geometry() {
    // Piecewise-linear tent in x, constant in z: a wall-sided "wedge" prism.
    // Volume = B L T / 2; wetted = 2 L T sqrt(1 + (B/L)²).
    let (l, b, t) = (8.0, 1.2, 0.9);
    let knots_x = vec![0.0, 0.0, l / 2.0, l, l];
    let knots_z = vec![0.0, 0.0, t, t];
    let control = vec![0.0, 0.0, b / 2.0, b / 2.0, 0.0, 0.0];
    let hull = Hull::new(BSplineSurface::new(1, 1, knots_x, knots_z, control).unwrap()).unwrap();
    let vol = b * l * t / 2.0;
    let wet = 2.0 * l * t * (1.0 + (b / l) * (b / l)).sqrt();
    assert!((hull.displaced_volume() - vol).abs() < 1e-12 * vol);
    assert!((hull.wetted_surface() - wet).abs() < 1e-10 * wet);
}

#[test]
fn froude_sweep_is_finite_positive_and_humped() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let mut cws = Vec::new();
    for i in 0..=20 {
        let fn_ = 0.15 + 0.35 * i as f64 / 20.0;
        let u = fn_ * (G * 10.0f64).sqrt();
        let r = michell::resistance(&hull, &Conditions::seawater(u)).unwrap();
        assert!(
            r.wave.resistance.is_finite() && r.wave.resistance >= 0.0,
            "Fn={fn_}: Rw = {}",
            r.wave.resistance
        );
        assert!(r.viscous.resistance > 0.0 && r.total > r.wave.resistance);
        cws.push(r.cw);
    }
    // The Cw(Fn) curve must not be monotone (humps and hollows).
    let increasing = cws.windows(2).all(|w| w[1] >= w[0]);
    let decreasing = cws.windows(2).all(|w| w[1] <= w[0]);
    assert!(!increasing && !decreasing, "Cw curve unexpectedly monotone");
}

#[test]
fn rejects_bad_hulls_and_conditions() {
    // Negative half-beam.
    let s = BSplineSurface::new(
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.5, 0.5, -0.5, 0.5],
    )
    .unwrap();
    assert!(Hull::new(s).is_err());
    // z-domain not starting at the waterline.
    let s = BSplineSurface::new(
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.5, 0.5, 1.0, 1.0],
        vec![0.5; 4],
    )
    .unwrap();
    assert!(Hull::new(s).is_err());
    // Nonsense speed.
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    assert!(michell::wave_resistance(&hull, &Conditions::seawater(0.0)).is_err());
    assert!(michell::wave_resistance(&hull, &Conditions::seawater(f64::NAN)).is_err());
}

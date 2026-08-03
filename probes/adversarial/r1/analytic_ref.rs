use michell::{wave_resistance_with, BSplineSurface, Conditions, Hull, WaveOptions, STANDARD_GRAVITY};
use std::f64::consts::PI;

fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let m = n.div_ceil(2);
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..m {
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        loop {
            let mut p0 = 1.0;
            let mut p1 = z;
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            let derivative = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            let next = z - p1 / derivative;
            if (next - z).abs() < 2.0e-16 {
                z = next;
                let weight = 2.0 / ((1.0 - z * z) * derivative * derivative);
                x[i] = -z;
                x[n - 1 - i] = z;
                w[i] = weight;
                w[n - 1 - i] = weight;
                break;
            }
            z = next;
        }
    }
    (x, w)
}

fn exact_inner_squared(lambda: f64, p: usize, nu: f64, length: f64, beam: f64, draft: f64) -> f64 {
    let h = 0.5 * length;
    let k = nu * lambda;
    let a = nu * lambda * lambda * draft;
    let decay = (-a).exp();
    let i0 = (1.0 - decay) / a;
    let i1 = (1.0 - decay * (1.0 + a)) / (a * a);
    let i2 = (2.0 - decay * (a * a + 2.0 * a + 2.0)) / (a * a * a);
    let z_transform = draft * (i0 - 0.5 * i1 - 0.5 * i2);
    let x_transform = 2.0 * ((k * h).sin() / (k * k) - h * (k * h).cos() / k);
    let amplitude = (-4.0 * beam * (1.0 - 1.0 / p as f64) / (length * length))
        * x_transform
        * z_transform;
    amplitude * amplitude
}

fn independent_reference(p: usize, fn_: f64, lambda_max: f64) -> (f64, usize) {
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let speed = fn_ * (STANDARD_GRAVITY * length).sqrt();
    let nu = STANDARD_GRAVITY / (speed * speed);
    let omega = nu * length;
    let du = PI / omega;
    let panels = ((lambda_max - 1.0) / du).ceil() as usize;
    let (nodes, weights) = gauss_legendre(32);
    let mut integral = 0.0;
    let mut correction = 0.0;
    for panel in 0..panels {
        let u0 = panel as f64 * du;
        let u1 = ((panel + 1) as f64 * du).min(lambda_max - 1.0);
        let t0 = u0.sqrt();
        let t1 = u1.sqrt();
        let half = 0.5 * (t1 - t0);
        let mid = 0.5 * (t1 + t0);
        for (&node, &weight) in nodes.iter().zip(&weights) {
            let t = mid + half * node;
            let lambda = 1.0 + t * t;
            let term = weight * half * 2.0 * lambda * lambda / (2.0 + t * t).sqrt()
                * exact_inner_squared(lambda, p, nu, length, beam, draft);
            let y = term - correction;
            let next = integral + y;
            correction = (next - integral) - y;
            integral = next;
        }
    }
    let coeff = 4.0 * 999.1 * STANDARD_GRAVITY * STANDARD_GRAVITY / (PI * speed * speed);
    (coeff * integral, panels * nodes.len())
}

fn hull(p: usize) -> Hull {
    let q = 2;
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let mut control = Vec::new();
    for i in 0..=p {
        let u = i as f64 / p as f64;
        let bx = 2.0 * beam * u * (1.0 - u);
        for j in 0..=q {
            let v = j as f64 / q as f64;
            control.push(bx * (1.0 - v * v));
        }
    }
    let knots_x = std::iter::repeat_n(-0.5 * length, p + 1)
        .chain(std::iter::repeat_n(0.5 * length, p + 1))
        .collect();
    let knots_z = std::iter::repeat_n(0.0, q + 1)
        .chain(std::iter::repeat_n(draft, q + 1))
        .collect();
    Hull::new(BSplineSurface::new(p, q, knots_x, knots_z, control).unwrap()).unwrap()
}

fn main() {
    let fn_ = 0.05;
    for p in [2, 48] {
        let h = hull(p);
        let speed = fn_ * (STANDARD_GRAVITY * h.length()).sqrt();
        let production = wave_resistance_with(&h, &Conditions::freshwater(speed), &WaveOptions::default()).unwrap();
        for cutoff in [1_000.0, 2_000.0, 4_000.0] {
            let (reference, evaluations) = independent_reference(p, fn_, cutoff);
            println!(
                "p={p} cutoff={cutoff:.0} independent={reference:.12e} evals={evaluations} production={:.12e} method={:?} est={:.3e} actual_rel={:.3e}",
                production.resistance,
                production.method,
                production.est_rel_error,
                (production.resistance - reference).abs() / reference,
            );
        }
    }
}

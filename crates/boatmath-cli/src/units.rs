//! Quantities on the command line with their units: `16in`, `8kn`, `1kN`.
//! A bare number is SI (m, m/s, N).

/// The longest leading number, and the unit after it.
fn split(s: &str) -> (f64, &str) {
    let t = s.trim();
    (1..=t.len())
        .rev()
        .filter(|&i| t.is_char_boundary(i))
        .find_map(|i| {
            t[..i]
                .trim()
                .parse::<f64>()
                .ok()
                .map(|v| (v, t[i..].trim()))
        })
        .unwrap_or((f64::NAN, t))
}

fn with(s: &str, what: &str, units: &[(&str, f64)]) -> Result<f64, String> {
    let (v, u) = split(s);
    let k = units
        .iter()
        .find(|(name, _)| name.eq_ignore_ascii_case(u))
        .map(|(_, k)| *k)
        .ok_or_else(|| {
            let names: Vec<&str> = units
                .iter()
                .map(|u| u.0)
                .filter(|u| !u.is_empty())
                .collect();
            format!("{s:?}: a {what} in {}", names.join(", "))
        })?;
    if !v.is_finite() {
        return Err(format!("{s:?}: expected a number"));
    }
    Ok(v * k)
}

/// A length [m]: m (default), mm, cm, in, ft.
pub fn length(s: &str) -> Result<f64, String> {
    with(
        s,
        "length",
        &[
            ("", 1.0),
            ("m", 1.0),
            ("mm", 1e-3),
            ("cm", 1e-2),
            ("in", 0.0254),
            ("ft", 0.3048),
        ],
    )
}

/// A speed [m/s]: m/s (default), kn, mph, km/h.
pub fn speed(s: &str) -> Result<f64, String> {
    with(
        s,
        "speed",
        &[
            ("", 1.0),
            ("m/s", 1.0),
            ("kn", 0.514444),
            ("kt", 0.514444),
            ("mph", 0.44704),
            ("km/h", 1.0 / 3.6),
        ],
    )
}

/// A force [N]: N (default), kN, kgf, lbf.
pub fn force(s: &str) -> Result<f64, String> {
    with(
        s,
        "force",
        &[
            ("", 1.0),
            ("N", 1.0),
            ("kN", 1e3),
            ("kgf", 9.81),
            ("lbf", 4.448222),
        ],
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantities() {
        assert!((length("16in").unwrap() - 0.4064).abs() < 1e-12);
        assert!((length("400mm").unwrap() - 0.4).abs() < 1e-12);
        assert!((length("0.4").unwrap() - 0.4).abs() < 1e-12);
        assert!((speed("8kn").unwrap() - 4.115552).abs() < 1e-9);
        assert!((force("1kN").unwrap() - 1000.0).abs() < 1e-12);
        assert!((force("1e3 N").unwrap() - 1000.0).abs() < 1e-12);
        assert!(length("16 furlongs").is_err());
        assert!(force("kN").is_err());
    }
}

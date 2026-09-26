//! Argument parsing helpers: lengths, roughness specs, units, ranges.

use michell::Roughness;

pub fn parse_length(s: &str) -> Result<f64, String> {
    let t = s.trim();
    let (num, scale) = [
        ("um", 1e-6),
        ("µm", 1e-6),
        ("mm", 1e-3),
        ("cm", 1e-2),
        ("ft", 0.3048),
        ("in", 0.0254),
        ("m", 1.0),
    ]
    .iter()
    .find_map(|(suf, k)| t.strip_suffix(suf).map(|n| (n, *k)))
    .unwrap_or((t, 1.0));
    num.trim()
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && *v >= 0.0)
        .map(|v| v * scale)
        .ok_or_else(|| format!("cannot parse length {t:?} (try 100um, 0.1mm, or 1e-4)"))
}

/// Parse a `--roughness` spec: `off`, `cf=VALUE` (a ΔC_F added outside the
/// form factor), or `ks=LENGTH` (equivalent sand-grain height).
pub fn parse_roughness(spec: &str) -> Result<Roughness, String> {
    let t = spec.trim();
    if matches!(t, "off" | "none" | "smooth") {
        return Ok(Roughness::None);
    }
    if let Some(v) = t.strip_prefix("cf=") {
        let c = v
            .trim()
            .parse::<f64>()
            .ok()
            .filter(|x| x.is_finite() && *x >= 0.0)
            .ok_or_else(|| format!("--roughness cf={v:?}: expected a non-negative number"))?;
        return Ok(Roughness::DeltaCf(c));
    }
    if let Some(v) = t.strip_prefix("ks=") {
        return Ok(Roughness::SandGrain(parse_length(v)?));
    }
    Err(format!(
        "--roughness {t:?}: expected off | cf=DELTA_CF | ks=HEIGHT (e.g. ks=100um)"
    ))
}

/// Parse a units name (or raw scale) to metres-per-unit.
pub fn parse_units(s: &str) -> Result<f64, String> {
    match s.trim() {
        "mm" => Ok(0.001),
        "cm" => Ok(0.01),
        "m" => Ok(1.0),
        "in" => Ok(0.0254),
        "ft" => Ok(0.3048),
        other => other
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && *v > 0.0)
            .ok_or_else(|| format!("--units {other:?}: expected mm|cm|m|in|ft or a scale")),
    }
}

/// Binary-STL detection: exact size match on the triangle count.
pub fn looks_binary_stl(bytes: &[u8]) -> bool {
    bytes.len() >= 84 && {
        let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        bytes.len() == 84 + 50 * n
    }
}

/// Parse `a`, or an inclusive range `a:b:step`.
pub fn parse_range(s: &str) -> Result<Vec<f64>, String> {
    let parts: Vec<&str> = s.split(':').collect();
    let num = |t: &str| {
        t.parse::<f64>()
            .map_err(|_| format!("cannot parse number {t:?} in {s:?}"))
    };
    match parts.as_slice() {
        [one] => Ok(vec![num(one)?]),
        [a, b, step] => {
            let (a, b, step) = (num(a)?, num(b)?, num(step)?);
            if !(step > 0.0 && b >= a) {
                return Err(format!("bad range {s:?}: need start <= end and step > 0"));
            }
            let n = ((b - a) / step + 1e-9).floor() as usize;
            Ok((0..=n).map(|i| a + i as f64 * step).collect())
        }
        _ => Err(format!(
            "bad range {s:?}: expected `value` or `start:end:step`"
        )),
    }
}

/// Parse `NxM` (e.g. `12x8`).
pub fn parse_pair(s: &str) -> Result<(usize, usize), String> {
    let (a, b) = s
        .split_once(['x', 'X'])
        .ok_or_else(|| format!("expected NxM, got {s:?}"))?;
    Ok((
        a.parse().map_err(|_| format!("bad integer {a:?}"))?,
        b.parse().map_err(|_| format!("bad integer {b:?}"))?,
    ))
}

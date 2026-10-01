//! Lists of numbers on the command line, `a,b,c` or `start:stop:step`, and
//! every combination of several of them.

/// Parse `a,b,c`, `start:stop:step` (stop included when the steps land on
/// it), or a mixture: `0.2:0.4:0.1,0.5`.
pub fn parse(s: &str) -> Result<Vec<f64>, String> {
    let num = |t: &str| {
        t.trim()
            .parse::<f64>()
            .ok()
            .filter(|x| x.is_finite())
            .ok_or_else(|| format!("{t:?}: expected a number"))
    };
    let mut out = Vec::new();
    for part in s.split(',').filter(|p| !p.trim().is_empty()) {
        let bits: Vec<&str> = part.split(':').collect();
        match bits[..] {
            [x] => out.push(num(x)?),
            [a, b, step] => {
                let (a, b, step) = (num(a)?, num(b)?, num(step)?);
                if step <= 0.0 || b < a {
                    return Err(format!(
                        "{part}: expected start:stop:step, start ≤ stop, step > 0"
                    ));
                }
                let n = ((b - a) / step + 1e-9).floor() as usize;
                if n > 10_000 {
                    return Err(format!("{part}: more than 10000 values"));
                }
                // Rounded to the step's precision, so 0.1 steps print as 0.3,
                // not 0.30000000000000004.
                let digits = (-step.log10().floor()).max(0.0) as i32 + 6;
                let scale = 10f64.powi(digits);
                out.extend((0..=n).map(|i| ((a + i as f64 * step) * scale).round() / scale));
            }
            _ => return Err(format!("{part}: expected a number or start:stop:step")),
        }
    }
    if out.is_empty() {
        return Err(format!("{s:?}: an empty list"));
    }
    Ok(out)
}

/// Every combination of the given lists, one value from each, the last
/// varying fastest; an axis left out (`None`) is `None` in each.
pub fn product(axes: &[Option<Vec<f64>>]) -> Vec<Vec<Option<f64>>> {
    let mut out: Vec<Vec<Option<f64>>> = vec![Vec::new()];
    for axis in axes {
        out = out
            .into_iter()
            .flat_map(|prefix| {
                let values: Vec<Option<f64>> = match axis {
                    Some(v) => v.iter().copied().map(Some).collect(),
                    None => vec![None],
                };
                values.into_iter().map(move |v| {
                    let mut p = prefix.clone();
                    p.push(v);
                    p
                })
            })
            .collect();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lists_and_ranges() {
        assert_eq!(parse("1,2.5").unwrap(), vec![1.0, 2.5]);
        assert_eq!(parse("0.2:0.5:0.1").unwrap(), vec![0.2, 0.3, 0.4, 0.5]);
        assert_eq!(parse("0:1:0.4,3").unwrap(), vec![0.0, 0.4, 0.8, 3.0]);
        assert!(parse("1:0:0.1").is_err());
        assert!(parse("x").is_err());
        assert!(parse("").is_err());
    }

    #[test]
    fn every_combination() {
        let p = product(&[Some(vec![1.0, 2.0]), None, Some(vec![3.0, 4.0])]);
        assert_eq!(p.len(), 4);
        assert_eq!(p[1], vec![Some(1.0), None, Some(4.0)]);
    }
}

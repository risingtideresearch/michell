//! Line plots as SVG: one series per `y` field and `--by` value, points
//! joined in order of x. The SVG carries its own light and dark colours, a
//! legend for two or more series, a label at each line's end for up to four,
//! and a tooltip on every point.

use std::fmt::Write;

/// The reference categorical palette, light and dark, in its fixed order.
const SERIES: [(&str, &str); 8] = [
    ("#2a78d6", "#3987e5"),
    ("#eb6834", "#d95926"),
    ("#1baf7a", "#199e70"),
    ("#eda100", "#c98500"),
    ("#e87ba4", "#d55181"),
    ("#008300", "#008300"),
    ("#4a3aa7", "#9085e9"),
    ("#e34948", "#e66767"),
];

pub struct Series {
    pub name: String,
    pub points: Vec<(f64, f64)>,
}

pub struct Plot {
    pub title: Option<String>,
    pub x_label: String,
    pub y_label: String,
    pub series: Vec<Series>,
}

/// About `n` round tick values covering `[lo, hi]`.
pub(crate) fn ticks(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    let span = (hi - lo).max(f64::MIN_POSITIVE);
    let raw = span / n as f64;
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 2.5, 5.0, 10.0]
        .iter()
        .map(|m| m * mag)
        .find(|s| *s >= raw * (1.0 - 1e-9))
        .unwrap_or(10.0 * mag);
    let first = (lo / step - 1e-9).ceil() as i64;
    let last = (hi / step + 1e-9).floor() as i64;
    (first..=last).map(|i| i as f64 * step).collect()
}

/// A tick label: no more digits than the step needs.
pub(crate) fn label(v: f64, step: f64) -> String {
    let digits = (-step.log10().floor()).max(0.0) as usize;
    let s = format!("{v:.digits$}");
    if s == "-0" {
        "0".into()
    } else {
        s
    }
}

/// A value for a tooltip: four significant figures.
pub(crate) fn sig(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return format!("{v}");
    }
    let digits = (3 - v.abs().log10().floor() as i32).max(0) as usize;
    let s = format!("{v:.digits$}");
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    }
}

pub(crate) fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Pad a range that is a single value, so it can be drawn.
fn range(vs: impl Iterator<Item = f64>) -> (f64, f64) {
    let (lo, hi) = vs.fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| {
        (a.min(v), b.max(v))
    });
    if lo == hi {
        let d = if lo == 0.0 { 1.0 } else { lo.abs() * 0.1 };
        (lo - d, hi + d)
    } else {
        (lo, hi)
    }
}

pub fn svg(p: &Plot) -> Result<String, String> {
    if p.series.len() > SERIES.len() {
        return Err(format!(
            "{} series; at most {} — narrow --by, or plot them in several charts",
            p.series.len(),
            SERIES.len()
        ));
    }
    let all = || p.series.iter().flat_map(|s| s.points.iter());
    if all().next().is_none() {
        return Err("nothing to plot: no record has numbers at those fields".into());
    }
    let n = p.series.len();
    let legend = n >= 2;
    let direct = n <= 4 && legend;
    let (w, h) = (720.0, 440.0);
    let top = if p.title.is_some() { 40.0 } else { 16.0 } + if legend { 28.0 } else { 0.0 };
    let (left, bottom) = (64.0, 48.0);
    let longest = p
        .series
        .iter()
        .map(|s| s.name.chars().count())
        .max()
        .unwrap_or(0);
    let right = if direct {
        (16.0 + 7.0 * longest as f64).clamp(24.0, 240.0)
    } else {
        24.0
    };
    let (pw, ph) = (w - left - right, h - top - bottom);

    let (x0, x1) = range(all().map(|q| q.0));
    let (y0, y1) = range(all().map(|q| q.1));
    let xt = ticks(x0, x1, 6);
    let yt = ticks(y0, y1, 5);
    // Widen the y range to whole ticks, so the plot starts and ends on a line.
    let ystep = yt.get(1).zip(yt.first()).map_or(y1 - y0, |(b, a)| b - a);
    let xstep = xt.get(1).zip(xt.first()).map_or(x1 - x0, |(b, a)| b - a);
    let (y0, y1) = (
        (y0 / ystep + 1e-9).floor() * ystep,
        (y1 / ystep - 1e-9).ceil() * ystep,
    );
    let yt = ticks(y0, y1, 5);
    let sx = |x: f64| left + (x - x0) / (x1 - x0) * pw;
    let sy = |y: f64| top + ph - (y - y0) / (y1 - y0) * ph;

    let mut s = String::new();
    let mut css =
        String::from(":root{--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--grid:#e4e3df;");
    let mut dark = String::from("--surface:#1a1a19;--ink:#ffffff;--ink2:#c3c2b7;--grid:#3a3a37;");
    for (i, (l, d)) in SERIES.iter().enumerate().take(n) {
        let _ = write!(css, "--s{i}:{l};");
        let _ = write!(dark, "--s{i}:{d};");
    }
    css.push('}');
    let _ = write!(
        s,
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" font-family="system-ui, -apple-system, sans-serif" font-size="12">
<style>{css}@media (prefers-color-scheme: dark){{:root{{{dark}}}}}
text{{fill:var(--ink2)}} .title{{fill:var(--ink);font-size:15px;font-weight:600}}
.pt:hover{{stroke-width:6}}</style>
<rect width="{w}" height="{h}" fill="var(--surface)"/>
"#
    );
    if let Some(t) = &p.title {
        let _ = writeln!(
            s,
            r#"<text class="title" x="{left}" y="24">{}</text>"#,
            esc(t)
        );
    }
    // Grid and axes, recessive.
    for &y in &yt {
        let _ = writeln!(
            s,
            r#"<line x1="{left}" x2="{}" y1="{y:.1}" y2="{y:.1}" stroke="var(--grid)"/><text x="{}" y="{:.1}" text-anchor="end">{}</text>"#,
            left + pw,
            left - 8.0,
            sy(y) + 4.0,
            label(y, ystep),
            y = sy(y)
        );
    }
    for &x in &xt {
        let _ = writeln!(
            s,
            r#"<line x1="{x:.1}" x2="{x:.1}" y1="{}" y2="{}" stroke="var(--grid)"/><text x="{x:.1}" y="{}" text-anchor="middle">{}</text>"#,
            top + ph,
            top + ph + 5.0,
            top + ph + 20.0,
            label(x, xstep),
            x = sx(x)
        );
    }
    let _ = writeln!(
        s,
        r#"<line x1="{left}" x2="{}" y1="{b}" y2="{b}" stroke="var(--ink2)"/>"#,
        left + pw,
        b = top + ph
    );
    let _ = writeln!(
        s,
        r#"<text x="{}" y="{}" text-anchor="middle">{}</text>"#,
        left + pw / 2.0,
        h - 10.0,
        esc(&p.x_label)
    );
    let _ = writeln!(
        s,
        r#"<text transform="translate(16 {}) rotate(-90)" text-anchor="middle">{}</text>"#,
        top + ph / 2.0,
        esc(&p.y_label)
    );
    // Legend, in one row above the plot.
    if legend {
        let mut x = left;
        let y = top - 20.0;
        for (i, ser) in p.series.iter().enumerate() {
            let _ = writeln!(
                s,
                r#"<rect x="{x}" y="{}" width="12" height="3" rx="1.5" fill="var(--s{i})"/><text x="{}" y="{}">{}</text>"#,
                y - 4.0,
                x + 16.0,
                y,
                esc(&ser.name)
            );
            x += 28.0 + 7.0 * ser.name.chars().count() as f64;
        }
    }
    for (i, ser) in p.series.iter().enumerate() {
        let mut pts = ser.points.clone();
        pts.sort_by(|a, b| a.0.total_cmp(&b.0));
        let line: Vec<String> = pts
            .iter()
            .map(|&(x, y)| format!("{:.1},{:.1}", sx(x), sy(y)))
            .collect();
        let _ = writeln!(
            s,
            r#"<polyline points="{}" fill="none" stroke="var(--s{i})" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>"#,
            line.join(" ")
        );
        for &(x, y) in &pts {
            let _ = writeln!(
                s,
                r#"<circle class="pt" cx="{:.1}" cy="{:.1}" r="4" fill="var(--s{i})" stroke="var(--surface)" stroke-width="2"><title>{}{} = {}, {} = {}</title></circle>"#,
                sx(x),
                sy(y),
                if n > 1 {
                    format!("{}: ", esc(&ser.name))
                } else {
                    String::new()
                },
                esc(&p.x_label),
                sig(x),
                esc(&p.y_label),
                sig(y)
            );
        }
    }
    // A label at each line's end, pushed apart so none overlap.
    if direct {
        let mut ends: Vec<(f64, f64, &str)> = p
            .series
            .iter()
            .filter_map(|ser| {
                let &(x, y) = ser.points.iter().max_by(|a, b| a.0.total_cmp(&b.0))?;
                Some((sx(x) + 8.0, sy(y) + 4.0, ser.name.as_str()))
            })
            .collect();
        ends.sort_by(|a, b| a.1.total_cmp(&b.1));
        for i in 1..ends.len() {
            ends[i].1 = ends[i].1.max(ends[i - 1].1 + 14.0);
        }
        for (x, y, name) in ends {
            let _ = writeln!(
                s,
                r#"<text x="{x:.1}" y="{y:.1}" style="fill:var(--ink)">{}</text>"#,
                esc(name)
            );
        }
    }
    s.push_str("</svg>\n");
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_ticks() {
        assert_eq!(
            ticks(0.0, 1.0, 5),
            vec![0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0]
        );
        assert_eq!(label(0.6000000000000001, 0.2), "0.6");
        assert_eq!(ticks(0.25, 0.55, 6).len(), 7);
    }
}

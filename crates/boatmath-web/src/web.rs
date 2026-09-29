//! The pages and their shared scripts, built into the binary. With
//! `BOATMATH_WEB_DIR` set (to `crates/boatmath-web/src/web`) they are read
//! from there on every request instead, so a page can be edited without a
//! rebuild.

/// `(name, contents)` of every file under `src/web/`.
const FILES: &[(&str, &str)] = &[
    ("common.css", include_str!("web/common.css")),
    ("common.js", include_str!("web/common.js")),
    ("viewer.js", include_str!("web/viewer.js")),
    ("chart.js", include_str!("web/chart.js")),
    ("hulls.html", include_str!("web/hulls.html")),
    ("hull-new.html", include_str!("web/hull-new.html")),
    ("hull.html", include_str!("web/hull.html")),
    ("cases.html", include_str!("web/cases.html")),
    ("case-new.html", include_str!("web/case-new.html")),
    ("case.html", include_str!("web/case.html")),
    ("studies.html", include_str!("web/studies.html")),
    ("studies-new.html", include_str!("web/studies-new.html")),
    ("study.html", include_str!("web/study.html")),
    ("queue.html", include_str!("web/queue.html")),
    ("plot.html", include_str!("web/plot.html")),
];

/// The page a path shows, if any.
fn page(path: &str) -> Option<&'static str> {
    let seg: Vec<&str> = path.trim_matches('/').split('/').collect();
    let is_id = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    Some(match seg[..] {
        ["hulls"] => "hulls.html",
        ["hulls", "new"] => "hull-new.html",
        ["hulls", id] if is_id(id) => "hull.html",
        ["cases"] => "cases.html",
        ["cases", "new"] => "case-new.html",
        ["cases", id] if is_id(id) => "case.html",
        ["studies"] => "studies.html",
        ["studies", "new"] => "studies-new.html",
        ["studies", id] if is_id(id) => "study.html",
        ["queue"] => "queue.html",
        // Results was the plot's first name.
        ["plot"] | ["results"] => "plot.html",
        _ => return None,
    })
}

/// The file a GET for `path` is answered with, and its content type:
/// a page, or `/static/<file>`.
pub fn file(path: &str) -> Option<(String, &'static str)> {
    let name = match path.strip_prefix("/static/") {
        Some(n) => n,
        None => page(path)?,
    };
    let kind = match name.rsplit_once('.').map(|(_, e)| e) {
        Some("html") => "text/html; charset=utf-8",
        Some("js") => "text/javascript; charset=utf-8",
        Some("css") => "text/css; charset=utf-8",
        _ => return None,
    };
    let built = FILES.iter().find(|(n, _)| *n == name)?.1;
    if let Ok(dir) = std::env::var("BOATMATH_WEB_DIR") {
        if let Ok(s) = std::fs::read_to_string(std::path::Path::new(&dir).join(name)) {
            return Some((s, kind));
        }
    }
    Some((built.to_string(), kind))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paths_find_their_pages() {
        assert_eq!(page("/hulls"), Some("hulls.html"));
        assert_eq!(page("/hulls/"), Some("hulls.html"));
        assert_eq!(page("/hulls/new"), Some("hull-new.html"));
        assert_eq!(page("/hulls/12"), Some("hull.html"));
        assert_eq!(page("/hulls/x1"), None);
        assert_eq!(page("/studies/new"), Some("studies-new.html"));
        assert_eq!(page("/cases/4"), Some("case.html"));
        assert_eq!(page("/studies/7"), Some("study.html"));
        assert_eq!(page("/results"), Some("plot.html"));
        assert_eq!(page("/plot"), Some("plot.html"));
        assert_eq!(page("/cases"), Some("cases.html"));
        assert_eq!(page("/cases/new"), Some("case-new.html"));
        assert_eq!(page("/studies"), Some("studies.html"));
        assert_eq!(page("/queue"), Some("queue.html"));
        assert!(file("/static/viewer.js").is_some());
        assert!(file("/static/../main.rs").is_none());
        assert!(file("/static/nope.js").is_none());
    }
}

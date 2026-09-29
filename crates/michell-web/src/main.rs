//! `michell-web [--port N] [--host ADDR] [--data DIR]` — serve the web
//! front end, and work through its queue of runs.
//!
//! Binds `127.0.0.1:8080` by default. When `$PORT` is set (as on a hosting
//! platform) it listens there on all interfaces instead. Bound to loopback
//! it is meant to sit behind `tailscale serve`, and takes who is asking
//! from the `Tailscale-User-Login` and `Tailscale-User-Name` headers that
//! sets (and on any other address ignores them: anyone could send them
//! there). Hulls,
//! configurations, runs and results are kept in `DIR` (default
//! `$MICHELL_DATA`, else `./michell-data`); see `api.rs` for the endpoints
//! and `web.rs` for the pages. Two stateless endpoints serve the upload
//! page's preview:
//!
//!   POST /api/loft?name=F&...      body = the file's bytes; the cut, as
//!                                  JSON (see `michell_web::loft`)
//!   POST /api/geometry?name=F&units=U
//!                                  body = the file's bytes; its display
//!                                  geometry (see `michell_web::geometry`)

use michell_web::{geometry, loft, LoftRequest, MAX_UPLOAD};
use std::io::Read;
use std::sync::Arc;
use tiny_http::{Header, Method, Request, Response, Server};

mod api;
mod web;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let flag = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let env_port = std::env::var("PORT").ok();
    let port = flag("--port").or(env_port.clone()).unwrap_or("8080".into());
    let host = flag("--host").unwrap_or(if env_port.is_some() {
        "0.0.0.0".into()
    } else {
        "127.0.0.1".into()
    });
    let data = flag("--data")
        .or_else(|| std::env::var("MICHELL_DATA").ok())
        .unwrap_or("michell-data".into());
    let store = match michell_web::store::Store::open(std::path::Path::new(&data)) {
        Ok(s) => Arc::new(s),
        Err(e) => {
            eprintln!("michell-web: cannot open the store in {data}: {e}");
            std::process::exit(1);
        }
    };
    let worker = michell_web::worker::Worker::new(Arc::clone(&store));
    worker.spawn();
    let app = Arc::new(api::App { store, worker });
    eprintln!(
        "michell-web: data in {data}, solver {}",
        michell_web::store::SOLVER_VERSION
    );
    let addr = format!("{host}:{port}");
    let server = match Server::http(&addr) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("michell-web: cannot listen on {addr}: {e}");
            std::process::exit(1);
        }
    };
    // Only a proxy on this machine can reach a loopback listener, so only
    // then are its identity headers believable.
    let trust = matches!(host.as_str(), "127.0.0.1" | "localhost" | "::1" | "[::1]");
    eprintln!(
        "michell-web: serving http://{addr}/{}",
        if trust {
            " (identity from tailscale serve's headers)"
        } else {
            ""
        }
    );
    for req in server.incoming_requests() {
        // Cuts and statics are CPU-bound and fan out across cores
        // themselves; a thread per request keeps a slow IGES import from
        // blocking the pages.
        let app = Arc::clone(&app);
        std::thread::spawn(move || handle(req, &app, trust));
    }
}

/// Who `tailscale serve` says is asking, if it says.
fn tailnet_user(req: &Request) -> Option<api::Who> {
    let get = |k: &str| {
        req.headers()
            .iter()
            .find(|h| h.field.as_str().as_str().eq_ignore_ascii_case(k))
            .map(|h| h.value.as_str().trim().to_string())
            .filter(|v| !v.is_empty())
    };
    let login = get("Tailscale-User-Login")?;
    Some(api::Who {
        name: get("Tailscale-User-Name").unwrap_or_default(),
        login,
    })
}

fn handle(mut req: Request, app: &api::App, trust: bool) {
    let url = req.url().to_string();
    let (path, query) = url.split_once('?').unwrap_or((&url, ""));
    let pairs = parse_query(query);
    let method = req.method().clone();
    let who = if trust { tailnet_user(&req) } else { None };
    let routed = api::route(app, &method, path, &pairs, who.as_ref(), &mut || {
        read_body(&mut req)
    });
    if let Some(r) = routed {
        let resp = match r {
            Ok(api::Reply::Json(v)) => json_response(200, v.to_string()),
            Ok(api::Reply::Gz(gz, kind)) => gz_response(&req, gz, kind),
            Err((code, e)) => {
                if code >= 500 {
                    eprintln!("{method} {path}: {e}");
                }
                json_response(code, serde_json::json!({ "error": e }).to_string())
            }
        };
        let _ = req.respond(resp);
        return;
    }
    if method == Method::Get {
        if let Some((body, kind)) = web::file(path) {
            let _ = req.respond(
                Response::from_string(body)
                    .with_header(header("Content-Type", kind))
                    .with_header(header("Cache-Control", "no-cache")),
            );
            return;
        }
    }
    let name = pairs
        .iter()
        .find(|(k, _)| k == "name")
        .map_or("upload", |(_, v)| v.as_str())
        .to_string();
    let answer = |r: Result<serde_json::Value, String>| match r {
        Ok(v) => json_response(200, v.to_string()),
        Err(e) => json_response(400, serde_json::json!({ "error": e }).to_string()),
    };
    let resp = match (&method, path) {
        (Method::Get, "/") => Response::from_string("")
            .with_status_code(302)
            .with_header(header("Location", "/hulls")),
        (Method::Post, "/api/geometry") => {
            let units = pairs
                .iter()
                .find(|(k, v)| k == "units" && !v.trim().is_empty())
                .map(|(_, v)| michell_cli::parse_units(v.trim()))
                .transpose();
            answer(units.and_then(|u| geometry(&name, read_body(&mut req)?, u)))
        }
        (Method::Post, "/api/loft") => answer(
            LoftRequest::from_query(&pairs)
                .and_then(|opts| loft(&name, read_body(&mut req)?, &opts)),
        ),
        _ => Response::from_string("not found").with_status_code(404),
    };
    let _ = req.respond(resp);
}

fn read_body(req: &mut Request) -> Result<Vec<u8>, String> {
    if req.body_length().is_some_and(|n| n > MAX_UPLOAD) {
        return Err(format!("file too large (limit {} MB)", MAX_UPLOAD >> 20));
    }
    let mut bytes = Vec::new();
    req.as_reader()
        .take(MAX_UPLOAD as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| format!("upload failed: {e}"))?;
    if bytes.len() > MAX_UPLOAD {
        return Err(format!("file too large (limit {} MB)", MAX_UPLOAD >> 20));
    }
    Ok(bytes)
}

fn json_response(status: u16, body: String) -> Response<std::io::Cursor<Vec<u8>>> {
    Response::from_string(body)
        .with_status_code(status)
        .with_header(header("Content-Type", "application/json"))
}

/// A stored (gzipped) blob, sent as it is to a client that takes gzip.
fn gz_response(req: &Request, gz: Vec<u8>, kind: &str) -> Response<std::io::Cursor<Vec<u8>>> {
    let takes_gzip = req
        .headers()
        .iter()
        .any(|h| h.field.equiv("Accept-Encoding") && h.value.as_str().contains("gzip"));
    if takes_gzip {
        return Response::from_data(gz)
            .with_header(header("Content-Type", kind))
            .with_header(header("Content-Encoding", "gzip"));
    }
    let mut out = Vec::new();
    match flate2::read::GzDecoder::new(&gz[..]).read_to_end(&mut out) {
        Ok(_) => Response::from_data(out).with_header(header("Content-Type", kind)),
        Err(e) => json_response(
            500,
            serde_json::json!({ "error": e.to_string() }).to_string(),
        ),
    }
}

fn header(k: &str, v: &str) -> Header {
    Header::from_bytes(k.as_bytes(), v.as_bytes()).expect("static header")
}

fn parse_query(q: &str) -> Vec<(String, String)> {
    q.split('&')
        .filter(|p| !p.is_empty())
        .map(|p| {
            let (k, v) = p.split_once('=').unwrap_or((p, ""));
            (percent_decode(k), percent_decode(v))
        })
        .collect()
}

fn percent_decode(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'+' => out.push(b' '),
            b'%' if i + 2 < b.len() => {
                match std::str::from_utf8(&b[i + 1..i + 3])
                    .ok()
                    .and_then(|h| u8::from_str_radix(h, 16).ok())
                {
                    Some(v) => {
                        out.push(v);
                        i += 2;
                    }
                    None => out.push(b'%'),
                }
            }
            c => out.push(c),
        }
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

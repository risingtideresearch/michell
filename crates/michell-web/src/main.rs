//! `michell-web [--port N] [--host ADDR]` — serve the web front end.
//!
//! Binds `127.0.0.1:8080` by default. When `$PORT` is set (as on a hosting
//! platform) it listens there on all interfaces instead.
//!
//! Endpoints:
//!   GET  /                      the page
//!   POST /api/loft?name=F&...   body = the file's bytes; returns the loft as
//!                               JSON (see `michell_web::loft`), or
//!                               `{"error": ...}` with status 400
//!   POST /api/flow?name=F&froude=Fn&closure=..&param=..&...
//!                               body = the file's bytes; the near-field
//!                               pressure, free surface and forces at that
//!                               speed (see `michell_web::flow`)

use michell_web::{flow, loft, FlowRequest, LoftRequest, MAX_UPLOAD};
use std::io::Read;
use tiny_http::{Header, Method, Request, Response, Server};

const PAGE: &str = include_str!("index.html");

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
    let addr = format!("{host}:{port}");
    let server = match Server::http(&addr) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("michell-web: cannot listen on {addr}: {e}");
            std::process::exit(1);
        }
    };
    eprintln!("michell-web: serving http://{addr}/");
    for req in server.incoming_requests() {
        // Lofts are CPU-bound and fan out across cores themselves; a thread
        // per request keeps a slow IGES import from blocking the page.
        std::thread::spawn(move || handle(req));
    }
}

fn handle(mut req: Request) {
    let url = req.url().to_string();
    let (path, query) = url.split_once('?').unwrap_or((&url, ""));
    let pairs = parse_query(query);
    let resp = match (req.method(), path) {
        (Method::Get, "/") => Response::from_string(PAGE)
            .with_header(header("Content-Type", "text/html; charset=utf-8")),
        (Method::Post, "/api/loft") => {
            let name = pairs
                .iter()
                .find(|(k, _)| k == "name")
                .map_or("upload", |(_, v)| v.as_str())
                .to_string();
            let result = read_body(&mut req)
                .and_then(|bytes| {
                    let opts = LoftRequest::from_query(&pairs)?;
                    loft(&name, bytes, &opts)
                })
                .map_err(|e| {
                    eprintln!("loft {name}: {e}");
                    e
                });
            match result {
                Ok(v) => {
                    eprintln!(
                        "loft {name}: ok ({:.2} s)",
                        v["seconds"].as_f64().unwrap_or(0.0)
                    );
                    json_response(200, v.to_string())
                }
                Err(e) => json_response(400, serde_json::json!({ "error": e }).to_string()),
            }
        }
        (Method::Post, "/api/flow") => {
            let name = pairs
                .iter()
                .find(|(k, _)| k == "name")
                .map_or("upload", |(_, v)| v.as_str())
                .to_string();
            let result = read_body(&mut req).and_then(|bytes| {
                let opts = FlowRequest::from_query(&pairs)?;
                flow(&name, bytes, &opts)
            });
            match result {
                Ok(v) => {
                    eprintln!(
                        "flow {name}: ok ({:.2} s)",
                        v["seconds"].as_f64().unwrap_or(0.0)
                    );
                    json_response(200, v.to_string())
                }
                Err(e) => {
                    eprintln!("flow {name}: {e}");
                    json_response(400, serde_json::json!({ "error": e }).to_string())
                }
            }
        }
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

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
//!                               speed (see `michell_web::flow`). With
//!                               `&progress=1` the answer streams as
//!                               newline-delimited JSON: `{"progress": ..}`
//!                               lines as the stages advance, then the
//!                               result (or `{"error": ..}`) as the last
//!                               line; a client that goes away stops it

use michell_web::{
    flow, flow_with_progress, loft, span_sweep_with_progress, FlowRequest, LoftRequest,
    CANCELLED, MAX_UPLOAD,
};
use std::io::{Read, Write};
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
    if *req.method() == Method::Post
        && path == "/api/flow"
        && pairs.iter().any(|(k, v)| k == "progress" && v != "0")
    {
        return stream_flow(req, &pairs);
    }
    if *req.method() == Method::Post && path == "/api/span_sweep" {
        return stream_sweep(req, &pairs);
    }
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

/// `/api/flow?progress=1`: the flow with its progress streamed as it runs.
///
/// The response is written straight to the connection, chunked, flushing
/// after every line — tiny_http's own chunked writer buffers ~8 kB, which
/// would hold the small progress lines back until the end. When a write
/// fails (the page aborted the request, or went away), the progress
/// callback says so and the flow stops at its next report.
fn stream_flow(mut req: Request, pairs: &[(String, String)]) {
    let name = pairs
        .iter()
        .find(|(k, _)| k == "name")
        .map_or("upload", |(_, v)| v.as_str())
        .to_string();
    let body = read_body(&mut req);
    let mut out = req.into_writer();
    let head = "HTTP/1.1 200 OK\r\nContent-Type: application/x-ndjson\r\n\
                Cache-Control: no-cache\r\nTransfer-Encoding: chunked\r\n\
                Connection: close\r\n\r\n";
    let chunk = |out: &mut Box<dyn Write + Send>, line: &str| -> std::io::Result<()> {
        let data = format!("{line}\n");
        write!(out, "{:x}\r\n", data.len())?;
        out.write_all(data.as_bytes())?;
        out.write_all(b"\r\n")?;
        out.flush()
    };
    if out.write_all(head.as_bytes()).and_then(|_| out.flush()).is_err() {
        return;
    }
    let t0 = std::time::Instant::now();
    let mut gone = false;
    let result = body.and_then(|bytes| {
        let opts = FlowRequest::from_query(pairs)?;
        flow_with_progress(&name, bytes, &opts, &mut |p| {
            let line = serde_json::json!({ "progress": {
                "stage": p.stage,
                "step": p.step,
                "steps": p.steps,
                "detail": p.detail,
                "fraction": p.fraction,
                "seconds": t0.elapsed().as_secs_f64(),
            }});
            // A write to a connection the client has closed succeeds once
            // (the kernel buffers it) and fails once the peer's reset is
            // back: lead with a blank keep-alive line, so the report after
            // an abort — not the one after that — finds out.
            gone = chunk(&mut out, "").is_err()
                || {
                    std::thread::sleep(std::time::Duration::from_millis(5));
                    chunk(&mut out, &line.to_string()).is_err()
                };
            !gone
        })
    });
    let last = match result {
        Ok(v) => {
            eprintln!(
                "flow {name}: ok ({:.2} s, streamed)",
                v["seconds"].as_f64().unwrap_or(0.0)
            );
            v.to_string()
        }
        Err(e) if e == CANCELLED || gone => {
            eprintln!(
                "flow {name}: cancelled after {:.1} s (client gone)",
                t0.elapsed().as_secs_f64()
            );
            return;
        }
        Err(e) => {
            eprintln!("flow {name}: {e}");
            serde_json::json!({ "error": e }).to_string()
        }
    };
    let _ = chunk(&mut out, &last).and_then(|_| {
        out.write_all(b"0\r\n\r\n")?;
        out.flush()
    });
}

/// `/api/span_sweep?spans=S1,S2,…&…` (the flow's other keys, without
/// `span`): a fast catamaran span sweep (see
/// `michell_web::span_sweep_with_progress`), streamed like `/api/flow`:
/// progress lines, a `{"span_index": i, "result": …}` line per span as it is
/// done, and a closing `{"done": true}` (or `{"error": …}`).
fn stream_sweep(mut req: Request, pairs: &[(String, String)]) {
    use std::cell::RefCell;
    let name = pairs
        .iter()
        .find(|(k, _)| k == "name")
        .map_or("upload", |(_, v)| v.as_str())
        .to_string();
    let body = read_body(&mut req);
    let out = RefCell::new(req.into_writer());
    let head = "HTTP/1.1 200 OK\r\nContent-Type: application/x-ndjson\r\n\
                Cache-Control: no-cache\r\nTransfer-Encoding: chunked\r\n\
                Connection: close\r\n\r\n";
    let chunk = |line: &str| -> std::io::Result<()> {
        let mut out = out.borrow_mut();
        let data = format!("{line}\n");
        write!(out, "{:x}\r\n", data.len())?;
        out.write_all(data.as_bytes())?;
        out.write_all(b"\r\n")?;
        out.flush()
    };
    {
        let mut o = out.borrow_mut();
        if o.write_all(head.as_bytes()).and_then(|_| o.flush()).is_err() {
            return;
        }
    }
    let t0 = std::time::Instant::now();
    let gone = std::cell::Cell::new(false);
    let result = body.and_then(|bytes| {
        let spans: Vec<f64> = pairs
            .iter()
            .find(|(k, _)| k == "spans")
            .map(|(_, v)| v.split(',').filter_map(|t| t.trim().parse().ok()).collect())
            .unwrap_or_default();
        let opts = FlowRequest::from_query(pairs)?;
        span_sweep_with_progress(
            &name,
            bytes,
            &opts,
            &spans,
            &mut |i, v| {
                let ok = chunk(&serde_json::json!({ "span_index": i, "result": v }).to_string()).is_ok();
                gone.set(gone.get() || !ok);
                ok
            },
            &mut |p| {
                let line = serde_json::json!({ "progress": {
                    "stage": p.stage,
                    "step": p.step,
                    "steps": p.steps,
                    "detail": p.detail,
                    "fraction": p.fraction,
                    "seconds": t0.elapsed().as_secs_f64(),
                }});
                let bad = chunk("").is_err() || {
                    std::thread::sleep(std::time::Duration::from_millis(5));
                    chunk(&line.to_string()).is_err()
                };
                gone.set(gone.get() || bad);
                !bad
            },
        )
    });
    let last = match result {
        Ok(()) => {
            eprintln!("span sweep {name}: ok ({:.2} s, streamed)", t0.elapsed().as_secs_f64());
            serde_json::json!({ "done": true }).to_string()
        }
        Err(e) if e == CANCELLED || gone.get() => {
            eprintln!(
                "span sweep {name}: cancelled after {:.1} s (client gone)",
                t0.elapsed().as_secs_f64()
            );
            return;
        }
        Err(e) => {
            eprintln!("span sweep {name}: {e}");
            serde_json::json!({ "error": e }).to_string()
        }
    };
    let _ = chunk(&last).and_then(|_| {
        let mut out = out.borrow_mut();
        out.write_all(b"0\r\n\r\n")?;
        out.flush()
    });
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

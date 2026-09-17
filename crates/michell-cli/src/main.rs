//! `michell` — thin-ship wave resistance from hull files.
//!
//! This binary is a thin wrapper: all of the command logic lives in the
//! `michell_cli` library (so the egui editor can loft and run sweeps
//! in-process, without shelling out to this binary). Progress is written to
//! stderr, exactly as before.

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    // The library reports progress as human-readable lines (optionally carrying
    // a machine-readable fraction); the CLI just echoes the line to stderr.
    let mut report = |line: &str, _frac: Option<(usize, usize)>| eprintln!("{line}");
    if let Err(e) = michell_cli::run(&args, &mut report) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

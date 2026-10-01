# boatmath CLI — design draft

A command-line boatmath: the web app's hulls → cases → studies, as JSON
records that small commands make, read and pass along a pipe.

```sh
boatmath hull guillemot.igs --waterline 0.12 > g.hull.json
boatmath case --mass 150 --vcg 0.25 < g.hull.json > g150.case.json
boatmath study --froude 0.2:0.6:0.05 < g150.case.json | boatmath run > calm.jsonl
boatmath plot -x froude -y forces.rt < calm.jsonl > rt.svg
```

## Principles

- **One record, one JSON object; streams are JSONL.** Every command reads
  records on stdin (or files named as arguments) and writes them on stdout.
  Diagnostics and progress go to stderr.
- **Every record carries `type` and `id`.** `id` is the SHA-256 of the
  record's canonical inputs (every default filled in, as `params.rs` does now),
  so asking for the same thing twice gives the same id.
- **Records refer to their parents by id.** A case names its hull
  (`"hull": "<id>"`), a study its case, a result its study. Every record a
  command writes is also saved in the store (below), so a downstream command
  can look its parents up. Field paths in `table`/`plot` follow ids, so
  `study.case.params.span` works on a result.
- **The CAD file is the geometry.** The solver cuts the hull afresh at every
  attitude, so a hull record names its file (path and SHA-256) and the import
  settings. It doesn't hold sections. `run` refuses a file whose hash has
  changed.
- **The store.** `$BOATMATH_HOME` (default `~/.boatmath`) holds every
  record by id (`records/ab/cdef….json`) and bulk data (free-surface fields,
  display meshes) as content-addressed gzipped blobs (`blobs/ab/cdef….gz`),
  the layout the web store uses. Records refer to blobs as
  `{"blob": "<sha256>"}`. A result's id is its study's, so re-running a
  pipeline finds the results already there. Records are plain files: copy a
  store to share it, delete it to start over. `boatmath get ID` prints one.
- **Generating and computing are separate.** `study` only writes requests,
  which is cheap and lets you edit them with jq. `run` does the work.
- **Plain records for plotting.** `table` flattens records to TSV/CSV for
  gnuplot, pandas or a spreadsheet. `plot` covers the common cases.

## Records

### hull

```json
{
  "type": "hull", "id": "…", "name": "guillemot",
  "source": { "path": "/abs/guillemot.igs", "sha256": "…", "kind": "iges" },
  "import": { "waterline": 0.12, "stations": 41, "rays": 64 },
  "parent": null,
  "summary": { "length": 4.6, "beam": 0.61, "draft": 0.14, "volume": 0.146, "lcb": 2.31 }
}
```

`import` is the existing `LoftRequest` (waterline, centerplane, stations,
rays, units, scale, scale_yz). A scaled hull is the same source with
`scale`/`scale_yz` set and `parent` the id of the hull it came from.

### case

```json
{
  "type": "case", "id": "…", "name": "",
  "hull": "<hull id>",
  "params": { "span": null, "mass": 150, "mass_by": "sinking", "lcg": null, "vcg": 0.25,
              "kxx": null, "kyy": null, "kzz": null, "roll_damping": 0.0 },
  "statics": {
    "mass": 150, "lcg": 2.3, "vcg": 0.25,
    "at_rest": { "sinkage": …, "trim_deg": … },
    "hydrostatics": { … }, "roll": { "gm_t": …, "period": … },
    "gz": { "heel_deg": [...], "gz": [...], "gm": …, "max_gz": …, "vanishing_deg": … }
  }
}
```

`params` is `CaseParams`. `statics` is what `platform::statics` returns now,
less `meshes`/`body_meshes`, which move to the cache.

### study

```json
{ "type": "study", "id": "…", "case": "<case id>",
  "params": { "froude": 0.35, "dynamic": true, "closure": { "type": "ballistic", "coeff": … },
              "grid": 640, "waves": null } }
```

`params` is `StudyParams`. A study in waves has `waves: {heading, lambdas, sea}`.

### result

```json
{ "type": "result", "id": "<study id>", "study": "<study id>",
  "solver_version": "…", "seconds": 4.2,
  "attitude": { "sinkage": …, "trim_deg": … },
  "forces": { "rw": …, "rv": …, "rt": …, "pe": …, "cw": …, "ct": …, … },
  "surface": { "blob": "…" } }
```

A result in waves has `seakeeping` (the RAOs at each λ/L and the sea
statistics) in place of `forces`/`surface`, plus the peak scalars that
`wave_scalars` computes today (`heave_peak`, `sea_accel_bow`, …).

## Commands

| command | in → out | does |
|---|---|---|
| `hull FILE… [--waterline --stations --rays --units --centerplane --name]` | files → hull | cut the file and summarise it |
| `scale [--by K] [--beam K]` | hull → hull | a derived hull |
| `case [--span --mass --mass-by --lcg --vcg --kxx --kyy --kzz --roll-damping --name]` | hull → case | load it and compute statics |
| `study [--froude LIST] [--hold] [--closure …] [--waves HEADINGS --lambdas LIST --sea jonswap:hs=1,tp=6]` | case → study* | requests for every combination |
| `run [-j N] [--force]` | study* → result* | compute, cached; a study in waves finds or computes its calm-water attitude first |
| `gz-waves --length L --height H` | case → gz | quasi-static GZ in a regular wave |
| `table [--explode PATH] FIELD…` | any* → TSV/CSV | dotted paths as columns; `--explode seakeeping.headings.0.points` gives a row per λ |
| `plot -x F -y F [--by F] [--explode PATH] [-o out.svg]` | any* → SVG/PNG | line plot, one series per `--by` value |
| `get ID…` | → record* | print stored records |
| `show` | hull/case/result → browser | open the record in the 3-D viewer |

`LIST` is `a,b,c` or `start:stop:step`. Any flag on `case`/`study` that takes
a list makes one record per combination, so a sweep is one command:

```sh
boatmath case --span 2.4:3.6:0.2 --mass 150,200 < g.hull.json \
  | boatmath study --froude 0.25:0.55:0.05 \
  | boatmath run -j 8 \
  | boatmath plot -x study.case.params.span -y forces.rt --by study.params.froude
```

Seakeeping, head seas, RAOs per λ:

```sh
boatmath study --froude 0.3 --waves 180 --sea jonswap:hs=0.5,tp=3 < g150.case.json \
  | boatmath run \
  | boatmath plot --explode seakeeping.headings.0.points -x lambda -y '|heave|' -y '|pitch|'
```

(`|…|` is the modulus of a `[re, im]` pair, the one bit of extra syntax in
field paths.)

## Code

- Platforms: a single hull, or a catamaran of two by `span`, as in the web
  app. Other multi-hull platforms come later.
- `plot` writes SVG itself, with no plotting dependency.
- A new library crate, **`boatmath`**, takes what the CLI and the web app
  share out of `boatmath-web`: `params.rs`, `LoftRequest`, `platform.rs`
  (statics, wave GZ, waves) and the calm-water flow. `boatmath-web` then
  depends on it, so both front ends compute the same way.
- A new binary crate, **`boatmath-cli`** (binary `boatmath`): the record
  types, the cache, the commands. It uses clap and serde_json, unlike the
  hand-rolled `michell-cli`. `michell-cli` stays for now; we can retire it
  once `boatmath` covers what you use.

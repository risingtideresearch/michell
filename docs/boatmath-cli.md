# The boatmath CLI

The web app's hulls → cases → studies on the command line, as JSON records
that small commands make, read and pass along a pipe.

```sh
boatmath hull guillemot.igs --waterline 0.12 > g.hull.json
boatmath case --mass 150 --vcg 0.25 < g.hull.json > g150.case.json
boatmath study --froude 0.2:0.6:0.05 < g150.case.json | boatmath run > calm.jsonl
boatmath plot -x study.params.froude -y forces.rt -o rt.svg < calm.jsonl
```

## Principles

- **One record, one JSON object; streams are JSONL.** Commands read records
  on stdin and write them on stdout, one per line. Progress and errors go to
  stderr. A command that fails on some records carries on with the rest and
  exits 1.
- **Every record carries `type` and `id`.** `id` is the SHA-256 of the
  record's canonical inputs (with every default filled in, as `params.rs`
  does), so asking for the same thing twice gives the same id.
- **Records refer to their parents by id.** A case names its hull
  (`"hull": "<id>"`), a study its case, a result its study. Every record a
  command writes is also saved in the store, so a downstream command can look
  its parents up.
- **A hull record holds its geometry.** The solver cuts the hull afresh at
  every attitude, so a hull record holds what it cuts from: the hull's
  B-spline patches (from IGES) or its triangles (from STL), inline as JSON,
  not sections. `hull` reads the file once, applies its import settings
  (waterline, units) and splits it into hulls. After that the file is not
  needed: `source` records where the geometry came from, for reference only.
- **The store.** `$BOATMATH_HOME` (default `~/.boatmath`):

  ```text
  records/<type>/ab/cdef….json    every record, by type and id
  blobs/ab/cdef….gz               results' fields (free surface, meshes),
                                  by SHA-256 of the uncompressed bytes
  ```

  Records refer to blobs as `{"blob": "<sha256>"}`. A result's id is its
  study's, so a re-run finds the results already there. Records and results
  are stamped with the solver version, and one from another version is
  computed again. Records are plain files: copy a store to share it, delete it
  to start over.
- **Generating and computing are separate.** `study` writes requests, which is
  cheap and lets you edit them with jq. `run` does the work.

## Records

### hull

```json
{ "type": "hull", "id": "…", "name": "e12",
  "source": { "path": "/abs/e12.igs", "sha256": "…", "waterline": -0.95, "units": null },
  "cut": { "stations": null, "rays": null, "centerplane": null },
  "parent": null,
  "geometry": { "kind": "nurbs", "hulls": [ { "patches": [ … ] } ] },
  "summary": { "hulls": [ { "length": 8.63, "beam": 0.84, "draft": 0.22, "displaced_volume": 0.707,
                            "wetted_surface": …, "lcb_x": …, "waterplane_area": 5.58,
                            "transom": true } ],
               "notes": [] },
  "solver_version": "…" }
```

`geometry` is in metres, x and y as in the file, z up, with the design
waterline at z = 0. It comes in two kinds:

```json
{ "kind": "nurbs",
  "hulls": [ { "patches": [ { "degree": [3, 3], "n_ctrl": [8, 6],
                              "knots_u": [ … ], "knots_v": [ … ],
                              "ctrl": [ [x, y, z], … ], "trim_uv": null } ] } ] }
{ "kind": "mesh",
  "hulls": [ { "vertices": [ [x, y, z], … ], "triangles": [ [0, 1, 2], … ] } ] }
```

Patches are polynomial (clamped, unweighted) B-splines. `ctrl` is row-major
with v fastest, and `trim_uv` is the `[u0, u1, v0, v1]` box of a bounded
surface. Each entry of `hulls` is one hull, so a file holding two demihulls
gives two. The solver re-poses this geometry exactly as it would the file: a
trim, sinkage or scale maps the control points or vertices. `boatmath::native`
reads and writes it, and it's checked on reading (knot counts, sizes, finite
numbers). For e12, the record is 258 KB against the 1.7 MB IGES, and its cases
and studies come out the same as from the file.

`cut` holds the settings of the cut (stations, rays, a centreplane override).
The id is that of `geometry` plus `cut`. A scaled hull (`scale`) has the old
geometry with its control points scaled about the design waterline, and
`parent` set to the id of the hull it came from.

### case

```json
{ "type": "case", "id": "…", "name": "", "hull": "<hull id>",
  "params": { "span": 1.6, "mass": null, "mass_by": "sinking", "lcg": null, "vcg": null,
              "kxx": null, "kyy": null, "kzz": null, "roll_damping": 0.0 },
  "statics": { "mass": …, "lcg": …, "vcg": …,
               "at_rest": { "sinkage": …, "trim_rad": …, "trim_deg": … },
               "hydrostatics": { … }, "roll": { "gm_t": …, "period": … },
               "gz": { "heel_deg": [ … ], "gz": [ … ], "gm": …, "max_gz": …, … } },
  "seconds": 2.5, "solver_version": "…" }
```

`params` is `CaseParams`, and `span` makes a catamaran. `statics` is what
`platform::statics` returns, minus the display meshes. A case whose statics
fail is still saved, with `error` instead of `statics`.

### study

```json
{ "type": "study", "id": "…", "case": "<case id>", "kind": "calm",
  "params": { "froude": 0.35, "dynamic": true, "closure": { "type": "ballistic", "coeff": 1.414 },
              "grid": 640, "waves": null } }
```

`params` is `StudyParams`. A study in waves has
`waves: { "heading": 180, "lambdas": [ … ], "sea": … }` and `kind: "waves"`.

### result

In calm water:

```json
{ "type": "result", "id": "<study id>", "study": "<study id>", "kind": "calm",
  "froude": 0.35, "speed": …, "transverse_wavelength": …, "seconds": 41.0,
  "forces": { "rw": …, "rv": …, "rt": …, "pe": …, "cw": …, "ct": …, "sinkage": …, "trim_deg": …, … },
  "field": { "blob": "…" }, "solver_version": "…" }
```

In waves:

```json
{ "type": "result", "id": "<study id>", "study": "<study id>", "kind": "waves",
  "froude": …, "speed": …, "attitude": { "sinkage": …, "trim_deg": … },
  "calm": "<the calm-water study it is held about>",
  "peaks": { "heave_peak": …, "heave_peak_lambda": …, "sea_accel_bow": …, … },
  "seakeeping": { "gm_t": …, "headings": [ { "heading": 180, "points": [ … ], "sea": … } ], … },
  "field": { "blob": "…" }, "solver_version": "…" }
```

Each of the `seakeeping.headings.0.points` holds `lambda`, `omega`,
`omega_e`, the complex RAOs `heave`, `pitch`, `sway`, `roll` and `yaw` as
`[re, im]`, and the added resistance `raw_gb`.

## Commands

| command | in → out | does |
|---|---|---|
| `hull FILE… [--waterline LIST --stations --rays --units --centerplane --name]` | files → hull* | cut each file and summarise it |
| `scale [--by LIST] [--beam LIST] [--name]` | hull* → hull* | hulls scaled about the design waterline |
| `case [--span --mass --lcg --vcg --kxx --kyy --kzz --roll-damping (LISTs)] [--mass-by --name --force]` | hull* → case* | load each hull and compute its statics |
| `study --froude LIST [--hold] [--closure C] [--grid N] [--waves LIST --lambdas LIST --sea S]` | case* → study* | requests for every combination |
| `run [-j N] [--force] [-q]` | study* → result* | compute, or find in the store |
| `table FIELD… [--explode PATH] [--csv] [--no-header]` | any* → TSV/CSV | a column per field path |
| `plot -x F -y F… [--by F…] [--explode PATH] [--title --xlabel --ylabel] [-o FILE]` | any* → SVG | line plot, a series per y field and `--by` value |
| `get ID…` | → record* | print stored records, by id or a prefix of one |
| `blob ID` | → bytes | print a stored blob (a result's `field`) |

`LIST` is `a,b,c` or `start:stop:step` (or a mixture). Every option given a
list makes one record per combination, so a sweep is one pipeline:

```sh
boatmath case --span 1.6:2.4:0.2 --mass 150,200 < g.hull.json \
  | boatmath study --froude 0.25:0.55:0.05 \
  | boatmath run -j 8 \
  | boatmath plot -x study.params.froude -y forces.rt --by study.case.params.span -o rt.svg
```

`--closure` is `ballistic[:COEFF]`, `fixed:LENGTH` or `off`. `--sea` is
`bretschneider:hs=H,tp=T` or `jonswap:hs=H,tp=T[,gamma=G]`.

`run` first computes the calm-water studies, both those asked for and those
the studies in waves are held about, then the studies in waves. It runs
`-j` studies at once (default: the number of cores) and writes each result
as soon as it's ready. The results of prerequisite calm-water studies are
saved but not printed.

### Field paths

`table` and `plot` name their columns by path:

- `forces.rt`, `seakeeping.headings.0.points`: keys and array indices.
- A parent's id is followed into its record, so on a result
  `study.case.params.span` reads the span of the result's case. The keys
  followed are `hull`, `case`, `study`, `parent` and `calm`.
- `|heave|`: the modulus of a complex `[re, im]` pair.
- `id8`: the record's id, shortened to eight characters.

With `--explode PATH`, each element of the array at `PATH` is a row. Fields
are read from the element first, then from the record:

```sh
boatmath study --froude 0.3 --waves 180 --sea jonswap:hs=0.5,tp=3 < g150.case.json \
  | boatmath run \
  | boatmath plot --explode seakeeping.headings.0.points -x lambda -y '|heave|' -y '|pitch|'
```

### plot

The SVG uses the reference categorical palette, with light and dark colours
that follow the viewer's colour scheme. It has a legend for two or more
series, a label at the end of each line for up to four, and a tooltip on
every point. It takes at most eight series; narrow `--by` beyond that.

## Code

- **`crates/boatmath`** holds what the CLI and the web app share: `params.rs`
  (`CaseParams`, `StudyParams`), `LoftRequest`, `platform.rs` (statics, wave
  GZ, waves), the calm-water flow and `SOLVER_VERSION`. `boatmath-web`
  re-exports it and keeps the SQLite store, the worker and the pages, so both
  front ends compute the same way.
- **`crates/boatmath-cli`** (binary `boatmath`): `store.rs`, `records.rs`
  (hull, case, study), `run.rs`, `path.rs`, `plot.rs`, `list.rs`.
- `crates/boatmath/src/native.rs`: geometry as JSON, read from a file or
  re-posed. `boatmath`'s computations take its bytes wherever they take a
  hull file's, so the web app could store it too.
- Platforms are single hulls and catamarans (`span`), as in the web app.

## Not yet

- `gz-waves`: quasi-static GZ in a regular wave (`platform::wave_gz`).
- `show`: open a hull, case or result in the 3-D viewer.
- Warm starts: the web worker starts each equilibrium from the nearest
  speed already solved. `run` starts every one from scratch.
- A study in waves is held about the calm-water study at the default grid.
  A calm study at another grid shares its attitude but is not looked for.

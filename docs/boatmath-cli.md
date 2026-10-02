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
fail is still saved, with `error` instead of `statics`. `sections` is the id
of its hulls cut at rest (see below).

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
  "field": { "blob": "…" }, "sections": "<sections id>", "solver_version": "…" }
```

`field` holds the free surface (`surface`: a grid `nx × ny` over
`[x0, x1] × [y0, y1]`, ζ as base64 little-endian f32) and each hull's
pressure (`hulls`: `x`, `depth`, `half_beam` and `cp` on stations × depths,
and its centreplane `y`). It keeps no meshes: a picture re-poses the hull's
geometry at the attitude instead.

In waves:

```json
{ "type": "result", "id": "<study id>", "study": "<study id>", "kind": "waves",
  "froude": …, "speed": …, "attitude": { "sinkage": …, "trim_deg": … },
  "calm": "<the calm-water study it is held about>",
  "peaks": { "heave_peak": …, "heave_peak_lambda": …, "sea_accel_bow": …, … },
  "seakeeping": { "gm_t": …, "headings": [ { "heading": 180, "points": [ … ], "sea": … } ], … },
  "sections": "<its calm-water study's>", "solver_version": "…" }
```

A calm-water result's `sections` is the id of its hulls cut at the attitude
it solved (or held). A result in waves has its calm-water study's.

Each of the `seakeeping.headings.0.points` holds `lambda`, `omega`,
`omega_e`, the complex RAOs `heave`, `pitch`, `sway`, `roll` and `yaw` as
`[re, im]`, and the added resistance `raw_gb`.

### sections

A case's hulls as the solver cut them at one attitude: written for every
equilibrium, by `case` (at rest) and by `run` (each calm-water study at its
solved or held attitude).

```json
{ "type": "sections", "id": "…", "case": "<case id>",
  "attitude": { "sinkage": 0.0085, "trim_rad": 0.0015, "trim_deg": 0.087 },
  "hulls": [ { "placement": { "x": 0.0, "y": 0.0 }, "transom": true,
               "stations": [ { "x": 0.18, "z0": 0.0, "beam": 0.186, "depth": 0.0087,
                               "radii": [1.0, 0.9998, …] }, … ] } ] }
```

Each station holds its section exactly as the solver integrates it: rays from
the section's top on the centreplane, at the `n = radii.len()`
Chebyshev–Lobatto angles below the horizontal

```text
θ_k = ¼π (1 − cos(πk / (n − 1))),    k = 0 … n−1,
```

reaching the shell at the scaled distances `radii[k]`. The section is the curve

```text
y(θ) = beam · R(θ) · cos θ,    z(θ) = z0 + depth · R(θ) · sin θ,
```

where `y` is the half-breadth from the hull's centreplane, `z` the depth below
the water, and `R` the polynomial of degree `n − 1` through the radii. It's
best evaluated by barycentric interpolation, as
`michell_geometry::iges::PolarSection::point` does. A station with no radii
lies past the hull's tip. With `transom`, the aft station is a transom. A
catamaran has both demihulls, each with its own placement.

The id is that of the case and the attitude, so a case at rest and every study
held there share one record. Sections only go one way: they can't be re-posed
(that needs the surface above the water), so the hull's `geometry` stays the
source. For e12 at 121 stations of 33 rays, a record is 83 KB.

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
| `wake [-o FILE] [--range M] [--title]` | result* → SVG | the free surface from above |
| `pressure [-o FILE] [--range CP] [--title]` | result* → SVG | the pressure on the hulls from below |
| `profile [-o FILE] [--wave-scale K] [--title]` | result* → SVG | the hull at its attitude, with the wave along its side |
| `prop --d-max L [--shafts --wake --thrust-deduction --blades LIST --depth --keller-k --top-froude F …]` | result* → prop* | the best B-series propeller for each result's speed and thrust |
| `prop --thrust T --speed V --d-max L …` | → prop | the same, for a thrust and speed given outright |
| `match [--rank-by power\|mass\|price] [--direct-only] [--max-mass --max-od --max-price --vendor --mapped-only …]` | prop* → drive* | the motors that can drive each prop, ranked |
| `motors [filters]` | → motor* | the motor database |
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
  followed are `hull`, `case`, `study`, `parent`, `calm` and `sections`, so
  `--explode sections.hulls.0.stations` on a result gives a row per station.
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

### Propellers and motors

`prop` and `match` are a port of propopt's web app (`crates/propeller`, from
its `web/propcore.js` and `web/motorcore.js`; see that crate's docs for the
models). Golden tests hold the port to the original's answers.

```sh
# the hull pipeline's resistance, through to motors
boatmath study --froude 0.4 < cat.case.json | boatmath run \
  | boatmath prop --d-max 12in --top-froude 0.5 \
  | boatmath match --rank-by mass --max-mass 25 \
  | boatmath table rank motor motor.vendor P_elec ratio motor.mass_kg

# or a thrust and speed outright, as on the web page
boatmath prop --thrust 1kN --speed 8kn --d-max 16in | boatmath match | head -5
```

**`prop`** takes each calm-water result's speed and the thrust its resistance
asks for, T = R_t / (1 − t). It splits that across `--shafts`, by default one
per hull, so a catamaran has two. For each rpm it finds the best diameter,
blade-area ratio and blade count, with pitch solved to hold the thrust, under
Keller's and Burrill's cavitation limits. The cheapest point on that curve is
the answer.

- `--wake` and `--thrust-deduction` default to 0.
- A second operating point the same propeller must reach comes from
  `--top-froude` (the same study's result at that Froude number, already run)
  or from `--top-speed` and `--top-thrust`. It's a constraint, not an
  objective.
- Quantities take units: `16in`, `400mm`, `8kn`, `1kN`, `50kgf`.

```json
{ "type": "prop", "id": "…", "result": "<result id> | null",
  "inputs": { "speed": 4.1, "thrust": 630, "shafts": 2, "wake": 0, "thrust_deduction": 0,
              "d_min": 0.04, "d_max": 0.3048, "blades": [2,3,4,5,6,7], "depth": 0.3,
              "keller_k": 0.2, "cavitation": true, "re_correct": true, "strict_ear": false,
              "top": null },
  "shaft": { "V_A": …, "T": 315, "top_T": null },
  "boat": { "R_T": …, "P_E": …, "eta_H": 1, "QPC": … },
  "feasible": true,
  "best": { "rpm": 1089, "D": 0.305, "PD": 0.87, "EAR": 0.30, "Z": 2, "eta0": 0.778,
            "P_shaft": 1490, "Q": …, "J": …, "K_T": …, "K_Q": …, "Cth": …, "etaIdeal": …,
            "sigma": …, "kellerOk": true, "burrillOk": true, "top": …, "perZ": [ … ], … },
  "best_unconstrained": null,
  "window": { "rpm_lo": …, "rpm_hi": … }, "feasible_rpm": { "lo": …, "hi": … },
  "curve": [ { "rpm": …, "ok": true, "P_shaft": …, "Z": …, "D": …, … }, { "rpm": …, "ok": false } ] }
```

The design fields keep propcore.js's names. `curve` is the best propeller at
each of 180 shaft speeds, from the cheapest out to 1.6× its power (or the
whole feasible range with `--full-range`), so `plot --explode curve -x rpm
-y P_shaft --by Z` draws it.

**`match`** walks each motor along a prop's curve. At each point it picks the
reduction that costs least at the battery, and keeps the motor's best point.
Continuous ratings apply at the design point (`--peak` for peak) and peak
ratings at the second point, through the same gearbox. It writes one
`drive` per motor that can do the job, best first. By default that's the best
winding of each family (`--all-windings` for every one), ranked by electrical
power, or by `--rank-by mass` or `price`, where a motor that doesn't publish
the figure goes last. The filters (`--vendor`, `--mapped-only`,
`--single-only`, `--max-mass`, `--max-od`, `--max-price`) are hard limits,
which a motor not publishing the figure fails.

```json
{ "type": "drive", "id": "…", "prop": "<prop id>", "motor": "zapi_gsm309_4a_48v", "rank": 1,
  "rpm": 1055, "P_shaft": 5645,
  "propeller": { "D": 0.406, "PD": 0.81, "EAR": 0.31, "Z": 2, "eta0": 0.729 },
  "P_elec": 6125, "eta": …, "etaWithController": …, "etaDrive": 0.922,
  "ratio": 1, "gearEta": 1, "motorRpm": …, "motorQ": …,
  "V_bus": 20.4, "I_dc": …, "I_arms": …, "load": …, "extrapolated": false,
  "top": null, "tier": "B", "settings": { … } }
```

**Motors** come from propopt's database, vendored unchanged in
`crates/propeller/data/motors.js`: 159 motors from nine makers, each with a
`tier` saying where its efficiency comes from (see that directory's README).
`boatmath motors` lists them as records, and a field path follows `motor`
into the database, so `motor.mass_kg` works on a drive. `match --motors FILE`
reads another database, or motor records (`boatmath motors | jq …`).

### Pictures

`wake`, `pressure` and `profile` draw calm-water results as SVG, in metres
with x forward, to scale. A picture too thin to read at true scale has its
short axis stretched, and the axis label says by how much. With several
results on stdin, `-o` is a pattern naming each picture by `{id8}` or
`{froude}` (`wake-{froude}.svg`).

- `wake`: ζ(x, y) from above, as a heatmap from trough (blue) to crest (red),
  neutral at still water, with the hulls' waterlines drawn over it. The
  colour scale is ±`--range` metres, by default the 99th percentile of |ζ|,
  so a spike at a bow doesn't wash it out.
- `pressure`: C_p on each hull's surface from below (port at the bottom),
  on the pressure solver's grid of stations × depths. The scale is likewise
  ±`--range`.
- `profile`: the first hull side-on at its attitude, its geometry re-posed
  there (topsides and all), with the still waterline and the wave along its
  side: ζ at the waterline half-breadth, and on the centreplane beyond the
  ends. `--wave-scale` exaggerates the wave only, and the picture says so.

## Code

- **`crates/boatmath`** holds what the CLI and the web app share: `params.rs`
  (`CaseParams`, `StudyParams`), `LoftRequest`, `platform.rs` (statics, wave
  GZ, waves), the calm-water flow and `SOLVER_VERSION`. `boatmath-web`
  re-exports it and keeps the SQLite store, the worker and the pages, so both
  front ends compute the same way.
- `crates/boatmath/src/sections.rs`: a case cut at an attitude, as JSON.
- **`crates/propeller`**: the B-series propeller search and the motor
  models, ported from propopt's web cores, with its motor database vendored.
- **`crates/boatmath-cli`** (binary `boatmath`): `store.rs`, `records.rs`
  (hull, case, study, sections), `run.rs`, `path.rs`, `plot.rs`, `views.rs`
  (wake, pressure, profile), `props.rs` (prop, drive, motor), `units.rs`,
  `list.rs`.
- `crates/boatmath/src/native.rs`: geometry as JSON, read from a file or
  re-posed. `boatmath`'s computations take its bytes wherever they take a
  hull file's, so the web app could store it too.
- Platforms are single hulls and catamarans (`span`), as in the web app.

## Not yet

- `gz-waves`: quasi-static GZ in a regular wave (`platform::wave_gz`).
- `show`: open a hull, case or result in the 3-D viewer.
- Pictures of results in waves (RAOs are a `plot --explode` away already).
- propopt's Python-only work: the BEM correction for an extended or scanned
  geometry (`--geometry`), and its own hull models (the hull pipeline
  replaces those).
- Estimating the wake fraction and thrust deduction from the hull.
- Warm starts: the web worker starts each equilibrium from the nearest
  speed already solved. `run` starts every one from scratch.
- A study in waves is held about the calm-water study at the default grid.
  A calm study at another grid shares its attitude but is not looked for.

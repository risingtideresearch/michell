# Seakeeping in `michell`: what was built, how it validates, and what it is good for

A record of the work of 26–27 September 2026: the split of the workspace into
geometry, thin-ship and seakeeping crates, the strip-theory seakeeping built on
top, and every comparison made against published data. Numbers here are the
ones the code produced at the time; the regression tests named below hold the
important ones.

## 1. Summary

- **Workspace.** `michell-geometry` (hulls cut from IGES/STL into sections,
  hydrostatics, equilibrium), `michell` (thin-ship theory: Michell resistance,
  squat, near field, spectrum, friction) and `michell-seakeeping` (strip
  theory). Seakeeping depends on geometry only; it shares with the thin-ship
  code the hull geometry and, through `--dynamic`, the speed-dependent attitude.
- **Seakeeping.** Linear strip theory in the frequency domain, all five
  unrestrained modes — sway, heave, roll, pitch, yaw — at forward speed, for a
  single hull or a rigid platform of several (catamaran, trimaran, proa);
  mean added resistance by two methods (Gerritsma–Beukelman radiated energy,
  Maruo far-field momentum); irregular-sea statistics. Command line:
  `michell seakeeping`.
- **Validation, in one line each.**
  - Heave: good against Journée's Wigley experiments (coefficients 10–15%,
    wave loads 5–15%, head-sea RAOs close including resonance peaks).
  - Pitch at speed: poor — the same discrepancies Journée reported for his own
    strip codes on these hulls.
  - Added resistance: GB 2–6× high at Fn 0.3–0.4; Maruo within 2× there but ~½
    at Fn 0.2 and on beamy hulls; neither ranks hulls reliably.
  - Sway/heave/roll section hydrodynamics: 1–6% against SEAWAY and 1–12%
    against Vugts' cylinder experiments.
  - Three-dimensional lateral (sway/roll/yaw) response: checked only by limits,
    reciprocity and symmetry — no experimental comparison yet.
- **For comparing hull forms.** Trust heave and vertical motion near the
  centre of gravity; treat pitch and bow motions as a rough guide for large
  differences; do not decide on added resistance above Fn ≈ 0.3. Multihull roll
  (driven by demihull heave) is expected to be sound; monohull roll resonance
  needs an empirical damping model (only a fraction-of-critical stand-in
  exists).

## 2. Workspace split

| crate | contents | depends on |
|---|---|---|
| `michell-geometry` | IGES/STL import, `SectionalHull`, hydrostatics (now including transverse waterplane inertia and centre-of-buoyancy depth), hydrostatic/dynamic equilibrium (`float`), B-splines, closed-form moments, `Conditions`, `Placement` | — |
| `michell` | thin-ship theory; the sectional amplitude and near-field transforms as extension traits `SectionalWave` / `NearFieldTransforms` | geometry |
| `michell-seakeeping` | strip-theory seakeeping | geometry |
| `michell-cli`, `michell-web` | front ends | all |

The thin-ship equilibrium solver reaches geometry through
`float::DynamicModel`; seakeeping uses the sectional hull's own transforms
(`SectionalHull::x_transform`) and section curves.

## 3. What the seakeeping crate computes

### 3.1 Section hydrodynamics (`green`, `section2d`)

- **2-D pulsating source**, deep water, frequency domain, in closed form
  through the complex exponential integral (principal branch, verified against
  direct quadrature of the principal-value integral). E₁ uses the power series
  (including far out, near the negative real axis, where its error is scaled
  away by `e^{Re w}`), the continued fraction elsewhere, and the asymptotic
  series beyond `|w| = 40`.
- **Frank close-fit panel method** on each station's section curve, symmetric
  (heave) and antisymmetric (sway, roll about the waterline centre) problems:
  added mass, damping, far-field wave amplitudes, radiation potentials, and the
  diffraction problem's sources.
- **Wave forces per section:** Froude–Krylov by panel integration, diffraction
  by the Haskind relation (checked against the solved diffraction problem).
- **Irregular frequencies removed** by an interior rigid lid of sources
  (Ohmatsu's method, interior vertical velocity zero). A lid holding the
  potential at zero was tried first: it removes them too but clashes with the
  hull's potential at the waterline corner and converges far too slowly. The
  damping-versus-radiated-energy check still runs on every solve and bridges a
  failure by interpolation.

### 3.2 Strip theory (`strip`)

The platform's motions about G — sway, heave, roll, pitch, yaw — are assembled
generically. Each station's local sway, heave and roll follow from the
platform's through a kinematic map `T(x)` (hull offset `y` and G's height
included), and its hydrodynamic force is the forward-speed operator
`D = iω − U∂ₓ` wrapped round its complex added mass `c = a − ib/ω`:

```text
K = ∫ (iωT + UT′)ᵀ c (iωT − UT′) dx + U T_Aᵀ c_A (iωT_A − UT′_A),   F_rad = −K η
```

with the aft station's section as the transom end term. This reproduces the
Salvesen–Tuck–Faltinsen (1970) coefficients; for transom hulls it also carries
`(U²/ω²)a_A` in A₃₅ and `(U²/ω²)b_A` in B₃₅, which the earlier hand-written
heave–pitch table (written from memory) lacked. On e12 at Fn 0.4 that moves
heave by about 2% at λ/L 1.2; transom-free hulls are unchanged. Neither form
was checked against a published derivation — see §7.

Excitation: closed-form Froude–Krylov for the vertical modes (from the
sectional hull's depth integrals, with the oblique-sea transverse variation as
a station-by-station correction), panel Froude–Krylov for the lateral ones,
Haskind diffraction at the encounter frequency with STF's speed terms.
Restoring from the waterplane integrals, GM_T for roll. Multihulls act as one
rigid platform without hull-to-hull wave interaction; asymmetric platforms
couple roll to heave and pitch automatically.

### 3.3 Added resistance and irregular seas

- **Gerritsma–Beukelman**: radiated energy of each strip's motion relative to
  the Smith-reduced wave.
- **Maruo far field**: momentum of the whole wave pattern, from a Kochin
  function built from the stations' sources (heave sources at each section's
  relative vertical velocity plus its diffraction sources), with the
  forward-speed `k₁`/`k₂` wave systems. The Kochin normalisation was checked by
  the radiation energy on a slender hull; the integrand stops at the stations'
  Nyquist wavenumber (beyond it the assembly only aliases). Head and following
  seas only.
- **Irregular seas**: Bretschneider and JONSWAP spectra; significant heave,
  pitch and point accelerations; mean added resistance by both methods.

### 3.4 Command line

`michell seakeeping <hull>[@x=,y=]... (--speed U | --froude F) [--heading DEG]
[--lambda A:B:STEP] [--vcg Z] [--kxx K] [--kzz K] [--roll-damping ZETA]
[--sea hs=H,tp=T[,gamma=G]] [--dynamic] [--csv]` prints heave, pitch, sway,
roll and yaw RAOs, both added-resistance estimates, GM_T and the natural roll
period (with added inertia), and sea-state statistics. `--dynamic` first floats
the platform at its thin-ship dynamic sinkage and trim.

Example results for e12 (single demihull unless noted):

| case | result |
|---|---|
| Fn 0.4 head seas, Bretschneider Hs 0.5 m Tp 4 s | significant heave 0.21 m, pitch 4.2°, bow acceleration 6.9 m/s²; mean added resistance 109 N (GB), 40 N (far field) |
| catamaran, 2.8 m span, VCG 0.8 m | GM_T 14.9 m (almost all from spacing) |
| demihull, VCG 0.2 m | GM_T 0.07 m, natural roll period 2.3 s with added inertia |
| demihull, beam seas at rest, at roll resonance | roll RAO 19 with potential damping only, 6.8 with 5% of critical added |

## 4. Validation

### 4.1 Internal checks (all in the test suite)

- Source: principal value against quadrature; gradient against finite
  differences; outgoing far field; E₁ series/fraction/asymptotic agreement.
- Sections: damping equals the energy the far field carries (heave, sway, roll,
  and the sway–roll cross damping); Haskind equals the solved diffraction force
  (heave; sway and roll in beam and oblique seas); a semicircle's exact
  infinite-frequency heave added mass `ρπR²/2` (first-order convergence);
  a semicircle rolling about its centre moves no water; reciprocity of the
  sway–roll coefficients; convergence with panel count; the lid converges to
  the plain method where that one is sound, and removes the irregular
  frequency where it is not.
- Hull: Froude–Krylov against brute-force quadrature of the Wigley hull (head,
  oblique, beam seas); long-wave limits (heave → 1, pitch → wave slope; in beam
  seas sway → 1, roll → slope, no yaw); zero-speed reciprocity of the five-mode
  coefficients; far-apart catamaran twins move exactly like one hull; mirrored
  proas and mirrored oblique seas give mirrored motions; Kochin damping equals
  ∫b dx on a slender hull; roll hydrostatics exact on the Wigley hull.

### 4.2 Journée's four Wigley hulls in head waves (report 0909, 1992)

Hulls `η = (1 − ζ²)(1 − ξ²)(1 + 0.2ξ²) + α ζ²(1 − ζ⁸)(1 − ξ²)⁴`, L = 3 m,
d = 0.1875 m, B = 0.3 m (I, III) or 0.6 m (II, IV), α = 1 (I, II) or 0 (III,
IV), k_yy = 0.75 m, Fn 0.2–0.4; built exactly from this formula (displacements
match the report to 0.5%). The report's pitch normalisation is
`θ_a/(2πζ_a/L)`, not `θ_a/(kζ_a)`.

| quantity | agreement |
|---|---|
| heave added mass A₃₃ | ~10–15% over mid frequencies |
| heave damping B₃₃ | good at low–mid frequency; about half the measurement at ω√(L/g) > 5 and speed-independent where the tank's rises with speed |
| heave force / pitch moment on the restrained hull | ~5–15% |
| heave RAO in head waves | close, resonance included (Wigley I, Fn 0.3, λ/L 1.25: 2.52 measured, 2.51 computed) |
| zero-speed heave / pitch | ~10% / 10–20% low |
| pitch added inertia A₅₅ at speed | ~30% low |
| pitch damping B₅₅ at speed | grows with U² where the tank shows none (Wigley III, ω′ 2.2: computed 0.108 / 0.130 / 0.162 at Fn 0.2 / 0.3 / 0.4; measured 0.074 / 0.073 / 0.066) |
| sway–pitch coupling B₃₅ | follows strip theory's `U·A₃₃`, which the measurements do not |

Journée found the same pitch and coupling discrepancies with his own Frank-
and Ursell-based strip codes (report 1275, 2001), and showed the measurements
self-consistent (motions computed from the measured coefficients and wave loads
reproduce the measured motions), so these are strip theory's limits on these
hulls rather than faults of this code.

**Added resistance** — peak `R_aw/(ρgζ²B²/L)`, measured / GB / Maruo:

| hull, Fn | measured | GB | Maruo |
|---|---|---|---|
| I, 0.2 | 31.8 | 35.0 | 16.4 |
| II, 0.2 | 11.1 | 11.1 | 4.8 |
| III, 0.2 | 22.2 | 24.0 | 12.0 |
| IV, 0.2 | 13.5 | 10.0 | 4.2 |
| I, 0.3 | 27.6 | 51.4 | 18.0 |
| III, 0.3 | 20.2 | 49.2 | 19.4 |
| IV, 0.3 | 19.3 | 15.8 | 4.3 |
| I, 0.4 | 14.9 (noisy) | 91.8 | 26.2 |
| III, 0.4 | 27.7 | 71.5 | 22.7 |

Speed trend on Wigley I: measured peaks fall (32 → 28 → 15 over Fn 0.2–0.4); GB
rises steeply (35 → 51 → 92); Maruo rises gently (16 → 18 → 26). The far field
counts no momentum in the forward-scattered wave and lets the sections' waves
interfere, which GB (a sum of strips' sideways-radiated energy) cannot.

**Ranking the four hulls** (does the model order them as the tank does?):

| Fn | measure | measured | GB | Maruo |
|---|---|---|---|---|
| 0.2 | heave peak | I > II > IV > III | same | same |
| 0.2 | pitch peak | I > III > IV > II | I > III > II > IV | — |
| 0.2 | added-resistance area | I > III > IV > II | same | I > III > II > IV |
| 0.3 | heave peak | I > III > IV | same | same |
| 0.3 | added-resistance peak | I > III > IV | same | III > I > IV |
| 0.3 | added-resistance area | IV > I > III | I > III > IV | I > III > IV |
| 0.4 | pitch / added resistance | III > I | I > III | I > III |

### 4.3 Vugts' cylinders in beam waves (via SEAWAY report 1213, §2.1)

Circle and rectangles of B/d 2, 4, 8, as tested by Vugts (1970) and plotted in
Journée's SEAWAY validation report; points and SEAWAY's curves were read off
the report's vector figures (`python/tools/vugts/`). Median differences
relative to each quantity's scale:

| section | a₂₂ | b₂₂ | a₃₃ | b₃₃ | F₂ | F₃ | F₄ | a₄₄ | b₄₄ |
|---|---|---|---|---|---|---|---|---|---|
| circle, vs tank / SEAWAY | 4% / 0% | 3 / 0 | 1 / 0 | 3 / 1 | 3 / 1 | 5 / 0 | 12 / 1 | (zero in theory) | (zero in theory) |
| B/d 2 | 2 / 1 | 7 / 1 | 7 / 1 | 6 / 0 | 3 / 1 | 2 / 1 | 1 / 1 | 30 / 1 | 16 / 4 |
| B/d 4 | 5 / 2 | 1 / 2 | 3 / 1 | 7 / 1 | 5 / 1 | 2 / 1 | 12 / 6 | 28 / 6 | 15 / 7 |
| B/d 8 | 8 / 1 | 4 / 2 | 2 / 0 | 9 / 0 | 7 / 1 | 4 / 0 | 10 / 1 | — | 13 / 1 |

Couplings a₂₄, b₂₄, a₄₂, b₄₂: B/d 2 (G on the waterline) 2–20% against the tank
and 1% against SEAWAY. Phases match SEAWAY to 2–9° (median) once its
conventions are applied. Deviations from the tank are mostly shared with
SEAWAY: low-frequency heave added mass (the tank was 2 m deep; this solver is
deep water) and the circle's high-frequency sway damping.

**Conventions found in report 1213** (not stated there; inferred by fitting):

- heave-force phase: as here; sway-force phase: this code's + 180° (opposite
  sway sign); roll-moment phase: 90° − this code's (against the wave slope).
- for B/d 4 and 8 (G above the waterline, OG/d = 1 and 3) the plotted sway–roll
  couplings match this code's taken **about O with reversed sign**, while the
  roll moment and roll added mass match this code's **about G**. Unresolved.

A regression test (`validation::sections_match_vugts_cylinders_in_beam_waves`)
holds the circle and the B/d 2 rectangle to Vugts' points within 15% of each
quantity's scale.

## 5. Problems found and fixed along the way

- **Self-influence of tiny panels.** A panel's own-midpoint normal offset came
  from roundoff, whose sign flipped the principal value by ±π: graded sections
  (millimetre panels at waterline and keel) got a diagonal of 0 or 1 instead of
  ½, and a curve-resampled semicircle's damping was five times too large. Now
  zero by construction; a test guards curve-resampled sections.
- **Irregular frequencies.** First bridged by the energy check; later found to
  leak spurious interior modes into the Kochin function at other wavenumbers
  (spikes to −210 in added resistance), which led to the rigid lid.
- **E₁ near the negative real axis.** The continued fraction crawled (550 ms per
  section at 24 rad/s); series and asymptotic branches brought it to 7 ms.
- **Kochin aliasing.** Maruo's k₁ branch reaches wavenumbers far past the
  stations' Nyquist limit; the integrand is now cut there.
- **Wrong expectations corrected.** Measured Wigley added-resistance peaks are
  20–30, not the 5–10 first recalled; the report's pitch normalisation is
  `2π/L`, not `k`; a first Kochin energy test failed because a wide hull's 2-D
  sources radiate fore and aft (the check needs `kL ≫ 1 ≫ kB`).

## 6. Using it to compare hull forms

- **Reliable:** heave and vertical motion near G (ordering and magnitude),
  multihull roll (demihull heave), the wave loads on the hull.
- **Rough guide:** pitch and bow motions and accelerations — look for large
  differences only; strip theory's pitch inertia and damping at speed are off.
- **Not for decisions above Fn ≈ 0.3:** added resistance (both methods);
  report both as a bracket.
- **Monohull roll resonance:** needs `--roll-damping` (a few percent of
  critical) until an empirical model exists; off resonance roll is sound.
- **Out of the validated regime:** e12-type speeds (Fn ≥ 0.4), transom sterns,
  catamaran gap interaction. Small variations of one hull are more trustworthy
  than comparisons across hull families, since systematic errors largely
  cancel.

## 7. Open questions and next steps

- Empirical roll damping (Ikeda-type) for monohulls.
- The transom end terms: this code's operator form and the earlier table
  differ by `(U²/ω²)a_A`, `(U²/ω²)b_A` in A₃₅, B₃₅; neither was checked against
  a published derivation (the SEAWAY theoretical manual, report 1370, was
  downloaded but its equations did not survive text extraction).
- The report-1213 coupling convention for sections with G above the waterline.
- Three-dimensional lateral response against experiment (report 1213's
  containership "Nedlloyd Dejima" roll case needs hull lines not available).
- A short-wave added-resistance correction; better pitch coefficients at speed
  (the root of the pitch and added-resistance errors); finite water depth;
  hull-to-hull interaction for catamarans; data for the high-Fn, transom,
  multihull regime of the boats this is meant for.

## 8. Reproducing

```sh
# Wigley comparison (data: 0909-DUT-92.zip, unzipped)
cargo run --release -p michell-seakeeping --example journee_wigley -- /path/to/wigley
# Vugts cylinders (report 1213 PDF; pip install pymupdf)
cd python/tools/vugts && python digitize.py 9 10 11 12 13 14 15 16
cargo run --release -p michell-seakeeping --example vugts > ours.tsv && python compare_vugts.py
# tests
cargo test --release -p michell-seakeeping
```

## 9. Sources

Read in full or in the parts used:

- J.M.J. Journée, *Experiments and Calculations on 4 Wigley Hull Forms in Head
  Waves*, Delft University of Technology, Ship Hydromechanics Laboratory,
  Report 0909, May 1992 (reprinted 2003), and its ASCII data file
  `0909-DUT-92.zip`. Archived:
  <https://web.archive.org/web/2007id_/http://www.ocp.tudelft.nl/mt/journee/Files/PapersReports/0909-DUT-92.pdf>,
  <https://web.archive.org/web/2007id_/http://www.ocp.tudelft.nl/mt/journee/Files/PapersReports/0909-DUT-92.zip>.
- J.M.J. Journée, *Discrepancies in Hydrodynamic Coefficients of Wigley Hull
  Forms*, MARIND 2001, Varna; Report 1275-P.
  <https://web.archive.org/web/2007id_/http://www.ocp.tudelft.nl/mt/journee/Files/PapersReports/1275-MARIND-01.pdf>.
- J.M.J. Journée, *Verification and Validation of Ship Motions Program SEAWAY*,
  Report 1213a, February 2001 (§2.1: Vugts' cylinders).
  <https://web.archive.org/web/2007id_/http://www.ocp.tudelft.nl/mt/journee/Files/PapersReports/1213-ValidationSEAWAY.pdf>.
- J.M.J. Journée and A.P. van 't Veer, *First Order Wave Loads in Beam Waves*,
  ISOPE 1995; Report 1027-P (method-versus-method only; not used as data).
  <https://web.archive.org/web/2007id_/http://www.ocp.tudelft.nl/mt/journee/Files/PapersReports/1027-ISOPE-95.pdf>.
- S. Liu, H. Liang, X. Chen, *Advancing the Understanding of Added Resistance in
  Waves Through Fourier-Kochin Theory*, 40th IWWWFB, 2025 — the statement of
  Maruo's forward-speed formula used here.
  <http://www.iwwwfb.org/Abstracts/iwwwfb40/IWWWFB40_30.pdf>.
- Journée's archived homepage, the index of all the above:
  <https://web.archive.org/web/20070101024210/http://www.ocp.tudelft.nl/mt/journee/>.

Downloaded but not used in substance: J.M.J. Journée and L.J.M. Adegeest,
*Theoretical Manual of Strip Theory Program "SEAWAY for Windows"*, Report 1370,
2003 (equations garbled in text extraction); J.M.J. Journée, *Quick Strip
Theory Calculations in Ship Design*, PRADS 1992, Report 0902-P. An attempt at
M. Kashiwagi, *Study on the Wave-Induced Steady Force and Moment* (J. Soc. Nav.
Arch. Japan 173, 1993, J-STAGE) failed (downloads truncated).

Cited from memory, not retrieved — the methods as implemented:

- N. Salvesen, E.O. Tuck, O. Faltinsen, *Ship Motions and Sea Loads*, Trans.
  SNAME 78, 1970 (strip theory, forward-speed and transom terms).
- W. Frank, *Oscillation of Cylinders in or below the Free Surface of Deep
  Fluids*, NSRDC Report 2375, 1967 (close-fit source method).
- S. Ohmatsu, *On the Irregular Frequencies in the Theory of Oscillating
  Bodies in a Free Surface*, Papers of the Ship Research Institute 48, 1975
  (the lid).
- J.H. Vugts, *The Hydrodynamic Forces and Ship Motions in Waves*, PhD thesis,
  Delft, 1970 (the cylinder experiments, via report 1213).
- J. Gerritsma, W. Beukelman, *Analysis of the Resistance Increase in Waves of
  a Fast Cargo Ship*, International Shipbuilding Progress 19, 1972.
- H. Maruo, *The Drift of a Body Floating on Waves*, J. Ship Research 4, 1960,
  and later forward-speed work; J.N. Newman, *The Exciting Forces on Fixed
  Bodies in Waves*, J. Ship Research 6, 1962 (Haskind–Newman relation).
- Wehausen & Laitone, *Surface Waves*, 1960 (the 2-D and 3-D pulsating
  sources); O.M. Faltinsen, *Sea Loads on Ships and Offshore Structures*, 1990.

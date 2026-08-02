# Insel-Wigley validation method record

This document fixes the experimental interpretation before any resistance values
are transcribed and before any `michell` prediction is computed. Page references
give the printed page followed by the PDF page in parentheses. The offset in the
Insel thesis is +10 PDF pages; the offset in Ship Science Report 72 is +2.

## Source correction and study boundary

The supplied study brief describes resistance tables beginning near thesis page
173 and describes Ship Science Report 72 as a thin-ship comparison with the same
Wigley experiments. Direct inspection of rendered pages shows that neither
description is correct:

- Insel's Tables 1--5, printed pages 173--175 (PDF 183--185), contain model
  particulars, the test matrix, and wave-analysis settings. They do not contain
  resistance ordinates. The C2 resistance results exist as plotted marker series
  in the figure appendix.
- Report 72 uses the round-bilge, transom-stern NPL models reported in Ship
  Science Report 71, not Insel's C2 Wigley hull. It remains useful precedent for
  the theory and for judging discrepancies, but it supplies no directly
  comparable C2 prediction series.

Consequently, the external-data study will digitize visibly plotted experimental
markers, report calibration and reading uncertainty, and retain figure-level
provenance. No curve will be represented as a tabulated source value. Report 72
will be discussed as methodological precedent only. This is an explicit source
limitation, not a reason to infer missing values.

## C2 hull identity and geometry

Insel designates the model **C2** when free to trim and sink and **C2-FX** when
fixed. He describes it as a glass-reinforced-plastic mathematical form with
parabolic waterlines and parabolic cross-sections, "known as a Wigley Model"
(printed 66--67, PDF 76--77). Table 1 gives the following particulars (printed
173, PDF 183):

| quantity | C2 / C2-FX |
|---|---:|
| waterline length, `L` | 1.800 m |
| `L/B` | 10.000 |
| `B/T` | 1.600 |
| `L / volume^(1/3)` | 7.116 |
| block coefficient, `C_B` | 0.444 |
| prismatic coefficient, `C_P` | 0.667 |
| midship coefficient, `C_M` | 0.667 |
| listed wetted surface, `WS` | 0.482 m^2 |

The arithmetic consequences are `B = 0.180 m` and `T = 0.1125 m`. These are
derived from the tabulated ratios, not separately printed values. The C2 body
plan is Figure 132 (printed 242, PDF 252), followed by model photographs in
Figures 133--134 (printed 243, PDF 253).

The thesis does not print an algebraic hull equation. Its description and
coefficients identify the ordinary product-parabola Wigley form

```text
y(x,z) = +/- (B/2) [1 - (2x/L)^2] [1 - (z/T)^2],
-L/2 <= x <= L/2,  -T <= z <= 0.
```

This equation yields `C_B = 4/9` and `C_P = C_M = 2/3`, matching Table 1 to
the shown precision. It is therefore a documented reconstruction, not a
verbatim equation from Insel. Xu et al. (2025) print the same equation for the
Insel Wigley geometry; Srinakaew (2017) later calls this family "Wigley III."
The original source's C2/C2-FX designation is used here to avoid importing a
later naming convention.

No prototype or model-scale ratio is given for C2. All calculations will use
the physical 1.800 m tank model. No full-scale ship or scale factor will be
invented.

## Separation and placements

`S` is the transverse distance between **demihull centrelines**, not the clear
gap between inner surfaces (thesis nomenclature, printed vii / PDF 9; test
description, printed 52 / PDF 62). Table 2 lists `S/L = 0.2, 0.3, 0.4, 0.5`
for both C2 and C2-FX (printed 173, PDF 183), equivalent to `S/B = 2, 3, 4, 5`
for this hull. The computational placements will therefore be `y = -S/2` and
`y = +S/2`.

Table 2 gives a documented C2 range `Fn = 0.2--0.95` and C2-FX range
`Fn = 0.2--0.8`. Narrative passages instead say `Fn = 0.1--0.95` (printed 67,
PDF 77) or approximately `0.1--0.9`, with some high-speed series curtailed
(printed 110, PDF 120). Digitization will follow the markers actually present,
while preserving these differing source statements.

## Fixed and free attitude

All hull families were tested free to trim and sink; C2 was additionally tested
fixed in trim and sink (printed 64, PDF 74). Every test restrained surge, sway,
yaw, and roll. C2-FX used vertical restraint at the tow post and dedicated fixed
trim fittings. Fixed-attitude C2-FX is the primary validation target because
the current linear Michell calculation does not solve the running attitude.
Free C2 data will be retained as a documented out-of-model comparison rather
than mixed into the fixed-attitude score.

Insel reports trim to about +/-0.05 degrees and sinkage to about +/-0.1 mm,
with positive trim bow-up and positive sinkage downward (printed 65, PDF 75).
The thesis says fixed/free C2 comparisons were intended to expose the effects
of running attitude and changing mean wetted surface (printed 110, PDF 120).

## Coefficient normalization and wetted surface

Insel defines the total-resistance coefficient with denominator
`rho * S_W * U^2 / 2` (printed 64, PDF 74). Table 1's `WS = 0.482 m^2` is the
isolated C2 hull value. The thesis does not contain an explicit sentence saying
that a catamaran uses exactly `2 * WS`, nor does the nomenclature label `WS`
as static or running area. The component/interference equations require a
common coefficient basis, and the conventional reading is total static wetted
area: `0.482 m^2` for the monohull and `0.964 m^2` for two identical demihulls.
That reading is recorded as an inference and will be sensitivity-checked rather
than presented as a direct quotation.

Report 72 is explicit for its separate NPL series: `A` is static wetted surface
area and resistance coefficients use `rho * A * v^2 / 2` (printed 2, PDF 4).
This supports, but does not independently prove, the C2 interpretation.

The prediction harness will report dimensional `R_w`, the area used, and the
resulting `C_W`, so the normalization remains auditable.

## Insel--Molland component decomposition

For an isolated demihull, Insel writes (printed 38, PDF 48)

```text
C_T = (1 + k) C_F + C_W.
```

For a catamaran, Equation 3.12 is (printed 37, PDF 47)

```text
(C_T)_cat = (1 + phi k) sigma C_F + tau C_W
          = (1 + beta k) C_F + tau C_W.
```

Here `sigma` is the frictional-interference factor, `phi` is the form-resistance
interference factor, and `tau` is the wave-resistance interference factor. The
combined viscous multiplier is defined by

```text
(1 + phi k) sigma = 1 + beta k.
```

For the experimental separation, Insel gives the operative ratio (printed 124,
PDF 134)

```text
tau = [C_T - (1 + beta k) C_F]_cat
      / [C_T - (1 + k) C_F]_mono.
```

Thus `tau` compares catamaran and isolated-demihull wave-resistance components
on the common coefficient convention. It is not `C_T(cat) / C_T(mono)`, and it
must not be confused with the separate wave-pattern factor denoted `mu` in the
direct wave-pattern analysis. Insel assumes `beta` constant with speed for a
given hull and separation and fits the smallest compatible viscous multiplier
using the wave-pattern separation (printed 107, PDF 117).

The `michell` library's reported multihull interference is a wave-resistance
ratio. It is therefore comparable to `tau` only after the experimental
normalization and numerator/denominator definitions have been made identical.
Michell contains no model for `beta`.

## Tank, turbulence stimulation, and corrections

The University of Southampton Ship Science towing tank was 60.0 m long, 3.7 m
wide, and 1.85 m deep. The carriage maximum was 4.6 m/s; tests reached about
4.2 m/s, with roughly 20 m for acceleration, 20 m at steady speed, and a
15.24 m measurement section (printed 62, PDF 72).

The largest model cross-section was below 0.5% of the tank section, so viscous
blockage was neglected. The estimated wall-interference effect for C2 was under
1%, so no wall correction was applied. Shallow-water corrections were also
omitted: theory suggested less than 2% for the monohull and 4% for the
catamaran, except near `0.95 < Fn_h < 1.02` (printed 63, PDF 73). This makes the
upper end of the Froude range a documented tank-effect caution.

Water temperature ranged from 13 to 18.5 degrees C. Results were standardized
to freshwater at 15 degrees C by

```text
C_T = C_T(TC) - C_F(TC) + C_F(15 C),
```

and the stated resistance accuracy was +/-0.02 N, considered satisfactory for
speeds above 1.0 m/s (printed 64, PDF 74). Report 72 supplies `rho = 1000 kg/m^3`,
`nu = 1.141e-6 m^2/s`, and `g = 9.80665 m/s^2` for 15-degree freshwater
(printed 2, PDF 4); the harness will use these values.

C2 turbulence studs were 3.2 mm in diameter, 2.5 mm high, 90 mm aft of the
stem, and spaced 25 mm apart. Insel made no explicit correction for stud drag
or the laminar length ahead of the studs, assuming the two effects cancelled
(printed 67, PDF 77).

Insel judges the useful experimental range mainly above `Fn = 0.25` and below
about `Fn = 0.9`; small forces degrade the low end, while acceleration,
shallow-water effects, wave breaking, and wave-pattern scatter affect the high
end (printed 120--121, PDF 130--131).

## Wave-pattern resistance measurement

`C_WP` was obtained with a multiple-longitudinal-cut method using four wave
probes and a simultaneous matrix solution, not by a transverse-cut integral
(printed 69--84, PDF 79--94). The method was chosen specifically to avoid the
infinite-trace assumption and explicit finite-trace truncation correction of a
direct Fourier transform. Four longitudinal traces were taken at transverse
fractions `y/W = 3/10, 1/3, 2/5, 4/9`; a 25 m record yielded about 15 m of
effective trace. A 40-harmonic representation was generally adequate above
1 m/s, but the useful harmonic count fell at higher speeds.

The source's defensible measurement limitations are finite recorded wave
length and harmonic/matrix resolution, viscous attenuation between hull and
probe, nonlinear or breaking waves, and poor resolution when C2 wave heights
fell below about 3 mm. The probe accuracy was approximately +/-0.05 mm. Insel
specifically associates low-speed C2 underprediction with the very small wave
height (printed 82, PDF 92).

The pre-registration must not describe this experiment as suffering
"transverse-cut truncation": that mechanism does not match the apparatus.
Xu et al. (2025) independently identify inadequate experimental longitudinal
wave record as a high-speed limitation when revisiting the C2 `S/L = 0.3`
case with CFD-generated wave cuts.

## Report 72 theoretical precedent

Report 72 uses linearized far-field Kelvin-source theory in a finite-width
channel. Thin-ship assumptions derive source strengths from geometry; hulls are
represented by parametric cubic splines and discretized into 20 depthwise by
50 lengthwise source panels. Point sources were retained after comparison with
constant-strength panels because a fine point-source mesh was faster with
similar results (printed 4, PDF 6).

Measured running trim and sinkage were inserted by translating and rotating the
theoretical hull and regenerating the panels. The waterline remained planar;
the local wave profile along the hull was not modelled. Several transom closure
treatments were tested, with a transverse sink line retained and only modest
benefit from alternatives (printed 5, PDF 7).

The experimental basis is Report 71's NPL round-bilge transom series, including
Model 5b catamarans at `S/L = 0.2, 0.3, 0.4, 0.5`, not the C2 Wigley series.
Agreement was reasonable and trends useful above approximately `Fn = 0.4`,
especially for the more slender models; discrepancies increased below that
speed before the transoms ran clear and for the fuller model. The authors treat
the method as qualitatively useful for parametric work while attributing the
remaining ceiling to linearized slender-body assumptions (printed 7--8,
PDF 9--10).

The exact Report 72 comparison plots are Figures 3--8 (printed 17--20,
PDF 19--22) and the virtual-appendage `C_WP` comparisons in Figures 12--14
(printed 22--24, PDF 24--26). They will not be overlaid on C2 results.

## Experimental series to digitize

Because the source has no resistance tables, each dataset below is a plotted
marker series. The digitization log will identify the axis calibration, blind
passes, mismatch resolution, and visual curve check for every figure.

Primary fixed-attitude C2-FX series:

| configuration | total-resistance figure | wave-pattern figure | printed page | PDF page |
|---|---:|---:|---:|---:|
| monohull | 135 | 136 | 244 | 254 |
| `S/L = 0.2` | 137 | 138 | 245 | 255 |
| `S/L = 0.3` | 139 | 140 | 246 | 256 |
| `S/L = 0.4` | 141 | 142 | 247 | 257 |
| `S/L = 0.5` | 143 | 144 | 248 | 258 |

The fixed-attitude separation summaries are Figures 146--147 (printed
249--250, PDF 259--260), and detailed fixed interference plots are Figures
155--158 (printed 254--255, PDF 264--265). Experimental `tau` will preferably
be reconstructed from the independently digitized component series at common
Froude numbers, rather than sampled from a smoothed summary curve.

Secondary free-attitude C2 series occupy Figures 161--192 (printed 258--273,
PDF 268--283). The separation summaries for `C_W` and `C_WP` are Figures
182--183 (printed 268--269, PDF 278--279), the detailed `tau` plots are Figures
189--192 (printed 272--273, PDF 282--283), and the running trim/sinkage summary
is Figures 193--194 (printed 274, PDF 284). These data will be archived, but
they are outside the primary fixed-attitude score unless the plotted markers
can be read without guessing.

Cross-model `tau` comparisons are Figures 347--350 (printed 352--353,
PDF 362--363). They are context, not an additional C2 dataset.

## Independent source audit

A second Codex task independently inspected rendered pages without OCR or text
extraction and without seeing any digitized CSV. It confirmed the geometry,
test matrix, coefficient definitions, fixed/free policy, tank treatment,
figure locations, lack of resistance tables, and Report 72 source mismatch.
The audit made no repository edits and computed no predictions.

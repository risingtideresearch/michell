# Insel C2 free-attitude digitization pass B

Blind independent pass B. Source: rendered thesis JPEGs only; no OCR, text extraction, repository digitized CSV, or solver output was used.

## Counts

- S/L=0.2, CT: 44 visible markers accepted
- S/L=0.2, CWP: 41 visible markers accepted
- S/L=0.2, sinkage_over_draught: 19 visible markers accepted
- S/L=0.2, trim_deg: 43 visible markers accepted
- S/L=0.3, CT: 48 visible markers accepted
- S/L=0.3, CWP: 37 visible markers accepted
- S/L=0.3, sinkage_over_draught: 19 visible markers accepted
- S/L=0.3, trim_deg: 34 visible markers accepted
- S/L=0.4, CT: 64 visible markers accepted
- S/L=0.4, CWP: 25 visible markers accepted
- S/L=0.4, sinkage_over_draught: 42 visible markers accepted
- S/L=0.4, trim_deg: 49 visible markers accepted
- S/L=0.5, CT: 32 visible markers accepted
- S/L=0.5, CWP: 40 visible markers accepted
- S/L=0.5, sinkage_over_draught: 30 visible markers accepted
- S/L=0.5, trim_deg: 51 visible markers accepted
- monohull, CT: 47 visible markers accepted
- monohull, CWP: 47 visible markers accepted
- monohull, sinkage_over_draught: 25 visible markers accepted
- monohull, trim_deg: 42 visible markers accepted

## Calibration and uncertainty

Axes were calibrated from the rendered plot-border intersections listed in `calibrations.csv`. Pixel-to-value conversion is affine on each printed axis. The stated uncertainty is digitization uncertainty only, not the experimental uncertainty reported by Insel. Marker centers use the centroid of the open-square white center (CT and trim) or the dense printed core (CWP and sinkage). Replicates are preserved; only white islands less than 7 px apart inside a single open marker were collapsed.

## Ambiguous/unreadable material

- Figures 161-179 contain no numeric tables; these are plot-marker readings.
- Several low-Froude CWP and sinkage markers merge with the fitted curve or neighboring replicates in the scan. A dense-core component larger than the 500-pixel inspection ceiling could not be separated without subjective tracing and is omitted rather than guessed. The raw overlays and candidate CSVs identify these merged regions.
- Figure 175 has multiple fitted curves passing through the measured series near the peak. Only filled-marker dense cores meeting the stated shape/area test are accepted; borderline thin curve fragments are omitted.
- Open-square centers split by a fitted line were collapsed only within 7 px. Distinct marker centers outside that radius remain separate, including same-speed vertical replicates.
- Figure 164/168/172/176/180 theory-comparison curves were not transcribed. Figure 180 is outside the requested Figure 161-179 range.

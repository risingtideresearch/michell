# Vugts' cylinders, from SEAWAY validation report 1213

Digitises the vector figures of J.M.J. Journée, *Verification and
Validation of Ship Motions Program SEAWAY* (Delft report 1213a, 2001, §2.1:
Vugts' 1970 experiments on horizontal cylinders in beam waves) and compares
them with this repository's 2-D section solver.

The report is freely distributed by its author; an archived copy is at
`https://web.archive.org/web/2007id_/http://www.ocp.tudelft.nl/mt/journee/Files/PapersReports/1213-ValidationSEAWAY.pdf`.

```sh
pip install pymupdf
python digitize.py 9 10 11 12 13 14 15 16          # -> vugts_digitised.json
cargo run --release -p michell-seakeeping --example vugts > ours.tsv
python compare_vugts.py
```

`digitize.py` finds each subplot from its tick labels, calibrates the axes,
and takes the experiment markers (small filled shapes) and SEAWAY's curve
(long polylines) in data coordinates. See the `vugts` example for the
report's phase and coupling conventions.

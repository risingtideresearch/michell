# Insel-Wigley validation sources

Primary documents for the experimental validation study. Source PDFs are
downloaded into `data/raw/`, which is ignored by git. Only factual table data
transcribed with the protocol in `DIGITIZATION.md` may be committed.

## Primary sources

1. Mustafa Insel (1990), *An Investigation into the Resistance Components of
   High Speed Displacement Catamarans*, PhD thesis, Department of Ship Science,
   University of Southampton.
   - Repository record: https://eprints.soton.ac.uk/462776/
   - Version of record: https://eprints.soton.ac.uk/462776/1/457354.pdf
   - Local raw filename: `data/raw/insel-1990-thesis.pdf`

2. A. F. Molland, J. F. Wellicome, and P. R. Couser (1994), *Resistance
   Experiments on a Systematic Series of High Speed Displacement Catamaran
   Forms: Variation of Length-Displacement Ratio and Breadth-Draught Ratio*,
   Ship Science Report 71, University of Southampton, Southampton, UK, 84 pp.,
   ISSN 0140-3818.
   - Repository record: https://eprints.soton.ac.uk/46442/
   - Version of record: https://eprints.soton.ac.uk/46442/1/071.pdf
   - Local raw filename: `data/raw/ship-science-71.pdf`

3. A. F. Molland, J. F. Wellicome, and P. R. Couser (1994), *Theoretical
   Prediction of the Wave Resistance of Slender Hull Forms in Catamaran
   Configurations*, Ship Science Report 72, University of Southampton,
   Southampton, UK, 36 pp., ISSN 0140-3818.
   - Repository record: https://eprints.soton.ac.uk/46441/
   - Version of record: https://eprints.soton.ac.uk/46441/1/072.pdf
   - Local raw filename: `data/raw/ship-science-72.pdf`

## Reproduction

Run `./fetch.sh` from this directory. The script downloads only missing files
and leaves partial downloads under temporary names until each transfer
completes successfully.

## Supplementary sources checked during method extraction

These sources are not substitutes for the primary experimental record. They
were fetched from the linked publisher or institutional repository and used to
adjudicate geometry naming, later reuse of the data, and limitations of the
wave-cut method.

1. S. Srinakaew (2017), *A Numerical Study of Resistance Components of High-
   Speed Catamarans and the Scale Effects on Form Factor*, PhD thesis,
   University of Southampton.
   - Repository record: https://eprints.soton.ac.uk/420755/
   - Version of record:
     https://eprints.soton.ac.uk/420755/1/Final_Thesis_Sarawuth_Srinakaew.pdf
   - Note: later uses the name "Wigley III" and reproduces some Insel total and
     residual-resistance data, but its `C_W` tables are friction-line residuals,
     not measured wave-pattern `C_WP`.

2. X. Xu, Z. Zou, and X. Chen (2025), "Mechanism analysis and prediction of
   longitudinal cut wave pattern resistance based on CFD simulation," *Journal
   of Ocean Engineering and Science* 10(2), 271--288.
   - Article: https://doi.org/10.1016/j.joes.2023.07.001
   - Note: prints the Wigley equation and revisits Insel's fixed-attitude
     `S/L = 0.3` case, including experimental wave-record limitations.

3. P. R. Couser, J. F. Wellicome, and A. F. Molland (2000), *An Improved Method
   for the Theoretical Prediction of the Wave Resistance of Transom-Stern Hulls
   Using a Slender Body Approach*, Ship Science Report 125, University of
   Southampton.
   - Version of record: https://eprints.soton.ac.uk/46408/1/125ShipScience_Report.pdf
   - Note: includes a later Wigley monohull thin-ship comparison; it does not
     provide the missing C2 catamaran tables.

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

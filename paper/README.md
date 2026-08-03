# Paper A — Journal of Ship Research manuscript

Compile the manuscript with the repository's chosen LaTeX engine:

```sh
mkdir -p ../output/pdf
tectonic main.tex --outdir ../output/pdf --keep-logs --keep-intermediates
```

Copy the generated `main.pdf` to
`endpoint-bickley-michell-jsr.pdf` for review. The manuscript follows the
public SNAME journal template's letter-paper, Times, 10-point, two-column,
unnumbered-section, title, abstract, and keyword conventions. The final PDF is
kept under `output/pdf/`; auxiliary files remain ignored. The Zenodo DOI and
corresponding-author submission details remain manual pre-submission items.

# Paper draft

Compile the manuscript with the repository's chosen LaTeX engine:

```sh
mkdir -p ../output/pdf
tectonic main.tex --outdir ../output/pdf --keep-logs --keep-intermediates
```

Rename the generated `main.pdf` to
`endpoint-bickley-michell-draft.pdf` for distribution.  The final PDF is kept
under `output/pdf/`; auxiliary files remain ignored. The author list, target
journal, and repository archive URL are placeholders that must be resolved
before submission.

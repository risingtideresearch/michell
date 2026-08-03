# Archive and DOI checklist

The repository is prepared for a Zenodo-backed release, but no remote state or
DOI has been created by this work.

1. Confirm that `https://github.com/risingtideresearch/michell` is public and
   that the intended repository-wide MIT licensing is documented to the
   foundation's satisfaction. The Rust package manifests already declare MIT;
   the repository currently has no root `LICENSE` file.
2. Confirm the final branch is clean and rerun the exact commands in
   `paper/COMPLETION_REPORT.md`.
3. Inspect the local annotated tag `paper-a-jsr-v1`; move or recreate it only
   if the manuscript changes after review.
4. Push the final branch and tag manually. No push is authorized or performed
   by the submission-preparation task.
5. Enable the GitHub repository in Zenodo, create a release from the tag, and
   verify that Zenodo imports `.zenodo.json` with Rob Story and Avi Bryant as
   creators and Rising Tide Research Foundation as their affiliation.
6. Reserve or mint the DOI. Replace `10.5281/zenodo.REPLACE-ME` in
   `paper/main.tex`, rebuild with Tectonic, inspect the final PDF, and commit
   that single substitution.
7. Verify the DOI resolves to an archive containing the tagged source,
   `paper/results/`, the measurement script, and the manuscript PDF before
   submitting to JSR.

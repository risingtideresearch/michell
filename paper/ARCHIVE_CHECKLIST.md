# Archive and DOI checklist

The repository is prepared for a Zenodo-backed release, but no remote state or
DOI has been created by this work.

1. Confirm that `https://github.com/risingtideresearch/michell` is public and
   that the intended repository-wide MIT licensing is documented to the
   foundation's satisfaction. The Rust package manifests already declare MIT;
   the repository currently has no root `LICENSE` file.
2. Confirm the final branch is clean and rerun the exact commands in
   `paper/COMPLETION_REPORT.md`.
3. Create the local annotated tag `paper-a-jsr-v2` only after the final source,
   tests, numerical evidence, completion report, and reviewed PDF are committed.
   Treat that tag as immutable: never move, force-update, or reuse it. If the
   paper changes, create a new monotonically numbered tag.
4. Export and checksum exactly that tag:

   ```sh
   mkdir -p output/archive
   git archive --format=tar.gz --prefix=michell-paper-a-jsr-v2/ \
     -o output/archive/michell-paper-a-jsr-v2.tar.gz paper-a-jsr-v2
   shasum -a 256 output/archive/michell-paper-a-jsr-v2.tar.gz \
     > output/archive/michell-paper-a-jsr-v2.tar.gz.sha256
   ```

   Confirm that the tarball contains no `.git` directory. Extract it into a
   fresh temporary directory and run `paper/reproduce_measurements.sh` there;
   `ARCHIVE_REVISION` records the exported commit when Git metadata is absent.
5. Push the final branch and tag manually. No push is authorized or performed
   by the submission-preparation task.
6. Enable the GitHub repository in Zenodo, create a release from the tag, and
   verify that Zenodo imports `.zenodo.json` with Rob Story and Avi Bryant as
   creators and Rising Tide Research Foundation as their affiliation.
7. Reserve or mint the DOI. Until that external action occurs, the manuscript
   must say that the archive *will* be deposited and must display
   `10.5281/zenodo.REPLACE-ME` only as a placeholder. After minting, replace the
   placeholder, rebuild with Tectonic, inspect the PDF, and create a new tag;
   do not move `paper-a-jsr-v2`.
8. Verify the DOI resolves to an archive containing the tagged source,
   `paper/results/`, the measurement script, and the manuscript PDF before
   submitting to JSR.

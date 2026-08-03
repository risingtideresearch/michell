#!/usr/bin/env python3
"""Attach audited identities to the existing blind C2 theory traces."""

from __future__ import annotations

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent
PASSES = HERE / "passes"
OUTPUT_A = PASSES / "final_c2_360_362_pass_a_attributed.tsv"
OUTPUT_B = PASSES / "final_c2_360_362_pass_b_attributed.tsv"
COMBINED = HERE / "theory_c2_attributed_359_362.csv"

POSITIONS = {
    359: "solid C2: oscillatory left branch; isolated upper broad crest near Fn 0.45",
    360: "solid C2: oscillatory left branch; leftmost principal crest near Fn 0.41",
    361: "solid C2: oscillatory left branch; principal crest near Fn 0.40",
    362: "solid C2: small left oscillations; sharp principal crest near Fn 0.40",
}


def read(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter=delimiter))


def write(path: Path, rows: list[dict[str, str]], delimiter: str = ",", comments: tuple[str, ...] = ()) -> None:
    with path.open("w", newline="") as handle:
        for comment in comments:
            handle.write(f"# {comment}\n")
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def attributed_pass(path: Path, pass_name: str) -> list[dict[str, str]]:
    rows = []
    for source in read(path, delimiter="\t"):
        figure = int(source["figure"])
        if figure not in (360, 361, 362):
            continue
        rows.append({
            "pass": pass_name,
            "figure": source["figure"],
            "pdf_page": source["pdf_page"],
            "printed_page": source["printed_page"],
            "separation_over_length": source["separation_over_length"],
            "model": "C2",
            "hull": "Wigley hull",
            "line_style": "solid",
            "relative_position_description": POSITIONS[figure],
            "fn_anchor": source["fn_anchor"],
            "pixel_x": source["pixel_x"],
            "pixel_y": source["pixel_y"],
            "pixel_x_left": source["pixel_x_left"],
            "pixel_x_right": source["pixel_x_right"],
            "pixel_y_top": source["pixel_y_top"],
            "pixel_y_bottom": source["pixel_y_bottom"],
            "identity_basis": "LEGEND-AUDIT-H1.md; identity added without changing archived pixel reading",
        })
    return rows


def main() -> None:
    pass_a = attributed_pass(PASSES / "pass_a_theory.tsv", "A")
    pass_b = attributed_pass(PASSES / "pass_b_theory.tsv", "B")
    write(OUTPUT_A, pass_a, delimiter="\t")
    write(OUTPUT_B, pass_b, delimiter="\t")

    combined: list[dict[str, str]] = []
    for source in read(HERE / "theory_family_359.csv"):
        if source["model"] != "C2":
            continue
        combined.append({
            "figure": source["figure"],
            "source_printed_page": source["source_printed_page"],
            "source_pdf_page": source["source_pdf_page"],
            "separation_over_length": source["separation_over_length"],
            "model": "C2", "hull": "Wigley hull", "line_style": "solid",
            "relative_position_description": POSITIONS[359],
            "fn_anchor": source["fn_anchor"], "fn": source["fn"], "tau": source["tau"],
            "fn_digitization_uncertainty": source["fn_digitization_uncertainty"],
            "tau_digitization_uncertainty": source["tau_digitization_uncertainty"],
            "pass_count": source["pass_count"],
            "source_archive": "theory_family_359.csv",
        })
    for source in read(HERE / "theory_interference.csv"):
        figure = int(source["figure"])
        if figure not in (360, 361, 362):
            continue
        combined.append({
            "figure": source["figure"],
            "source_printed_page": source["source_printed_page"],
            "source_pdf_page": source["source_pdf_page"],
            "separation_over_length": source["separation_over_length"],
            "model": "C2", "hull": "Wigley hull", "line_style": "solid",
            "relative_position_description": POSITIONS[figure],
            "fn_anchor": source["fn_anchor"], "fn": source["fn"], "tau": source["tau"],
            "fn_digitization_uncertainty": source["fn_digitization_uncertainty"],
            "tau_digitization_uncertainty": source["tau_digitization_uncertainty"],
            "pass_count": source["pass_count"],
            "source_archive": "theory_interference.csv",
        })
    combined.sort(key=lambda row: (int(row["figure"]), float(row["fn_anchor"])))
    write(
        COMBINED,
        combined,
        comments=(
            "source: Mustafa Insel, 1990 PhD thesis, Figures 359--362",
            "source_location: printed pages 358--359; PDF pages 368--369",
            "identity: solid C2 (Wigley hull), independently audited before this attribution",
        ),
    )
    print(f"wrote {len(pass_a)} pass-A and {len(pass_b)} pass-B attributed rows")
    print(f"wrote {len(combined)} reconciled C2 rows to {COMBINED}")


if __name__ == "__main__":
    main()

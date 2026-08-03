#!/usr/bin/env python3
"""Reconcile the two source-only Insel C2 theory-curve traces."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
PASSES = HERE / "passes"
OUTPUT = HERE / "theory_interference.csv"
MISMATCHES = PASSES / "theory_mismatches.csv"

MAX_FN_DIFFERENCE = 0.015
MAX_TAU_DIFFERENCE = 0.08
BASE_FN_UNCERTAINTY = 0.003
BASE_TAU_UNCERTAINTY = 0.02


def read_pass(path: Path) -> dict[tuple[int, str], dict[str, str]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {(int(row["figure"]), row["fn_anchor"]): row for row in rows}
    assert len(result) == len(rows), f"duplicate anchors in {path}"
    return result


def calibrate(row: dict[str, str]) -> tuple[float, float]:
    x = float(row["pixel_x"])
    y = float(row["pixel_y"])
    x_left = float(row["pixel_x_left"])
    x_right = float(row["pixel_x_right"])
    y_top = float(row["pixel_y_top"])
    y_bottom = float(row["pixel_y_bottom"])
    fn = 0.1 + (x - x_left) * 0.9 / (x_right - x_left)
    tau = (y_bottom - y) * 2.5 / (y_bottom - y_top)
    return fn, tau


def main() -> None:
    pass_a = read_pass(PASSES / "pass_a_theory.tsv")
    pass_b = read_pass(PASSES / "pass_b_theory.tsv")
    assert pass_a.keys() == pass_b.keys(), "the passes must use common anchors"

    admitted = []
    omitted = []
    for key in sorted(pass_a):
        a = pass_a[key]
        b = pass_b[key]
        for field in ("figure", "pdf_page", "printed_page", "separation_over_length", "fn_anchor"):
            assert a[field] == b[field], (key, field, a[field], b[field])
        fn_a, tau_a = calibrate(a)
        fn_b, tau_b = calibrate(b)
        fn_difference = abs(fn_a - fn_b)
        tau_difference = abs(tau_a - tau_b)
        common = {
            "figure": a["figure"],
            "source_printed_page": a["printed_page"],
            "source_pdf_page": a["pdf_page"],
            "separation_over_length": a["separation_over_length"],
            "fn_anchor": a["fn_anchor"],
            "fn_a": f"{fn_a:.7f}",
            "tau_a": f"{tau_a:.7f}",
            "fn_b": f"{fn_b:.7f}",
            "tau_b": f"{tau_b:.7f}",
            "fn_pass_difference": f"{fn_difference:.7f}",
            "tau_pass_difference": f"{tau_difference:.7f}",
            "pass_a_x_px": a["pixel_x"],
            "pass_a_y_px": a["pixel_y"],
            "pass_b_x_px": b["pixel_x"],
            "pass_b_y_px": b["pixel_y"],
        }
        if fn_difference <= MAX_FN_DIFFERENCE and tau_difference <= MAX_TAU_DIFFERENCE:
            admitted.append({
                **common,
                "fn": f"{(fn_a + fn_b) / 2:.7f}",
                "tau": f"{(tau_a + tau_b) / 2:.7f}",
                "fn_digitization_uncertainty": f"{max(BASE_FN_UNCERTAINTY, fn_difference / 2):.7f}",
                "tau_digitization_uncertainty": f"{max(BASE_TAU_UNCERTAINTY, tau_difference / 2):.7f}",
                "pass_count": "2",
                "adjudication": "two-pass source-only match; ordinate averaged",
            })
        else:
            omitted.append({
                **common,
                "resolution": "omitted_after_source_only_reinspection",
                "reason": "pass_difference_exceeds_preregistered_admission_limit",
            })

    output_fields = [
        "figure", "source_printed_page", "source_pdf_page", "separation_over_length",
        "fn_anchor", "fn", "tau", "fn_digitization_uncertainty",
        "tau_digitization_uncertainty", "pass_count", "adjudication", "fn_a", "tau_a",
        "fn_b", "tau_b", "fn_pass_difference", "tau_pass_difference", "pass_a_x_px",
        "pass_a_y_px", "pass_b_x_px", "pass_b_y_px",
    ]
    with OUTPUT.open("w", newline="") as handle:
        handle.write("# source: Mustafa Insel, 1990 PhD thesis, Figures 359--362\n")
        handle.write("# source_location: printed pages 358--359; PDF pages 368--369\n")
        writer = csv.DictWriter(handle, fieldnames=output_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(admitted)

    mismatch_fields = [
        "figure", "source_printed_page", "source_pdf_page", "separation_over_length",
        "fn_anchor", "reason", "resolution", "fn_a", "tau_a", "fn_b", "tau_b",
        "fn_pass_difference", "tau_pass_difference", "pass_a_x_px", "pass_a_y_px",
        "pass_b_x_px", "pass_b_y_px",
    ]
    with MISMATCHES.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=mismatch_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(omitted)

    counts = Counter(row["figure"] for row in admitted)
    scoring_counts = Counter(
        row["figure"] for row in admitted if 0.20 <= float(row["fn"]) <= 0.80
    )
    print(f"admitted {len(admitted)} of {len(pass_a)} anchors; omitted {len(omitted)}")
    for figure in sorted(counts, key=int):
        print(f"Figure {figure}: {counts[figure]} admitted, {scoring_counts[figure]} in Fn 0.20--0.80")


if __name__ == "__main__":
    main()

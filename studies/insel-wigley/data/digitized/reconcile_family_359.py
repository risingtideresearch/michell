#!/usr/bin/env python3
"""Reconcile the two source-only full-family traces of Insel Figure 359."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
PASSES = HERE / "passes"
OUTPUT = HERE / "theory_family_359.csv"
MISMATCHES = PASSES / "h1_figure_359_mismatches.csv"
METRICS = HERE / "h1_figure_359_family_metrics.csv"
CHECKS = HERE / "h1_figure_359_family_checks.csv"

MAX_FN_DIFFERENCE = 0.015
MAX_TAU_DIFFERENCE = 0.08
BASE_FN_UNCERTAINTY = 0.003
BASE_TAU_UNCERTAINTY = 0.02


def read_pass(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {(row["model"], row["fn_anchor"]): row for row in rows}
    assert len(result) == len(rows), f"duplicate model/anchor in {path}"
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


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]], comments: list[str] = []) -> None:
    with path.open("w", newline="") as handle:
        for comment in comments:
            handle.write(f"# {comment}\n")
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    pass_a = read_pass(PASSES / "h1_figure_359_pass_a.tsv")
    pass_b = read_pass(PASSES / "h1_figure_359_pass_b.tsv")
    assert pass_a.keys() == pass_b.keys(), "the passes must use common model/anchor keys"

    admitted: list[dict[str, str]] = []
    omitted: list[dict[str, str]] = []
    for key in sorted(pass_a):
        a = pass_a[key]
        b = pass_b[key]
        for field in (
            "figure", "pdf_page", "printed_page", "separation_over_length",
            "model", "hull", "line_style", "fn_anchor",
        ):
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
            "model": a["model"],
            "hull": a["hull"],
            "line_style": a["line_style"],
            "relative_position_pass_a": a["relative_position_description"],
            "relative_position_pass_b": b["relative_position_description"],
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
                "reason": "pass_difference_exceeds_preregistered_admission_limit",
                "resolution": "omitted; no solver overlay consulted",
            })

    output_fields = [
        "figure", "source_printed_page", "source_pdf_page", "separation_over_length",
        "model", "hull", "line_style", "relative_position_pass_a",
        "relative_position_pass_b", "fn_anchor", "fn", "tau",
        "fn_digitization_uncertainty", "tau_digitization_uncertainty", "pass_count",
        "adjudication", "fn_a", "tau_a", "fn_b", "tau_b", "fn_pass_difference",
        "tau_pass_difference", "pass_a_x_px", "pass_a_y_px", "pass_b_x_px",
        "pass_b_y_px",
    ]
    write_rows(
        OUTPUT,
        output_fields,
        admitted,
        [
            "source: Mustafa Insel, 1990 PhD thesis, Figure 359",
            "source_location: printed page 358; PDF page 368",
            "source-only two-pass trace; no michell output used in extraction or reconciliation",
        ],
    )

    mismatch_fields = [
        "figure", "source_printed_page", "source_pdf_page", "separation_over_length",
        "model", "hull", "line_style", "relative_position_pass_a",
        "relative_position_pass_b", "fn_anchor", "reason", "resolution", "fn_a",
        "tau_a", "fn_b", "tau_b", "fn_pass_difference", "tau_pass_difference",
        "pass_a_x_px", "pass_a_y_px", "pass_b_x_px", "pass_b_y_px",
    ]
    write_rows(MISMATCHES, mismatch_fields, omitted)

    by_model: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in admitted:
        if 0.35 <= float(row["fn"]) <= 0.55:
            by_model[row["model"]].append(row)
    metric_rows: list[dict[str, str]] = []
    for model in ("C2", "C3", "C4", "C5"):
        peak = max(by_model[model], key=lambda row: float(row["tau"]))
        amplitude = float(peak["tau"]) - 1.0
        metric_rows.append({
            "figure": "359",
            "separation_over_length": "0.2",
            "model": model,
            "hull": peak["hull"],
            "line_style": peak["line_style"],
            "amplitude_window_fn_min": "0.35",
            "amplitude_window_fn_max": "0.55",
            "peak_fn": peak["fn"],
            "peak_tau": peak["tau"],
            "amplitude_max_tau_minus_one": f"{amplitude:.7f}",
            "fn_digitization_uncertainty": peak["fn_digitization_uncertainty"],
            "tau_digitization_uncertainty": peak["tau_digitization_uncertainty"],
        })
    write_rows(METRICS, list(metric_rows[0]), metric_rows)

    metrics = {row["model"]: row for row in metric_rows}
    check_rows: list[dict[str, str]] = []
    for left, right in (("C3", "C4"), ("C4", "C2"), ("C2", "C5")):
        left_value = float(metrics[left]["amplitude_max_tau_minus_one"])
        right_value = float(metrics[right]["amplitude_max_tau_minus_one"])
        left_uncertainty = float(metrics[left]["tau_digitization_uncertainty"])
        right_uncertainty = float(metrics[right]["tau_digitization_uncertainty"])
        if left_value + left_uncertainty < right_value - right_uncertainty:
            status = "VIOLATION"
        elif left_value - left_uncertainty >= right_value + right_uncertainty:
            status = "PASS"
        else:
            status = "UNCERTAINTY_OVERLAP_TIE"
        check_rows.append({
            "check": "broad_hump_amplitude_order",
            "registered_relation": f"{left} >= {right}",
            "left_model": left,
            "left_value": f"{left_value:.7f}",
            "left_uncertainty": f"{left_uncertainty:.7f}",
            "right_model": right,
            "right_value": f"{right_value:.7f}",
            "right_uncertainty": f"{right_uncertainty:.7f}",
            "central_margin": f"{left_value - right_value:.7f}",
            "status": status,
        })

    c2_fn = float(metrics["C2"]["peak_fn"])
    c2_uncertainty = float(metrics["C2"]["fn_digitization_uncertainty"])
    first_rbh = min(
        (metrics[model] for model in ("C3", "C4", "C5")),
        key=lambda row: float(row["peak_fn"]),
    )
    rbh_fn = float(first_rbh["peak_fn"])
    rbh_uncertainty = float(first_rbh["fn_digitization_uncertainty"])
    if c2_fn - c2_uncertainty > rbh_fn + rbh_uncertainty:
        hump_status = "VIOLATION"
    elif c2_fn + c2_uncertainty <= rbh_fn - rbh_uncertainty:
        hump_status = "PASS"
    else:
        hump_status = "UNCERTAINTY_OVERLAP_TIE"
    check_rows.append({
        "check": "broad_hump_fn_order",
        "registered_relation": "C2 <= earliest(C3,C4,C5)",
        "left_model": "C2",
        "left_value": f"{c2_fn:.7f}",
        "left_uncertainty": f"{c2_uncertainty:.7f}",
        "right_model": first_rbh["model"],
        "right_value": f"{rbh_fn:.7f}",
        "right_uncertainty": f"{rbh_uncertainty:.7f}",
        "central_margin": f"{rbh_fn - c2_fn:.7f}",
        "status": hump_status,
    })
    write_rows(CHECKS, list(check_rows[0]), check_rows)

    print(f"admitted {len(admitted)} of {len(pass_a)} anchors; omitted {len(omitted)}")
    for row in metric_rows:
        print(
            f"{row['model']}: A={row['amplitude_max_tau_minus_one']} "
            f"at Fn={row['peak_fn']} +/- {row['fn_digitization_uncertainty']}"
        )
    for row in check_rows:
        print(f"{row['registered_relation']}: {row['status']}")


if __name__ == "__main__":
    main()

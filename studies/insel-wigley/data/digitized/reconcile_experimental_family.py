#!/usr/bin/env python3
"""Reconcile the two source-only traces of experimental Figures 347--350."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
PASSES = HERE / "passes"
OUTPUT = HERE / "experimental_family_347_350.csv"
MISMATCHES = PASSES / "final_experimental_family_mismatches.csv"
METRICS = HERE / "final_experimental_family_metrics.csv"
CHECKS = HERE / "final_experimental_family_checks.csv"

MAX_FN_DIFFERENCE = 0.015
MAX_TAU_DIFFERENCE = 0.08
BASE_FN_UNCERTAINTY = 0.003
BASE_TAU_UNCERTAINTY = 0.02


def read_pass(path: Path) -> dict[tuple[int, str, str], dict[str, str]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {(int(row["figure"]), row["model"], row["fn_anchor"]): row for row in rows}
    assert len(result) == len(rows), f"duplicate figure/model/anchor in {path}"
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


def write(path: Path, rows: list[dict[str, str]], comments: tuple[str, ...] = ()) -> None:
    with path.open("w", newline="") as handle:
        for comment in comments:
            handle.write(f"# {comment}\n")
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    pass_a = read_pass(PASSES / "final_experimental_family_pass_a.tsv")
    pass_b = read_pass(PASSES / "final_experimental_family_pass_b.tsv")
    assert pass_a.keys() == pass_b.keys()
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
            "figure": a["figure"], "source_printed_page": a["printed_page"],
            "source_pdf_page": a["pdf_page"],
            "separation_over_length": a["separation_over_length"],
            "model": a["model"], "hull": a["hull"], "line_style": a["line_style"],
            "relative_position_pass_a": a["relative_position_description"],
            "relative_position_pass_b": b["relative_position_description"],
            "fn_anchor": a["fn_anchor"], "fn_a": f"{fn_a:.7f}",
            "tau_a": f"{tau_a:.7f}", "fn_b": f"{fn_b:.7f}",
            "tau_b": f"{tau_b:.7f}", "fn_pass_difference": f"{fn_difference:.7f}",
            "tau_pass_difference": f"{tau_difference:.7f}",
            "pass_a_x_px": a["pixel_x"], "pass_a_y_px": a["pixel_y"],
            "pass_b_x_px": b["pixel_x"], "pass_b_y_px": b["pixel_y"],
        }
        if fn_difference <= MAX_FN_DIFFERENCE and tau_difference <= MAX_TAU_DIFFERENCE:
            admitted.append({
                **common, "fn": f"{(fn_a + fn_b) / 2:.7f}",
                "tau": f"{(tau_a + tau_b) / 2:.7f}",
                "fn_digitization_uncertainty": f"{max(BASE_FN_UNCERTAINTY, fn_difference / 2):.7f}",
                "tau_digitization_uncertainty": f"{max(BASE_TAU_UNCERTAINTY, tau_difference / 2):.7f}",
                "pass_count": "2", "adjudication": "two-pass source-only match; ordinate averaged",
            })
        else:
            omitted.append({
                **common, "reason": "pass_difference_exceeds_preregistered_admission_limit",
                "resolution": "omitted; no theory or solver curve consulted",
            })
    write(
        OUTPUT,
        admitted,
        (
            "source: Mustafa Insel, 1990 PhD thesis, Figures 347--350",
            "source_location: printed pages 352--353; PDF pages 362--363",
            "principal-hump window only; source-only two-pass trace",
        ),
    )
    write(MISMATCHES, omitted)

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in admitted:
        grouped[(row["separation_over_length"], row["model"])].append(row)
    metric_rows: list[dict[str, str]] = []
    for separation in ("0.2", "0.3", "0.4", "0.5"):
        for model in ("C2", "C3", "C4", "C5"):
            peak = max(grouped[(separation, model)], key=lambda row: float(row["tau"]))
            metric_rows.append({
                "figure": peak["figure"], "separation_over_length": separation,
                "model": model, "hull": peak["hull"], "line_style": peak["line_style"],
                "peak_fn": peak["fn"], "peak_tau": peak["tau"],
                "amplitude_max_tau_minus_one": f"{float(peak['tau']) - 1.0:.7f}",
                "fn_digitization_uncertainty": peak["fn_digitization_uncertainty"],
                "tau_digitization_uncertainty": peak["tau_digitization_uncertainty"],
            })
    write(METRICS, metric_rows)

    metrics = {(row["separation_over_length"], row["model"]): row for row in metric_rows}
    check_rows: list[dict[str, str]] = []
    for separation in ("0.2", "0.3"):
        c2 = metrics[(separation, "C2")]
        npl = max(
            (metrics[(separation, model)] for model in ("C3", "C4", "C5")),
            key=lambda row: float(row["amplitude_max_tau_minus_one"]),
        )
        c2_value = float(c2["amplitude_max_tau_minus_one"])
        npl_value = float(npl["amplitude_max_tau_minus_one"])
        c2_uncertainty = float(c2["tau_digitization_uncertainty"])
        npl_uncertainty = float(npl["tau_digitization_uncertainty"])
        resolved_c2_above = c2_value - c2_uncertainty > npl_value + npl_uncertainty
        resolved_npl_above = npl_value - npl_uncertainty > c2_value + c2_uncertainty
        difference = c2_value - npl_value
        ratio = c2_value / npl_value if c2_value > 0.0 and npl_value > 0.0 else float("nan")
        material = resolved_c2_above and (difference >= 0.10 or ratio >= 1.15)
        comparable = resolved_c2_above and (difference >= 0.20 or ratio >= 1.34)
        check_rows.append({
            "separation_over_length": separation, "c2_amplitude": f"{c2_value:.7f}",
            "c2_uncertainty": f"{c2_uncertainty:.7f}", "npl_envelope_model": npl["model"],
            "npl_envelope_amplitude": f"{npl_value:.7f}",
            "npl_uncertainty": f"{npl_uncertainty:.7f}",
            "c2_minus_npl": f"{difference:.7f}", "c2_over_npl": f"{ratio:.7f}",
            "resolved_c2_above": str(resolved_c2_above).lower(),
            "resolved_npl_above": str(resolved_npl_above).lower(),
            "material_c2_excess": str(material).lower(),
            "comparable_c2_excess": str(comparable).lower(),
        })
    if any(row["resolved_npl_above"] == "true" for row in check_rows) or not any(
        row["material_c2_excess"] == "true" for row in check_rows
    ):
        outcome = "ANOMALOUS-VS-OWN-EXPERIMENT"
    elif sum(row["comparable_c2_excess"] == "true" for row in check_rows) >= 1 and all(
        row["material_c2_excess"] == "true" for row in check_rows
    ):
        outcome = "DISSOLVES"
    else:
        outcome = "MIXED"
    for row in check_rows:
        row["registered_outcome"] = outcome
    write(CHECKS, check_rows)

    print(f"admitted {len(admitted)} of {len(pass_a)} anchors; omitted {len(omitted)}")
    for row in check_rows:
        print(
            f"S/L={row['separation_over_length']}: C2={row['c2_amplitude']}, "
            f"NPL={row['npl_envelope_model']} {row['npl_envelope_amplitude']}, "
            f"D={row['c2_minus_npl']}"
        )
    print(f"registered outcome: {outcome}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Apply the unchanged theory gates to the bounded historical-variant probe."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np


STUDY = Path(__file__).resolve().parents[1]
SOURCE = STUDY / "data/digitized/theory_c2_attributed_359_362.csv"
PREDICTIONS = STUDY / "data/predictions/historical_variant_predictions.csv"
OUTPUT = STUDY / "data/analysis/historical_variant_scores.csv"
EPSILON = 1e-12

FIGURES = {0.2: 359, 0.3: 360, 0.4: 361, 0.5: 362}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def curve(rows: list[dict[str, str]], fn_field: str, tau_field: str) -> tuple[np.ndarray, np.ndarray]:
    ordered = sorted(rows, key=lambda row: float(row[fn_field]))
    return (
        np.array([float(row[fn_field]) for row in ordered]),
        np.array([float(row[tau_field]) for row in ordered]),
    )


def interpolate(data: tuple[np.ndarray, np.ndarray], query: np.ndarray) -> np.ndarray:
    values = np.interp(query, data[0], data[1])
    values[(query < data[0][0]) | (query > data[0][-1])] = np.nan
    return values


def grade(median: float, p90: float, hump_error: float, amplitude_ratio: float) -> str:
    components = [
        "agreement" if median <= 0.05 + EPSILON and p90 <= 0.10 + EPSILON
        else "partial" if median <= 0.10 + EPSILON and p90 <= 0.20 + EPSILON
        else "disagreement",
        "agreement" if hump_error <= 0.010 + EPSILON
        else "partial" if hump_error <= 0.020 + EPSILON
        else "disagreement",
        "agreement" if 0.90 - EPSILON <= amplitude_ratio <= 1.10 + EPSILON
        else "partial" if 0.80 - EPSILON <= amplitude_ratio <= 1.20 + EPSILON
        else "disagreement",
    ]
    if all(component == "agreement" for component in components):
        return "agreement"
    if any(component == "disagreement" for component in components):
        return "disagreement"
    return "partial"


def main() -> None:
    source_rows: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in read_csv(SOURCE):
        if row["model"] == "C2":
            source_rows[float(row["separation_over_length"])].append(row)

    prediction_rows: dict[tuple[str, float], list[dict[str, str]]] = defaultdict(list)
    for row in read_csv(PREDICTIONS):
        prediction_rows[(row["variant"], float(row["separation_over_length"]))].append(row)

    source = {
        separation: curve(rows, "fn_anchor", "tau")
        for separation, rows in source_rows.items()
    }
    predictions = {
        key: curve(rows, "fn", "interference")
        for key, rows in prediction_rows.items()
    }
    scoring_grid = np.round(np.arange(0.20, 0.8001, 0.01), 3)
    hump_grid = np.round(np.arange(0.35, 0.5501, 0.001), 3)
    scores: list[dict[str, object]] = []
    for (variant, separation), predicted_curve in sorted(predictions.items()):
        source_scoring = interpolate(source[separation], scoring_grid)
        predicted_scoring = interpolate(predicted_curve, scoring_grid)
        admitted = np.isfinite(source_scoring) & np.isfinite(predicted_scoring)
        points = int(np.count_nonzero(admitted))
        assert points >= 50
        signed_error = predicted_scoring[admitted] - source_scoring[admitted]
        absolute_error = np.abs(signed_error)

        source_hump = interpolate(source[separation], hump_grid)
        predicted_hump = interpolate(predicted_curve, hump_grid)
        source_index = int(np.nanargmax(source_hump))
        predicted_index = int(np.nanargmax(predicted_hump))
        source_hump_fn = float(hump_grid[source_index])
        predicted_hump_fn = float(hump_grid[predicted_index])
        hump_error = abs(predicted_hump_fn - source_hump_fn)
        source_hump_tau = float(source_hump[source_index])
        predicted_hump_tau = float(predicted_hump[predicted_index])
        amplitude_ratio = predicted_hump_tau / source_hump_tau
        median = float(np.median(absolute_error))
        p90 = float(np.quantile(absolute_error, 0.90))
        scores.append({
            "variant": variant,
            "source_figure": FIGURES[separation],
            "separation_over_length": separation,
            "points": points,
            "median_absolute_error": median,
            "p90_absolute_error": p90,
            "rms_error": float(np.sqrt(np.mean(signed_error**2))),
            "maximum_absolute_error": float(np.max(absolute_error)),
            "median_signed_error": float(np.median(signed_error)),
            "source_hump_fn": source_hump_fn,
            "predicted_hump_fn": predicted_hump_fn,
            "hump_fn_error": hump_error,
            "source_hump_tau": source_hump_tau,
            "predicted_hump_tau": predicted_hump_tau,
            "hump_amplitude_ratio": amplitude_ratio,
            "panel_verdict": grade(median, p90, hump_error, amplitude_ratio),
        })

    variants = sorted({str(row["variant"]) for row in scores})
    study_verdict = {
        variant: "hit"
        if all(
            row["panel_verdict"] == "agreement"
            for row in scores
            if row["variant"] == variant
        )
        else "miss"
        for variant in variants
    }
    for row in scores:
        row["four_panel_verdict"] = study_verdict[str(row["variant"])]

    OUTPUT.parent.mkdir(exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        handle.write("# gates: unchanged CRITERIA-THEORY.md; a hit requires agreement on every component of all four panels\n")
        writer = csv.DictWriter(handle, fieldnames=scores[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(scores)

    for variant in variants:
        panels = [row for row in scores if row["variant"] == variant]
        details = ", ".join(
            f"{row['separation_over_length']:.1f}:{row['panel_verdict']}"
            for row in panels
        )
        print(f"{variant}: {study_verdict[variant]} ({details})")


if __name__ == "__main__":
    main()

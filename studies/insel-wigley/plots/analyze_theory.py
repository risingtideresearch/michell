#!/usr/bin/env python3
"""Apply CRITERIA-THEORY.md and regenerate the Insel C2 theory comparison."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STUDY = Path(__file__).resolve().parents[1]
SOURCE = STUDY / "data/digitized/theory_interference.csv"
PREDICTIONS = STUDY / "data/predictions/theory_predictions.csv"
ANALYSIS = STUDY / "data/analysis"
PLOTS = Path(__file__).resolve().parent
CRITERIA_COMMIT = "7f588d8"
DIGITIZATION_COMMIT = "325aecf"
COMPARISON_EPSILON = 1e-12

CONFIGURATIONS = {
    "s_l_0_2": (0.2, 359),
    "s_l_0_3": (0.3, 360),
    "s_l_0_4": (0.4, 361),
    "s_l_0_5": (0.5, 362),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def source_curves() -> dict[float, dict[str, np.ndarray]]:
    grouped: dict[float, list[dict[str, str]]] = {
        separation: [] for separation, _ in CONFIGURATIONS.values()
    }
    for row in read_csv(SOURCE):
        grouped[float(row["separation_over_length"])].append(row)
    result = {}
    for separation, rows in grouped.items():
        rows.sort(key=lambda row: float(row["fn"]))
        result[separation] = {
            "fn": np.array([float(row["fn_anchor"]) for row in rows]),
            "tau": np.array([float(row["tau"]) for row in rows]),
            "uncertainty": np.array(
                [float(row["tau_digitization_uncertainty"]) for row in rows]
            ),
        }
    return result


def prediction_curves() -> dict[str, dict[str, np.ndarray]]:
    grouped: dict[str, list[dict[str, str]]] = {
        configuration: [] for configuration in CONFIGURATIONS
    }
    for row in read_csv(PREDICTIONS):
        grouped[row["configuration"]].append(row)
    result = {}
    for configuration, rows in grouped.items():
        rows.sort(key=lambda row: float(row["fn"]))
        assert len(rows) == 151
        assert all(
            row["outcome"] == "converged" and row["solo_outcome"] == "converged"
            for row in rows
        )
        result[configuration] = {
            "fn": np.array([float(row["fn"]) for row in rows]),
            "tau": np.array([float(row["interference"]) for row in rows]),
        }
    return result


def interpolate(curve: dict[str, np.ndarray], query: np.ndarray) -> np.ndarray:
    values = np.interp(query, curve["fn"], curve["tau"])
    outside = (query < curve["fn"][0]) | (query > curve["fn"][-1])
    values[outside] = np.nan
    return values


def pointwise_grade(median: float, p90: float) -> str:
    if median <= 0.05 + COMPARISON_EPSILON and p90 <= 0.10 + COMPARISON_EPSILON:
        return "agreement"
    if median <= 0.10 + COMPARISON_EPSILON and p90 <= 0.20 + COMPARISON_EPSILON:
        return "partial"
    return "disagreement"


def position_grade(error: float) -> str:
    if error <= 0.010 + COMPARISON_EPSILON:
        return "agreement"
    if error <= 0.020 + COMPARISON_EPSILON:
        return "partial"
    return "disagreement"


def amplitude_grade(ratio: float) -> str:
    if 0.90 - COMPARISON_EPSILON <= ratio <= 1.10 + COMPARISON_EPSILON:
        return "agreement"
    if 0.80 - COMPARISON_EPSILON <= ratio <= 1.20 + COMPARISON_EPSILON:
        return "partial"
    return "disagreement"


def verdict(grades: list[str]) -> str:
    if all(grade == "agreement" for grade in grades):
        return "agreement"
    if any(grade == "disagreement" for grade in grades):
        return "disagreement"
    return "partial"


def score(
    source: dict[float, dict[str, np.ndarray]],
    predictions: dict[str, dict[str, np.ndarray]],
) -> list[dict[str, object]]:
    scoring_grid = np.round(np.arange(0.20, 0.8001, 0.01), 3)
    hump_grid = np.round(np.arange(0.35, 0.5501, 0.001), 3)
    rows = []
    for configuration, (separation, figure) in CONFIGURATIONS.items():
        source_scoring = interpolate(source[separation], scoring_grid)
        predicted_scoring = interpolate(predictions[configuration], scoring_grid)
        valid = np.isfinite(source_scoring) & np.isfinite(predicted_scoring)
        assert np.count_nonzero(valid) >= 50
        signed_error = predicted_scoring[valid] - source_scoring[valid]
        absolute_error = np.abs(signed_error)

        source_hump = interpolate(source[separation], hump_grid)
        predicted_hump = interpolate(predictions[configuration], hump_grid)
        source_index = int(np.nanargmax(source_hump))
        predicted_index = int(np.nanargmax(predicted_hump))
        source_hump_fn = float(hump_grid[source_index])
        predicted_hump_fn = float(hump_grid[predicted_index])
        source_hump_tau = float(source_hump[source_index])
        predicted_hump_tau = float(predicted_hump[predicted_index])
        hump_fn_error = abs(predicted_hump_fn - source_hump_fn)
        hump_amplitude_ratio = predicted_hump_tau / source_hump_tau

        point_grade = pointwise_grade(
            float(np.median(absolute_error)), float(np.quantile(absolute_error, 0.90))
        )
        fn_grade = position_grade(hump_fn_error)
        amp_grade = amplitude_grade(hump_amplitude_ratio)
        rows.append({
            "configuration": configuration,
            "s_over_l": separation,
            "source_figure": figure,
            "points": int(np.count_nonzero(valid)),
            "median_absolute_error": float(np.median(absolute_error)),
            "p90_absolute_error": float(np.quantile(absolute_error, 0.90)),
            "rms_error": float(np.sqrt(np.mean(signed_error**2))),
            "maximum_absolute_error": float(np.max(absolute_error)),
            "median_signed_error": float(np.median(signed_error)),
            "pointwise_grade": point_grade,
            "source_hump_fn": source_hump_fn,
            "predicted_hump_fn": predicted_hump_fn,
            "hump_fn_error": hump_fn_error,
            "hump_position_grade": fn_grade,
            "source_hump_tau": source_hump_tau,
            "predicted_hump_tau": predicted_hump_tau,
            "hump_amplitude_ratio": hump_amplitude_ratio,
            "hump_amplitude_grade": amp_grade,
            "verdict": verdict([point_grade, fn_grade, amp_grade]),
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        handle.write(f"# criteria_commit: {CRITERIA_COMMIT}\n")
        handle.write(f"# digitization_commit: {DIGITIZATION_COMMIT}\n")
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_overlays(
    source: dict[float, dict[str, np.ndarray]],
    predictions: dict[str, dict[str, np.ndarray]],
) -> None:
    for configuration, (separation, figure) in CONFIGURATIONS.items():
        source_curve = source[separation]
        predicted_curve = predictions[configuration]
        fig, axis = plt.subplots(figsize=(8.0, 4.8))
        axis.fill_between(
            source_curve["fn"],
            source_curve["tau"] - source_curve["uncertainty"],
            source_curve["tau"] + source_curve["uncertainty"],
            color="#2b6cb0",
            alpha=0.16,
            linewidth=0,
            label="digitization uncertainty",
        )
        axis.plot(
            source_curve["fn"], source_curve["tau"], color="#2b6cb0",
            linewidth=1.8, label=f"Insel C2 theory (Figure {figure})",
        )
        axis.plot(
            predicted_curve["fn"], predicted_curve["tau"], color="#c53030",
            linewidth=1.6, linestyle="--", label="michell unbounded-water theory",
        )
        axis.axhline(1.0, color="0.45", linewidth=0.8, linestyle=":")
        axis.axvspan(0.20, 0.80, color="0.4", alpha=0.035, label="scoring range")
        axis.set(xlabel="Froude number", ylabel="wave-resistance interference ratio")
        axis.set_title(f"C2 theory comparison, S/L = {separation:.1f}")
        axis.set_xlim(0.14, 0.96)
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(PLOTS / f"theory-overlay-s-l-{separation:.1f}.png", dpi=180)
        plt.close(fig)


def main() -> None:
    ANALYSIS.mkdir(exist_ok=True)
    source = source_curves()
    predictions = prediction_curves()
    scores = score(source, predictions)
    write_csv(ANALYSIS / "theory_scores.csv", scores)
    plot_overlays(source, predictions)

    study_verdict = verdict([str(row["verdict"]) for row in scores])
    print(f"theory-to-theory study verdict: {study_verdict}")
    for row in scores:
        print(
            f"S/L={row['s_over_l']:.1f}: {row['verdict']}; "
            f"median={row['median_absolute_error']:.4f}, "
            f"p90={row['p90_absolute_error']:.4f}, "
            f"hump dFn={row['hump_fn_error']:.3f}, "
            f"amplitude ratio={row['hump_amplitude_ratio']:.3f}"
        )


if __name__ == "__main__":
    main()

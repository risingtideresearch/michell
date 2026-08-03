#!/usr/bin/env python3
"""Apply CRITERIA-SEPARATION.md and plot the registered separation fits."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STUDY = Path(__file__).resolve().parents[1]
SOURCE = STUDY / "data/digitized/theory_c2_attributed_359_362.csv"
PREDICTIONS = STUDY / "data/predictions/separation_grid_predictions.csv"
ANALYSIS = STUDY / "data/analysis"
PLOTS = Path(__file__).resolve().parent
CRITERIA_COMMIT = "f878a9a"
EPSILON = 1e-12

CONFIGURATIONS = {
    359: 0.2,
    360: 0.3,
    361: 0.4,
    362: 0.5,
}

TRANSFORMATIONS = {
    "gap_minus_beam": lambda label: label - 0.1,
    "half_spacing": lambda label: label / 2.0,
    "null_label": lambda label: label,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def source_curves() -> dict[int, dict[str, np.ndarray]]:
    grouped = {figure: [] for figure in CONFIGURATIONS}
    for row in read_csv(SOURCE):
        grouped[int(row["figure"])].append(row)
    result = {}
    for figure, rows in grouped.items():
        rows.sort(key=lambda row: float(row["fn_anchor"]))
        result[figure] = {
            "fn": np.array([float(row["fn_anchor"]) for row in rows]),
            "tau": np.array([float(row["tau"]) for row in rows]),
            "uncertainty": np.array([float(row["tau_digitization_uncertainty"]) for row in rows]),
        }
    return result


def prediction_curves() -> dict[float, dict[str, np.ndarray]]:
    grouped: dict[float, list[dict[str, str]]] = {}
    for row in read_csv(PREDICTIONS):
        assert row["pair_outcome"] == "converged"
        assert row["solo_outcome"] == "converged"
        grouped.setdefault(float(row["separation_over_length"]), []).append(row)
    result = {}
    for separation, rows in grouped.items():
        rows.sort(key=lambda row: float(row["fn"]))
        assert len(rows) == 171
        result[round(separation, 3)] = {
            "fn": np.array([float(row["fn"]) for row in rows]),
            "tau": np.array([float(row["interference"]) for row in rows]),
        }
    assert len(result) == 95
    return result


def interpolation(curve: dict[str, np.ndarray], query: np.ndarray) -> np.ndarray:
    return np.interp(query, curve["fn"], curve["tau"])


def objective(source: dict[str, np.ndarray], prediction: dict[str, np.ndarray]) -> float:
    predicted = interpolation(prediction, source["fn"])
    return float(np.median(np.abs(source["tau"] - predicted)))


def accuracy(source: dict[str, np.ndarray], prediction: dict[str, np.ndarray]) -> dict[str, object]:
    anchors = source["fn"]
    exact_hundredth = np.abs(anchors * 100.0 - np.round(anchors * 100.0)) <= 1e-7
    scoring = (anchors >= 0.20 - EPSILON) & (anchors <= 0.80 + EPSILON) & exact_hundredth
    query = anchors[scoring]
    source_values = source["tau"][scoring]
    predicted_values = interpolation(prediction, query)
    errors = np.abs(predicted_values - source_values)
    hump_grid = np.round(np.arange(0.35, 0.5501, 0.001), 3)
    source_hump = interpolation(source, hump_grid)
    predicted_hump = interpolation(prediction, hump_grid)
    source_index = int(np.argmax(source_hump))
    predicted_index = int(np.argmax(predicted_hump))
    source_peak = float(source_hump[source_index])
    predicted_peak = float(predicted_hump[predicted_index])
    count_pass = len(query) >= 50
    median_error = float(np.median(errors))
    p90_error = float(np.quantile(errors, 0.90))
    hump_difference = abs(float(hump_grid[predicted_index] - hump_grid[source_index]))
    peak_ratio = predicted_peak / source_peak
    return {
        "scoring_points": len(query),
        "median_absolute_error": median_error,
        "p90_absolute_error": p90_error,
        "source_hump_fn": float(hump_grid[source_index]),
        "source_hump_tau": source_peak,
        "prediction_hump_fn": float(hump_grid[predicted_index]),
        "prediction_hump_tau": predicted_peak,
        "absolute_hump_fn_error": hump_difference,
        "prediction_source_peak_ratio": peak_ratio,
        "count_pass": count_pass,
        "median_pass": median_error <= 0.05 + EPSILON,
        "p90_pass": p90_error <= 0.10 + EPSILON,
        "hump_position_pass": hump_difference <= 0.010 + EPSILON,
        "hump_amplitude_pass": 0.90 - EPSILON <= peak_ratio <= 1.10 + EPSILON,
        "agreement": (
            count_pass
            and median_error <= 0.05 + EPSILON
            and p90_error <= 0.10 + EPSILON
            and hump_difference <= 0.010 + EPSILON
            and 0.90 - EPSILON <= peak_ratio <= 1.10 + EPSILON
        ),
    }


def write(path: Path, rows: list[dict[str, object]], comments: tuple[str, ...]) -> None:
    ANALYSIS.mkdir(exist_ok=True)
    with path.open("w", newline="") as handle:
        for comment in comments:
            handle.write(f"# {comment}\n")
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    sources = source_curves()
    predictions = prediction_curves()
    separations = np.array(sorted(predictions))
    fit_rows: list[dict[str, object]] = []
    fit_intervals: dict[int, tuple[float, float]] = {}
    best_separations: dict[int, float] = {}
    objectives: dict[int, np.ndarray] = {}
    for figure, label in CONFIGURATIONS.items():
        source = sources[figure]
        values = np.array([objective(source, predictions[float(s)]) for s in separations])
        objectives[figure] = values
        best_index = int(np.argmin(values))
        best = float(separations[best_index])
        admitted = values <= values[best_index] + 0.010 + EPSILON
        lower = best_index
        upper = best_index
        while lower > 0 and admitted[lower - 1]:
            lower -= 1
        while upper + 1 < len(separations) and admitted[upper + 1]:
            upper += 1
        grid_lower = float(separations[lower])
        grid_upper = float(separations[upper])
        reported_lower = max(0.08, grid_lower - 0.0025)
        reported_upper = min(0.55, grid_upper + 0.0025)
        fit_intervals[figure] = (reported_lower, reported_upper)
        best_separations[figure] = best
        label_accuracy = accuracy(source, predictions[round(label, 3)])
        best_accuracy = accuracy(source, predictions[round(best, 3)])
        fit_rows.append({
            "figure": figure,
            "label_separation_over_length": f"{label:.3f}",
            "best_separation_over_length": f"{best:.3f}",
            "fit_grid_lower": f"{grid_lower:.3f}",
            "fit_grid_upper": f"{grid_upper:.3f}",
            "fit_reported_lower": f"{reported_lower:.4f}",
            "fit_reported_upper": f"{reported_upper:.4f}",
            "minimum_median_absolute_error": f"{values[best_index]:.7f}",
            "label_median_absolute_error": f"{label_accuracy['median_absolute_error']:.7f}",
            "best_median_absolute_error": f"{best_accuracy['median_absolute_error']:.7f}",
            "best_p90_absolute_error": f"{best_accuracy['p90_absolute_error']:.7f}",
            "best_absolute_hump_fn_error": f"{best_accuracy['absolute_hump_fn_error']:.7f}",
            "best_peak_ratio": f"{best_accuracy['prediction_source_peak_ratio']:.7f}",
            "best_agreement": str(best_accuracy["agreement"]).lower(),
        })

    transformation_rows: list[dict[str, object]] = []
    transformation_counts: dict[str, int] = {}
    for transformation, mapping in TRANSFORMATIONS.items():
        count = 0
        for figure, label in CONFIGURATIONS.items():
            expected = round(mapping(label), 3)
            lower, upper = fit_intervals[figure]
            interval_consistent = lower - EPSILON <= expected <= upper + EPSILON
            result = accuracy(sources[figure], predictions[expected])
            panel_consistent = interval_consistent and bool(result["agreement"])
            count += int(panel_consistent)
            transformation_rows.append({
                "transformation": transformation,
                "figure": figure,
                "label_separation_over_length": f"{label:.3f}",
                "predicted_separation_over_length": f"{expected:.3f}",
                "best_separation_over_length": f"{best_separations[figure]:.3f}",
                "fit_reported_lower": f"{lower:.4f}",
                "fit_reported_upper": f"{upper:.4f}",
                "fit_interval_consistent": str(interval_consistent).lower(),
                "scoring_points": result["scoring_points"],
                "median_absolute_error": f"{result['median_absolute_error']:.7f}",
                "p90_absolute_error": f"{result['p90_absolute_error']:.7f}",
                "source_hump_fn": f"{result['source_hump_fn']:.3f}",
                "prediction_hump_fn": f"{result['prediction_hump_fn']:.3f}",
                "absolute_hump_fn_error": f"{result['absolute_hump_fn_error']:.3f}",
                "prediction_source_peak_ratio": f"{result['prediction_source_peak_ratio']:.7f}",
                "count_pass": str(result["count_pass"]).lower(),
                "median_pass": str(result["median_pass"]).lower(),
                "p90_pass": str(result["p90_pass"]).lower(),
                "hump_position_pass": str(result["hump_position_pass"]).lower(),
                "hump_amplitude_pass": str(result["hump_amplitude_pass"]).lower(),
                "agreement": str(result["agreement"]).lower(),
                "panel_consistent": str(panel_consistent).lower(),
            })
        transformation_counts[transformation] = count
    winner, winner_count = max(transformation_counts.items(), key=lambda item: item[1])
    if winner_count == 4:
        outcome = "CONFIRMED"
    elif winner_count == 3:
        outcome = "PARTIAL"
    else:
        outcome = "REFUTED"
    for row in transformation_rows:
        row["transformation_consistent_panels"] = transformation_counts[str(row["transformation"])]
        row["registered_winner"] = winner
        row["registered_outcome"] = outcome

    comments = (
        f"criteria_commit: {CRITERIA_COMMIT}",
        "source: attributed Insel Figures 359--362",
        "prediction: exact C2 library solver, all rows converged",
    )
    write(ANALYSIS / "separation_fit_scores.csv", fit_rows, comments)
    write(ANALYSIS / "separation_transformation_scores.csv", transformation_rows, comments)

    for figure, label in CONFIGURATIONS.items():
        source = sources[figure]
        best = best_separations[figure]
        figure_object, axis = plt.subplots(figsize=(8.0, 4.8), constrained_layout=True)
        axis.fill_between(
            source["fn"], source["tau"] - source["uncertainty"],
            source["tau"] + source["uncertainty"], color="#2b6cb0", alpha=0.15,
            linewidth=0, label="digitization uncertainty",
        )
        axis.plot(source["fn"], source["tau"], color="#2b6cb0", linewidth=1.8,
                  label=f"Insel Figure {figure}, label S/L={label:.1f}")
        label_curve = predictions[round(label, 3)]
        best_curve = predictions[round(best, 3)]
        axis.plot(label_curve["fn"], label_curve["tau"], color="#c53030",
                  linestyle="--", linewidth=1.4, label=f"library at label {label:.3f}")
        axis.plot(best_curve["fn"], best_curve["tau"], color="#2f855a",
                  linewidth=1.5, label=f"library best fit {best:.3f}")
        axis.axhline(1.0, color="#555555", linewidth=0.8, linestyle=":")
        axis.axvspan(0.20, 0.80, color="#718096", alpha=0.035, label="scoring range")
        axis.set_xlim(0.15, 0.95)
        axis.set_ylim(0.55, 2.10)
        axis.set_xlabel("Froude number")
        axis.set_ylabel("wave-resistance interference ratio")
        axis.set_title(f"Separation-definition fit, Figure {figure}")
        axis.grid(alpha=0.22)
        axis.legend(loc="best", fontsize=8)
        figure_object.savefig(PLOTS / f"separation-fit-s-l-{label:.1f}.png", dpi=180)
        plt.close(figure_object)

    print(f"separation-definition outcome: {outcome}; winner={winner} ({winner_count}/4)")
    for row in fit_rows:
        print(
            f"Figure {row['figure']}: label={row['label_separation_over_length']}, "
            f"best={row['best_separation_over_length']} "
            f"[{row['fit_reported_lower']}, {row['fit_reported_upper']}], "
            f"median={row['best_median_absolute_error']}"
        )
    print(f"transformation counts: {transformation_counts}")


if __name__ == "__main__":
    main()

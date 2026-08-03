#!/usr/bin/env python3
"""Apply CRITERIA-CANAL.md and regenerate the finite-canal attribution study."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STUDY = Path(__file__).resolve().parents[1]
SOURCE = STUDY / "data/digitized/theory_interference.csv"
UNBOUNDED = STUDY / "data/predictions/theory_predictions.csv"
CANAL = STUDY / "data/predictions/canal_predictions.csv"
DIGITIZED = STUDY / "data/digitized"
ANALYSIS = STUDY / "data/analysis"
PLOTS = Path(__file__).resolve().parent
CRITERIA_COMMIT = "79e267d"
IMPLEMENTATION_COMMIT = "b0e70c4"
DIGITIZATION_COMMIT = "325aecf"
EPSILON = 1e-12

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
    grouped = {separation: [] for separation, _ in CONFIGURATIONS.values()}
    for row in read_csv(SOURCE):
        grouped[float(row["separation_over_length"])].append(row)
    result = {}
    for separation, rows in grouped.items():
        rows.sort(key=lambda row: float(row["fn_anchor"]))
        result[separation] = {
            "fn": np.array([float(row["fn_anchor"]) for row in rows]),
            "tau": np.array([float(row["tau"]) for row in rows]),
            "uncertainty": np.array(
                [float(row["tau_digitization_uncertainty"]) for row in rows]
            ),
        }
    return result


def prediction_curves(path: Path) -> dict[str, dict[str, np.ndarray]]:
    grouped = {configuration: [] for configuration in CONFIGURATIONS}
    for row in read_csv(path):
        if row["configuration"] in grouped:
            grouped[row["configuration"]].append(row)
    result = {}
    expected_grid = [round(0.20 + 0.005 * index, 3) for index in range(151)]
    for configuration, rows in grouped.items():
        rows.sort(key=lambda row: float(row["fn"]))
        assert len(rows) == 151, (path, configuration, len(rows))
        assert [float(row["fn"]) for row in rows] == expected_grid
        assert all(row["outcome"] == "converged" for row in rows)
        if "solo_outcome" in rows[0]:
            assert all(row["solo_outcome"] == "converged" for row in rows)
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


def experimental(configuration: str) -> dict[str, np.ndarray]:
    rows = [
        row
        for row in read_csv(DIGITIZED / f"fixed_{configuration}.csv")
        if row["observable"] == "cwp"
    ]
    return {
        "fn": np.array([float(row["fn"]) for row in rows]),
        "value": np.array([float(row["value"]) for row in rows]),
    }


def smooth(data: dict[str, np.ndarray], query: np.ndarray) -> np.ndarray:
    result = np.full(query.shape, np.nan, dtype=float)
    for index, q in enumerate(query):
        delta = data["fn"] - q
        selected = np.abs(delta) <= 0.0375 + EPSILON
        if np.count_nonzero(selected) < 2:
            continue
        weights = np.exp(-0.5 * (delta[selected] / 0.0125) ** 2)
        result[index] = np.sum(weights * data["value"][selected]) / np.sum(weights)
    return result


def experimental_tau(configuration: str, query: np.ndarray) -> np.ndarray:
    monohull = smooth(experimental("monohull"), query)
    catamaran = smooth(experimental(configuration), query)
    result = catamaran / monohull
    result[(~np.isfinite(monohull)) | (monohull < 0.0005)] = np.nan
    return result


def score(
    source: dict[float, dict[str, np.ndarray]],
    unbounded: dict[str, dict[str, np.ndarray]],
    canal: dict[str, dict[str, np.ndarray]],
) -> tuple[list[dict[str, object]], str]:
    hump_grid = np.round(np.arange(0.35, 0.5501, 0.001), 3)
    rows = []
    for configuration, (separation, figure) in CONFIGURATIONS.items():
        source_curve = source[separation]
        scoring = (source_curve["fn"] >= 0.20 - EPSILON) & (
            source_curve["fn"] <= 0.95 + EPSILON
        )
        query = source_curve["fn"][scoring]
        source_values = source_curve["tau"][scoring]
        canal_values = interpolate(canal[configuration], query)
        assert np.all(np.isfinite(canal_values))
        signed_error = canal_values - source_values
        absolute_error = np.abs(signed_error)

        source_hump = interpolate(source_curve, hump_grid)
        source_index = int(np.nanargmax(source_hump))
        source_hump_fn = float(hump_grid[source_index])
        source_hump_tau = float(source_hump[source_index])
        common_query = np.array([source_hump_fn])
        unbounded_at_hump = float(interpolate(unbounded[configuration], common_query)[0])
        canal_at_hump = float(interpolate(canal[configuration], common_query)[0])
        denominator = source_hump_tau - unbounded_at_hump
        assert abs(denominator) >= 0.02 - EPSILON, (configuration, denominator)
        gap_closure = (canal_at_hump - unbounded_at_hump) / denominator
        hump_relative_error = abs(canal_at_hump - source_hump_tau) / abs(source_hump_tau)

        unbounded_hump = interpolate(unbounded[configuration], hump_grid)
        canal_hump = interpolate(canal[configuration], hump_grid)
        unbounded_index = int(np.nanargmax(unbounded_hump))
        canal_index = int(np.nanargmax(canal_hump))
        rows.append(
            {
                "configuration": configuration,
                "s_over_l": separation,
                "source_figure": figure,
                "points": len(query),
                "median_absolute_error": float(np.median(absolute_error)),
                "p90_absolute_error": float(np.quantile(absolute_error, 0.90)),
                "rms_error": float(np.sqrt(np.mean(signed_error**2))),
                "maximum_absolute_error": float(np.max(absolute_error)),
                "median_signed_error": float(np.median(signed_error)),
                "source_hump_fn": source_hump_fn,
                "source_hump_tau": source_hump_tau,
                "unbounded_tau_at_source_hump": unbounded_at_hump,
                "canal_tau_at_source_hump": canal_at_hump,
                "gap_closure_fraction": gap_closure,
                "canal_hump_relative_error": hump_relative_error,
                "unbounded_own_hump_fn": float(hump_grid[unbounded_index]),
                "unbounded_own_hump_tau": float(unbounded_hump[unbounded_index]),
                "canal_own_hump_fn": float(hump_grid[canal_index]),
                "canal_own_hump_tau": float(canal_hump[canal_index]),
            }
        )

    reproduced = all(
        float(row["median_absolute_error"]) <= 0.02 + EPSILON
        and float(row["canal_hump_relative_error"]) <= 0.05 + EPSILON
        for row in rows
    )
    close_rows = [row for row in rows if float(row["s_over_l"]) in (0.2, 0.3)]
    if reproduced:
        outcome = "REPRODUCED"
    elif any(float(row["gap_closure_fraction"]) < 0.5 - EPSILON for row in close_rows):
        outcome = "UNEXPLAINED"
    else:
        outcome = "PARTIAL"
    for row in rows:
        row["study_outcome"] = outcome
    return rows, outcome


def write_scores(rows: list[dict[str, object]]) -> None:
    ANALYSIS.mkdir(exist_ok=True)
    path = ANALYSIS / "canal_scores.csv"
    with path.open("w", newline="") as handle:
        handle.write(f"# criteria_commit: {CRITERIA_COMMIT}\n")
        handle.write(f"# implementation_commit: {IMPLEMENTATION_COMMIT}\n")
        handle.write(f"# digitization_commit: {DIGITIZATION_COMMIT}\n")
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_overlays(
    source: dict[float, dict[str, np.ndarray]],
    unbounded: dict[str, dict[str, np.ndarray]],
    canal: dict[str, dict[str, np.ndarray]],
) -> None:
    for configuration, (separation, figure_number) in CONFIGURATIONS.items():
        source_curve = source[separation]
        figure, axis = plt.subplots(figsize=(8.0, 4.8), constrained_layout=True)
        axis.fill_between(
            source_curve["fn"],
            source_curve["tau"] - source_curve["uncertainty"],
            source_curve["tau"] + source_curve["uncertainty"],
            color="#2b6cb0",
            alpha=0.15,
            linewidth=0,
            label="digitization uncertainty",
        )
        axis.plot(
            source_curve["fn"],
            source_curve["tau"],
            color="#2b6cb0",
            linewidth=1.8,
            label=f"Insel canal theory (Figure {figure_number})",
        )
        axis.plot(
            unbounded[configuration]["fn"],
            unbounded[configuration]["tau"],
            color="#c53030",
            linestyle="--",
            linewidth=1.5,
            label="library unbounded theory",
        )
        axis.plot(
            canal[configuration]["fn"],
            canal[configuration]["tau"],
            color="#2f855a",
            linewidth=1.5,
            label="independent canal reference",
        )
        axis.axhline(1.0, color="0.45", linewidth=0.8, linestyle=":")
        axis.set(
            xlabel="Froude number",
            ylabel="wave-resistance interference ratio",
            title=f"C2 water-geometry attribution, S/L = {separation:.1f}",
            xlim=(0.14, 0.96),
        )
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8, ncol=2)
        figure.savefig(PLOTS / f"canal-theory-overlay-s-l-{separation:.1f}.png", dpi=180)
        plt.close(figure)

        query = np.round(np.arange(0.20, 0.8001, 0.001), 3)
        measured = experimental_tau(configuration, query)
        figure, axis = plt.subplots(figsize=(8.0, 4.8), constrained_layout=True)
        axis.plot(source_curve["fn"], source_curve["tau"], color="#2b6cb0", label="Insel canal theory")
        axis.plot(unbounded[configuration]["fn"], unbounded[configuration]["tau"], color="#c53030", linestyle="--", label="library unbounded theory")
        axis.plot(canal[configuration]["fn"], canal[configuration]["tau"], color="#2f855a", label="independent canal reference")
        axis.plot(query, measured, color="black", linewidth=1.4, label=r"measured $\tau_{WP}$")
        axis.axhline(1.0, color="0.45", linewidth=0.8, linestyle=":")
        axis.set(
            xlabel="Froude number",
            ylabel="wave-interference ratio",
            title=f"C2 theory and experiment, S/L = {separation:.1f}",
            xlim=(0.20, 0.80),
        )
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8, ncol=2)
        figure.savefig(PLOTS / f"canal-experiment-overlay-s-l-{separation:.1f}.png", dpi=180)
        plt.close(figure)


def main() -> None:
    source = source_curves()
    unbounded = prediction_curves(UNBOUNDED)
    canal = prediction_curves(CANAL)
    rows, outcome = score(source, unbounded, canal)
    write_scores(rows)
    plot_overlays(source, unbounded, canal)
    print(f"finite-canal attribution outcome: {outcome}")
    for row in rows:
        print(
            f"S/L={row['s_over_l']:.1f}: median={row['median_absolute_error']:.4f}, "
            f"hump error={row['canal_hump_relative_error']:.3%}, "
            f"gap closure={row['gap_closure_fraction']:.3f}"
        )


if __name__ == "__main__":
    main()

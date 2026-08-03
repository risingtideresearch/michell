#!/usr/bin/env python3
"""Apply CRITERIA.md verbatim and regenerate Insel-Wigley comparison plots."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STUDY = Path(__file__).resolve().parents[1]
DIGITIZED = STUDY / "data/digitized"
PREDICTIONS = STUDY / "data/predictions/predictions.csv"
ANALYSIS = STUDY / "data/analysis"
PLOTS = Path(__file__).resolve().parent
CRITERIA_COMMIT = "3f81c2a"
COMPARISON_EPSILON = 1e-12

CONFIGURATIONS = ["monohull", "s_l_0_2", "s_l_0_3", "s_l_0_4", "s_l_0_5"]
SEPARATIONS = ["s_l_0_2", "s_l_0_3", "s_l_0_4", "s_l_0_5"]
S_OVER_L = {
    "s_l_0_2": 0.2,
    "s_l_0_3": 0.3,
    "s_l_0_4": 0.4,
    "s_l_0_5": 0.5,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def experimental(attitude: str, configuration: str) -> dict[str, np.ndarray]:
    rows = [
        row
        for row in read_csv(DIGITIZED / f"{attitude}_{configuration}.csv")
        if row["observable"] == "cwp"
    ]
    return {
        "fn": np.array([float(row["fn"]) for row in rows]),
        "value": np.array([float(row["value"]) for row in rows]),
        "fn_uncertainty": np.array(
            [float(row["fn_digitization_uncertainty"]) for row in rows]
        ),
        "value_uncertainty": np.array(
            [float(row["value_digitization_uncertainty"]) for row in rows]
        ),
    }


def smooth(data: dict[str, np.ndarray], query: np.ndarray) -> np.ndarray:
    result = np.full(query.shape, np.nan, dtype=float)
    for index, q in enumerate(query):
        delta = data["fn"] - q
        selected = np.abs(delta) <= 0.0375 + 1e-12
        if np.count_nonzero(selected) < 2:
            continue
        weights = np.exp(-0.5 * (delta[selected] / 0.0125) ** 2)
        result[index] = np.sum(weights * data["value"][selected]) / np.sum(weights)
    return result


def prediction_groups() -> dict[str, dict[str, np.ndarray]]:
    grouped: dict[str, list[dict[str, str]]] = {name: [] for name in CONFIGURATIONS}
    for row in read_csv(PREDICTIONS):
        grouped[row["configuration"]].append(row)
    result = {}
    for configuration, rows in grouped.items():
        rows.sort(key=lambda row: float(row["fn"]))
        result[configuration] = {
            "fn": np.array([float(row["fn"]) for row in rows]),
            "cw": np.array([float(row["cw"]) for row in rows]),
            "interference": np.array([float(row["interference"]) for row in rows]),
            "est_rel_error": np.array([float(row["est_rel_error"]) for row in rows]),
            "converged": np.array(
                [row["outcome"] == "converged" and row["solo_outcome"] == "converged" for row in rows]
            ),
        }
    return result


def interpolate(data: dict[str, np.ndarray], field: str, query: np.ndarray) -> np.ndarray:
    result = np.interp(query, data["fn"], data[field])
    outside = (query < data["fn"][0]) | (query > data["fn"][-1])
    result[outside] = np.nan
    return result


def quantile(values: np.ndarray, probability: float) -> float:
    finite = values[np.isfinite(values)]
    return float(np.quantile(finite, probability)) if len(finite) else math.nan


def threshold_grade(
    first: float,
    second: float,
    agreement: tuple[float, float],
    partial: tuple[float, float],
) -> str:
    if first <= agreement[0] + COMPARISON_EPSILON and second <= agreement[1] + COMPARISON_EPSILON:
        return "agreement"
    if first <= partial[0] + COMPARISON_EPSILON and second <= partial[1] + COMPARISON_EPSILON:
        return "partial"
    return "disagreement"


def upper_grade(value: float, agreement: float, partial: float) -> str:
    if value <= agreement + COMPARISON_EPSILON:
        return "agreement"
    if value <= partial + COMPARISON_EPSILON:
        return "partial"
    return "disagreement"


def lower_grade(value: float, agreement: float, partial: float) -> str:
    if value >= agreement - COMPARISON_EPSILON:
        return "agreement"
    if value >= partial - COMPARISON_EPSILON:
        return "partial"
    return "disagreement"


def verdict(grades: list[str]) -> str:
    applicable = [grade for grade in grades if grade != "not_applicable"]
    if applicable and all(grade == "agreement" for grade in applicable):
        return "agreement"
    if any(grade == "disagreement" for grade in applicable):
        return "disagreement"
    return "partial" if applicable else "not_applicable"


def coefficient_scores(
    predictions: dict[str, dict[str, np.ndarray]],
    fixed: dict[str, dict[str, np.ndarray]],
) -> list[dict[str, object]]:
    scoring = np.arange(0.25, 0.8001, 0.025)
    dense = np.arange(0.38, 0.5801, 0.001)
    rows = []
    for configuration in CONFIGURATIONS:
        experiment = smooth(fixed[configuration], scoring)
        predicted = interpolate(predictions[configuration], "cw", scoring)
        valid = np.isfinite(experiment) & np.isfinite(predicted)
        discrepancy = np.abs(predicted[valid] - experiment[valid]) / np.maximum(
            np.abs(experiment[valid]), 0.0005
        )
        median = quantile(discrepancy, 0.5)
        p80 = quantile(discrepancy, 0.8)
        point_grade = threshold_grade(median, p80, (0.25, 0.45), (0.40, 0.75))

        exp_dense = smooth(fixed[configuration], dense)
        pred_dense = interpolate(predictions[configuration], "cw", dense)
        exp_valid = np.isfinite(exp_dense)
        pred_valid = np.isfinite(pred_dense)
        exp_index = int(np.nanargmax(np.where(exp_valid, exp_dense, np.nan)))
        pred_index = int(np.nanargmax(np.where(pred_valid, pred_dense, np.nan)))
        exp_fn = float(dense[exp_index])
        pred_fn = float(dense[pred_index])
        exp_amplitude = float(exp_dense[exp_index])
        pred_amplitude = float(pred_dense[pred_index])
        fn_error = abs(pred_fn - exp_fn)
        amplitude_error = abs(pred_amplitude - exp_amplitude) / max(
            abs(exp_amplitude), 0.0005
        )
        fn_grade = upper_grade(fn_error, 0.020, 0.040)
        amplitude_grade = upper_grade(amplitude_error, 0.25, 0.45)
        rows.append(
            {
                "configuration": configuration,
                "points": int(np.count_nonzero(valid)),
                "median_d": median,
                "p80_d": p80,
                "pointwise_grade": point_grade,
                "experimental_hump_fn": exp_fn,
                "predicted_hump_fn": pred_fn,
                "hump_fn_error": fn_error,
                "hump_position_grade": fn_grade,
                "experimental_hump_cwp": exp_amplitude,
                "predicted_hump_cw": pred_amplitude,
                "hump_relative_amplitude_error": amplitude_error,
                "hump_amplitude_grade": amplitude_grade,
                "verdict": verdict([point_grade, fn_grade, amplitude_grade]),
            }
        )
    return rows


def tau_curves(
    predictions: dict[str, dict[str, np.ndarray]],
    fixed: dict[str, dict[str, np.ndarray]],
    query: np.ndarray,
    configuration: str,
) -> tuple[np.ndarray, np.ndarray]:
    mono = smooth(fixed["monohull"], query)
    cat = smooth(fixed[configuration], query)
    experimental_tau = cat / mono
    experimental_tau[(~np.isfinite(mono)) | (mono < 0.0005)] = np.nan
    predicted_tau = interpolate(predictions[configuration], "interference", query)
    return experimental_tau, predicted_tau


def feature_scores(
    query: np.ndarray,
    experimental_tau: np.ndarray,
    predicted_tau: np.ndarray,
    kind: str,
) -> dict[str, object]:
    if kind == "hump":
        window = (query >= 0.38) & (query <= 0.52)
        exp_index = int(np.nanargmax(np.where(window, experimental_tau, np.nan)))
        pred_index = int(np.nanargmax(np.where(window, predicted_tau, np.nan)))
        exp_prominence = experimental_tau[exp_index] - 1.0
        predicted_exists = predicted_tau[pred_index] > 1.05
    else:
        window = (query >= 0.30) & (query <= 0.42)
        exp_index = int(np.nanargmin(np.where(window, experimental_tau, np.nan)))
        pred_index = int(np.nanargmin(np.where(window, predicted_tau, np.nan)))
        exp_prominence = 1.0 - experimental_tau[exp_index]
        predicted_exists = predicted_tau[pred_index] < 0.95
    if exp_prominence < 0.10:
        return {
            f"{kind}_experimental_fn": float(query[exp_index]),
            f"{kind}_predicted_fn": float(query[pred_index]),
            f"{kind}_fn_error": math.nan,
            f"{kind}_position_grade": "not_applicable",
            f"{kind}_experimental_ratio": float(experimental_tau[exp_index]),
            f"{kind}_predicted_ratio": float(predicted_tau[pred_index]),
            f"{kind}_amplitude_error": math.nan,
            f"{kind}_amplitude_grade": "not_applicable",
        }
    fn_error = abs(float(query[pred_index] - query[exp_index]))
    amplitude_error = abs(float(predicted_tau[pred_index] - experimental_tau[exp_index]))
    if predicted_exists:
        position_grade = upper_grade(fn_error, 0.020, 0.040)
        amplitude_grade = upper_grade(amplitude_error, 0.25, 0.50)
    else:
        position_grade = "disagreement"
        amplitude_grade = "disagreement"
    return {
        f"{kind}_experimental_fn": float(query[exp_index]),
        f"{kind}_predicted_fn": float(query[pred_index]),
        f"{kind}_fn_error": fn_error,
        f"{kind}_position_grade": position_grade,
        f"{kind}_experimental_ratio": float(experimental_tau[exp_index]),
        f"{kind}_predicted_ratio": float(predicted_tau[pred_index]),
        f"{kind}_amplitude_error": amplitude_error,
        f"{kind}_amplitude_grade": amplitude_grade,
    }


def interference_scores(
    predictions: dict[str, dict[str, np.ndarray]],
    fixed: dict[str, dict[str, np.ndarray]],
) -> list[dict[str, object]]:
    scoring = np.arange(0.25, 0.5501, 0.025)
    dense = np.arange(0.25, 0.5501, 0.001)
    rows = []
    for configuration in SEPARATIONS:
        exp_tau, pred_tau = tau_curves(predictions, fixed, scoring, configuration)
        valid = np.isfinite(exp_tau) & np.isfinite(pred_tau)
        level_error = np.abs(pred_tau[valid] - exp_tau[valid])
        median = quantile(level_error, 0.5)
        p80 = quantile(level_error, 0.8)
        level_grade = threshold_grade(median, p80, (0.25, 0.45), (0.40, 0.70))

        nonneutral = valid & (np.abs(exp_tau - 1.0) > 0.10)
        exp_sign = np.sign(exp_tau[nonneutral] - 1.0)
        pred_sign = np.where(
            pred_tau[nonneutral] > 1.05,
            1.0,
            np.where(pred_tau[nonneutral] < 0.95, -1.0, 0.0),
        )
        sign_fraction = (
            float(np.mean(exp_sign == pred_sign)) if len(exp_sign) else math.nan
        )
        sign_grade = (
            lower_grade(sign_fraction, 0.80, 0.60)
            if math.isfinite(sign_fraction)
            else "not_applicable"
        )

        exp_dense, pred_dense = tau_curves(predictions, fixed, dense, configuration)
        hump = feature_scores(dense, exp_dense, pred_dense, "hump")
        hollow = feature_scores(dense, exp_dense, pred_dense, "hollow")
        grades = [
            level_grade,
            sign_grade,
            str(hump["hump_position_grade"]),
            str(hump["hump_amplitude_grade"]),
            str(hollow["hollow_position_grade"]),
            str(hollow["hollow_amplitude_grade"]),
        ]
        rows.append(
            {
                "configuration": configuration,
                "s_over_l": S_OVER_L[configuration],
                "points": int(np.count_nonzero(valid)),
                "median_absolute_error": median,
                "p80_absolute_error": p80,
                "level_grade": level_grade,
                "nonneutral_points": len(exp_sign),
                "sign_agreement_fraction": sign_fraction,
                "sign_grade": sign_grade,
                **hump,
                **hollow,
                "verdict": verdict(grades),
            }
        )
    return rows


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2 + 1
        start = end
    return ranks


def spearman(first: np.ndarray, second: np.ndarray) -> float:
    first_ranks = average_ranks(first)
    second_ranks = average_ranks(second)
    return float(np.corrcoef(first_ranks, second_ranks)[0, 1])


def trend_scores(
    predictions: dict[str, dict[str, np.ndarray]],
    fixed: dict[str, dict[str, np.ndarray]],
) -> list[dict[str, object]]:
    scoring = np.arange(0.25, 0.5501, 0.025)
    exp_by_separation = []
    pred_by_separation = []
    for configuration in SEPARATIONS:
        experimental_tau, predicted_tau = tau_curves(
            predictions, fixed, scoring, configuration
        )
        exp_by_separation.append(experimental_tau)
        pred_by_separation.append(predicted_tau)
    exp_matrix = np.array(exp_by_separation)
    pred_matrix = np.array(pred_by_separation)
    correlations = []
    for index in range(len(scoring)):
        if np.all(np.isfinite(exp_matrix[:, index])) and np.all(
            np.isfinite(pred_matrix[:, index])
        ):
            correlations.append(spearman(exp_matrix[:, index], pred_matrix[:, index]))
    median_correlation = float(np.median(correlations))

    envelope_grid = np.arange(0.40, 0.5501, 0.025)
    exp_amplitudes = []
    pred_amplitudes = []
    for configuration in SEPARATIONS:
        experimental_tau, predicted_tau = tau_curves(
            predictions, fixed, envelope_grid, configuration
        )
        valid = np.isfinite(experimental_tau) & np.isfinite(predicted_tau)
        exp_amplitudes.append(float(np.sqrt(np.mean((experimental_tau[valid] - 1.0) ** 2))))
        pred_amplitudes.append(float(np.sqrt(np.mean((predicted_tau[valid] - 1.0) ** 2))))
    envelope_correlation = spearman(
        np.array(exp_amplitudes), np.array(pred_amplitudes)
    )
    return [
        {
            "metric": "pointwise_separation_ranking",
            "points": len(correlations),
            "spearman": median_correlation,
            "grade": lower_grade(median_correlation, 0.70, 0.40),
        },
        {
            "metric": "rms_interference_envelope_ranking",
            "points": 4,
            "spearman": envelope_correlation,
            "grade": lower_grade(envelope_correlation, 0.70, 0.40),
        },
    ]


def trust_envelope(
    predictions: dict[str, dict[str, np.ndarray]],
    fixed: dict[str, dict[str, np.ndarray]],
) -> list[dict[str, object]]:
    scoring = np.arange(0.25, 0.5501, 0.025)
    mono_exp = smooth(fixed["monohull"], scoring)
    mono_pred = interpolate(predictions["monohull"], "cw", scoring)
    mono_d = np.abs(mono_pred - mono_exp) / np.maximum(np.abs(mono_exp), 0.0005)
    rows = []
    for configuration in SEPARATIONS:
        cat_exp = smooth(fixed[configuration], scoring)
        cat_pred = interpolate(predictions[configuration], "cw", scoring)
        cat_d = np.abs(cat_pred - cat_exp) / np.maximum(np.abs(cat_exp), 0.0005)
        exp_tau, pred_tau = tau_curves(predictions, fixed, scoring, configuration)
        for index, fn_ in enumerate(scoring):
            valid = all(
                math.isfinite(value)
                for value in (mono_d[index], cat_d[index], exp_tau[index], pred_tau[index])
            )
            neutral = valid and abs(exp_tau[index] - 1.0) <= 0.10
            sign_match = valid and (
                neutral
                or (
                    exp_tau[index] > 1.0
                    and pred_tau[index] > 1.05
                )
                or (
                    exp_tau[index] < 1.0
                    and pred_tau[index] < 0.95
                )
            )
            pair_error = abs(pred_tau[index] - exp_tau[index]) if valid else math.nan
            passed = (
                valid
                and mono_d[index] <= 0.35
                and cat_d[index] <= 0.35
                and pair_error <= 0.35
                and sign_match
            )
            rows.append(
                {
                    "configuration": configuration,
                    "s_over_l": S_OVER_L[configuration],
                    "fn": fn_,
                    "monohull_d": mono_d[index],
                    "catamaran_d": cat_d[index],
                    "experimental_tau_wp": exp_tau[index],
                    "predicted_interference": pred_tau[index],
                    "interference_absolute_error": pair_error,
                    "experimental_neutral": neutral,
                    "sign_match": sign_match,
                    "passed": passed,
                }
            )
    return rows


def free_summary(
    predictions: dict[str, dict[str, np.ndarray]],
    free: dict[str, dict[str, np.ndarray]],
) -> list[dict[str, object]]:
    scoring = np.arange(0.25, 0.8001, 0.025)
    rows = []
    for configuration in CONFIGURATIONS:
        experiment = smooth(free[configuration], scoring)
        predicted = interpolate(predictions[configuration], "cw", scoring)
        discrepancy = np.abs(predicted - experiment) / np.maximum(np.abs(experiment), 0.0005)
        valid = np.isfinite(discrepancy)
        below = discrepancy[valid & (scoring <= 0.4)]
        above = discrepancy[valid & (scoring > 0.4)]
        rows.append(
            {
                "configuration": configuration,
                "below_or_at_fn_0_4_points": len(below),
                "below_or_at_fn_0_4_median_d": quantile(below, 0.5),
                "above_fn_0_4_points": len(above),
                "above_fn_0_4_median_d": quantile(above, 0.5),
            }
        )
    return rows


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        handle.write(f"# criteria_commit: {CRITERIA_COMMIT}\n")
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def plot_coefficients(
    predictions: dict[str, dict[str, np.ndarray]],
    observations: dict[str, dict[str, np.ndarray]],
    attitude: str,
) -> None:
    for configuration in CONFIGURATIONS:
        predicted = predictions[configuration]
        observed = observations[configuration]
        figure, axis = plt.subplots(figsize=(7.2, 4.3), constrained_layout=True)
        band = predicted["cw"] * predicted["est_rel_error"]
        axis.fill_between(
            predicted["fn"],
            predicted["cw"] - band,
            predicted["cw"] + band,
            color="C0",
            alpha=0.25,
            label="reported numerical band",
        )
        axis.plot(predicted["fn"], predicted["cw"], color="C0", label="Michell $C_W$")
        axis.errorbar(
            observed["fn"],
            observed["value"],
            xerr=observed["fn_uncertainty"],
            yerr=observed["value_uncertainty"],
            fmt="x",
            color="black",
            markersize=4,
            linewidth=0.6,
            label="Insel measured $C_{WP}$",
        )
        axis.set_xlabel(r"$F_n$")
        axis.set_ylabel("wave coefficient")
        axis.set_title(f"{attitude} {configuration.replace('_', ' ')}")
        axis.grid(alpha=0.2)
        axis.legend()
        figure.savefig(PLOTS / f"overlay-{attitude}-{configuration.replace('_', '-')}.png", dpi=180)
        plt.close(figure)


def plot_interference(
    predictions: dict[str, dict[str, np.ndarray]],
    fixed: dict[str, dict[str, np.ndarray]],
) -> None:
    query = np.arange(0.25, 0.5501, 0.001)
    for configuration in SEPARATIONS:
        experimental_tau, predicted_tau = tau_curves(
            predictions, fixed, query, configuration
        )
        figure, axis = plt.subplots(figsize=(7.2, 4.3), constrained_layout=True)
        axis.axhspan(0.9, 1.1, color="0.9", label="experimental neutral band")
        axis.axhline(1.0, color="0.4", linewidth=0.8)
        axis.plot(query, experimental_tau, color="black", label=r"measured $\tau_{WP}$")
        axis.plot(query, predicted_tau, color="C1", label="library interference")
        axis.set_xlabel(r"$F_n$")
        axis.set_ylabel("wave-interference ratio")
        axis.set_title(f"fixed {configuration.replace('_', ' ')}")
        axis.grid(alpha=0.2)
        axis.legend()
        figure.savefig(PLOTS / f"interference-{configuration.replace('_', '-')}.png", dpi=180)
        plt.close(figure)


def main() -> None:
    predictions = prediction_groups()
    assert all(np.all(values["converged"]) for values in predictions.values())
    fixed = {configuration: experimental("fixed", configuration) for configuration in CONFIGURATIONS}
    free = {configuration: experimental("free", configuration) for configuration in CONFIGURATIONS}

    coefficient = coefficient_scores(predictions, fixed)
    interference = interference_scores(predictions, fixed)
    trends = trend_scores(predictions, fixed)
    trust = trust_envelope(predictions, fixed)
    secondary = free_summary(predictions, free)
    write_table(ANALYSIS / "coefficient_scores.csv", coefficient)
    write_table(ANALYSIS / "interference_scores.csv", interference)
    write_table(ANALYSIS / "trend_scores.csv", trends)
    write_table(ANALYSIS / "trust_envelope.csv", trust)
    write_table(ANALYSIS / "free_attitude_summary.csv", secondary)
    plot_coefficients(predictions, fixed, "fixed")
    plot_coefficients(predictions, free, "free")
    plot_interference(predictions, fixed)
    print(
        f"wrote {len(coefficient)} coefficient scores, {len(interference)} interference scores, "
        f"{len(trust)} trust cells, 14 comparison plots; criteria {CRITERIA_COMMIT}"
    )


if __name__ == "__main__":
    main()

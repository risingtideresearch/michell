#!/usr/bin/env python3
"""Apply the preregistered critical-Froude consistency definition."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


STUDY = Path(__file__).resolve().parents[1]
PREDICTIONS = STUDY / "data/predictions/theory_predictions.csv"
OUTPUT = STUDY / "data/analysis/critical_froude.csv"
CRITERIA_COMMIT = "7f588d8"

CONFIGURATIONS = {
    "s_l_0_2": 0.2,
    "s_l_0_3": 0.3,
    "s_l_0_4": 0.4,
    "s_l_0_5": 0.5,
}


def read_predictions() -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with PREDICTIONS.open() as handle:
        rows = csv.DictReader(line for line in handle if not line.startswith("#"))
        for row in rows:
            grouped[row["configuration"]].append(row)
    assert set(grouped) == set(CONFIGURATIONS)
    for configuration, rows in grouped.items():
        rows.sort(key=lambda row: float(row["fn"]))
        assert len(rows) == 151, (configuration, len(rows))
        assert [float(row["fn"]) for row in rows] == [
            round(0.20 + 0.005 * index, 3) for index in range(151)
        ]
        assert all(
            row["outcome"] == "converged" and row["solo_outcome"] == "converged"
            for row in rows
        )
    return grouped


def main() -> None:
    predictions = read_predictions()
    source_statements = {"s_l_0_2": "about 0.8", "s_l_0_5": "about 0.55"}
    output_rows = []
    for configuration, separation in CONFIGURATIONS.items():
        rows = predictions[configuration]
        within = [abs(float(row["interference"]) - 1.0) < 0.05 for row in rows]
        critical_index = next(
            (index for index in range(len(within)) if all(within[index:])), None
        )
        output_rows.append({
            "configuration": configuration,
            "s_over_l": separation,
            "critical_fn": (
                "not_reached"
                if critical_index is None
                else f"{float(rows[critical_index]['fn']):.3f}"
            ),
            "insel_stated_fn": source_statements.get(configuration, ""),
            "definition": "first 0.005-grid Fn with abs(tau-1)<0.05 through Fn=0.95",
            "source_location": "Insel printed page 131; PDF page 141",
        })

    with OUTPUT.open("w", newline="") as handle:
        handle.write(f"# criteria_commit: {CRITERIA_COMMIT}\n")
        writer = csv.DictWriter(
            handle, fieldnames=output_rows[0].keys(), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print("critical Fn: " + ", ".join(
        f"S/L={row['s_over_l']:.1f}: {row['critical_fn']}" for row in output_rows
    ))


if __name__ == "__main__":
    main()

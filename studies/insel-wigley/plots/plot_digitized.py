#!/usr/bin/env python3
"""Plot reconciled Insel C2 marker data without reading raw PDFs or passes."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


STUDY = Path(__file__).resolve().parents[1]
DATA = STUDY / "data/digitized"
OUTPUT = Path(__file__).resolve().parent

LABELS = {
    "ct": r"$C_T$",
    "cwp": r"$C_{WP}$",
    "trim": "trim (degree)",
    "sinkage_over_draught": "sinkage / draught",
}


def read(path: Path) -> dict[str, list[tuple[float, float]]]:
    result: dict[str, list[tuple[float, float]]] = defaultdict(list)
    with path.open() as handle:
        rows = csv.DictReader(line for line in handle if not line.startswith("#"))
        for row in rows:
            result[row["observable"]].append((float(row["fn"]), float(row["value"])))
    return result


def main() -> None:
    for path in sorted(DATA.glob("*.csv")):
        if not path.name.startswith(("fixed_", "free_")):
            continue
        series = read(path)
        observables = [name for name in ("ct", "cwp", "trim", "sinkage_over_draught") if name in series]
        figure, axes = plt.subplots(
            len(observables),
            1,
            figsize=(8, 2.7 * len(observables)),
            sharex=True,
            constrained_layout=True,
        )
        if len(observables) == 1:
            axes = [axes]
        for axis, observable in zip(axes, observables, strict=True):
            points = series[observable]
            axis.scatter(
                [point[0] for point in points],
                [point[1] for point in points],
                marker="s" if observable in ("ct", "trim") else "x",
                facecolors="none" if observable in ("ct", "trim") else None,
                color="black",
                s=22,
                linewidths=0.9,
            )
            axis.set_ylabel(LABELS[observable])
            axis.grid(alpha=0.2)
        axes[-1].set_xlabel(r"$F_n$")
        figure.suptitle(path.stem.replace("_", " "))
        destination = OUTPUT / f"digitized-{path.stem.replace('_', '-')}.png"
        figure.savefig(destination, dpi=180)
        plt.close(figure)


if __name__ == "__main__":
    main()

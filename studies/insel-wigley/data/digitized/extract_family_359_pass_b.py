#!/usr/bin/env python3
"""Independent source-only trace of every curve in Insel Figure 359.

This pass uses a separately read calibration and checkpoint set.  Unlike pass
A's compact local search, it scores a wider horizontal strip to bridge printed
dash and dot gaps, then uses a vertical ink centroid.  It reads no solver data
and no output from pass A.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / "tmp/pdfs/insel-h1-family/theory-highres/page-368.jpg"
OVERLAY = ROOT / "tmp/pdfs/insel-h1-family/figure-359-pass-b.png"
OUTPUT = HERE / "passes/h1_figure_359_pass_b.tsv"


@dataclass(frozen=True)
class Curve:
    model: str
    hull: str
    line_style: str
    relative_position: str
    fn_first: float
    fn_last: float
    checkpoints: tuple[tuple[float, float], ...]


CURVES = (
    Curve(
        "C2", "Wigley hull", "solid",
        "oscillatory left branch and uniquely high broad crest",
        0.15, 0.80,
        (
            (0.15, 1.19), (0.18, 1.13), (0.20, 1.34), (0.23, 0.88),
            (0.25, 1.24), (0.28, 0.85), (0.30, 1.49), (0.33, 0.97),
            (0.35, 0.76), (0.38, 1.36), (0.40, 1.79), (0.43, 1.97),
            (0.45, 1.98), (0.48, 1.88), (0.50, 1.75), (0.55, 1.49),
            (0.60, 1.29), (0.65, 1.13), (0.70, 0.99), (0.75, 0.89),
            (0.80, 0.80),
        ),
    ),
    Curve(
        "C3", "round-bilge hull, L/B 7", "long dashed",
        "highest RBH line on the left shoulder; depressed middle shoulder; upper RBH crest",
        0.20, 0.64,
        (
            (0.20, 1.07), (0.25, 1.17), (0.30, 1.22), (0.34, 1.38),
            (0.37, 1.43), (0.40, 1.38), (0.43, 1.32), (0.45, 1.31),
            (0.48, 1.41), (0.50, 1.53), (0.53, 1.58), (0.55, 1.56),
            (0.58, 1.50), (0.60, 1.43), (0.64, 1.29),
        ),
    ),
    Curve(
        "C4", "round-bilge hull, L/B 9", "dash-dot",
        "near-horizontal lower RBH start; central rising branch; longest descending non-dotted branch",
        0.20, 0.84,
        (
            (0.20, 1.03), (0.25, 1.04), (0.30, 1.06), (0.35, 1.19),
            (0.40, 1.32), (0.45, 1.37), (0.48, 1.46), (0.50, 1.54),
            (0.53, 1.57), (0.55, 1.53), (0.58, 1.46), (0.60, 1.39),
            (0.65, 1.25), (0.70, 1.11), (0.75, 0.98), (0.80, 0.86),
            (0.84, 0.79),
        ),
    ),
    Curve(
        "C5", "round-bilge hull, L/B 11", "dotted",
        "bottom RBH line through the broad hump; separated shallow trough and upturn at high Fn",
        0.20, 1.00,
        (
            (0.20, 1.00), (0.25, 1.01), (0.30, 1.07), (0.35, 1.12),
            (0.40, 1.22), (0.45, 1.39), (0.48, 1.49), (0.50, 1.51),
            (0.53, 1.46), (0.55, 1.40), (0.60, 1.27), (0.65, 1.12),
            (0.70, 0.98), (0.75, 0.91), (0.80, 0.89), (0.85, 0.89),
            (0.90, 0.91), (0.95, 0.96), (1.00, 1.08),
        ),
    ),
)

X_LEFT = 536
X_RIGHT = 2402
Y_TOP = 315
Y_BOTTOM = 1581


def pixel_x(fn: float) -> float:
    return X_LEFT + (fn - 0.1) / 0.9 * (X_RIGHT - X_LEFT)


def pixel_y(tau: float) -> float:
    return Y_BOTTOM - tau / 2.5 * (Y_BOTTOM - Y_TOP)


def expected_tau(curve: Curve, fn: float) -> float:
    return float(np.interp(
        fn,
        [point[0] for point in curve.checkpoints],
        [point[1] for point in curve.checkpoints],
    ))


def locate_line(gray: np.ndarray, x: float, expected_y: float) -> float:
    xi = round(x)
    centre = round(expected_y)
    candidates = range(max(0, centre - 18), min(gray.shape[0], centre + 19))
    best_y = centre
    best_score = -float("inf")
    for y in candidates:
        patch = gray[max(0, y - 4):y + 5, max(0, xi - 16):xi + 17]
        darkness = np.clip(184.0 - patch.astype(float), 0.0, None)
        vertical_response = darkness.sum(axis=1)
        score = float(vertical_response.sum() + 0.5 * vertical_response.max())
        score -= 6.0 * (y - expected_y) ** 2
        if score > best_score:
            best_score = score
            best_y = y
    strip = gray[max(0, best_y - 5):best_y + 6, max(0, xi - 8):xi + 9]
    weights = np.clip(202.0 - strip.astype(float), 0.0, None).sum(axis=1)
    pixels = np.arange(max(0, best_y - 5), best_y + 6, dtype=float)
    return float((weights * pixels).sum() / weights.sum()) if weights.sum() else float(best_y)


def main() -> None:
    image = Image.open(SOURCE).convert("RGB")
    gray = np.asarray(image.convert("L"))
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    colours = {"C2": (220, 0, 0), "C3": (0, 100, 230), "C4": (0, 160, 70), "C5": (190, 0, 190)}
    rows: list[dict[str, str | int]] = []
    for curve in CURVES:
        count = round((curve.fn_last - curve.fn_first) / 0.005)
        for index in range(count + 1):
            fn = round(curve.fn_first + 0.005 * index, 3)
            x = pixel_x(fn)
            y = locate_line(gray, x, pixel_y(expected_tau(curve, fn)))
            rows.append({
                "pass": "B",
                "figure": 359,
                "pdf_page": 368,
                "printed_page": 358,
                "separation_over_length": "0.2",
                "model": curve.model,
                "hull": curve.hull,
                "line_style": curve.line_style,
                "relative_position_description": curve.relative_position,
                "fn_anchor": f"{fn:.3f}",
                "pixel_x": f"{x:.3f}",
                "pixel_y": f"{y:.3f}",
                "pixel_x_left": X_LEFT,
                "pixel_x_right": X_RIGHT,
                "pixel_y_top": Y_TOP,
                "pixel_y_bottom": Y_BOTTOM,
            })
            draw.rectangle((x - 2, y - 2, x + 2, y + 2), outline=colours[curve.model], width=1)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    overlay.save(OVERLAY)
    print(f"wrote {len(rows)} source-only points to {OUTPUT}")


if __name__ == "__main__":
    main()

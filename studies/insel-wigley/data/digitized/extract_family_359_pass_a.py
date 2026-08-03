#!/usr/bin/env python3
"""First source-only trace of every theoretical curve in Insel Figure 359.

Curve identities come only from the printed legend and relative positions in
the source panel.  The hand-read guide ordinates keep the pixel search on the
intended line; a local darkness score centres each common-Fn sample on ink.
The optional overlay is solely a source-identification check.
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
OVERLAY = ROOT / "tmp/pdfs/insel-h1-family/figure-359-pass-a.png"
OUTPUT = HERE / "passes/h1_figure_359_pass_a.tsv"


@dataclass(frozen=True)
class Curve:
    model: str
    hull: str
    line_style: str
    relative_position: str
    fn_first: float
    fn_last: float
    guides: tuple[tuple[float, float], ...]


CURVES = (
    Curve(
        "C2", "Wigley hull", "solid",
        "oscillatory below Fn 0.39; isolated upper broad hump near Fn 0.45",
        0.15, 0.80,
        (
            (0.15, 1.20), (0.175, 1.10), (0.20, 1.35), (0.225, 0.82),
            (0.25, 1.25), (0.275, 0.80), (0.30, 1.50), (0.325, 1.05),
            (0.35, 0.75), (0.375, 1.25), (0.40, 1.78), (0.425, 1.96),
            (0.45, 1.98), (0.475, 1.90), (0.50, 1.76), (0.525, 1.62),
            (0.55, 1.50), (0.60, 1.30), (0.65, 1.14), (0.70, 1.00),
            (0.75, 0.90), (0.80, 0.80),
        ),
    ),
    Curve(
        "C3", "round-bilge hull, L/B 7", "long dashed",
        "upper RBH branch before Fn 0.40; lower shoulder near Fn 0.45; upper main-hump shoulder",
        0.20, 0.64,
        (
            (0.20, 1.08), (0.25, 1.18), (0.30, 1.22), (0.35, 1.40),
            (0.375, 1.43), (0.40, 1.38), (0.425, 1.33), (0.45, 1.31),
            (0.475, 1.40), (0.50, 1.54), (0.525, 1.58), (0.55, 1.56),
            (0.575, 1.51), (0.60, 1.43), (0.625, 1.34), (0.64, 1.29),
        ),
    ),
    Curve(
        "C4", "round-bilge hull, L/B 9", "dash-dot",
        "near-unity initial RBH branch; middle rising branch; long descending branch after the common hump",
        0.20, 0.84,
        (
            (0.20, 1.03), (0.25, 1.05), (0.30, 1.06), (0.35, 1.20),
            (0.40, 1.32), (0.45, 1.37), (0.475, 1.45), (0.50, 1.55),
            (0.525, 1.57), (0.55, 1.53), (0.575, 1.47), (0.60, 1.39),
            (0.65, 1.25), (0.70, 1.12), (0.75, 0.98), (0.80, 0.86),
            (0.84, 0.79),
        ),
    ),
    Curve(
        "C5", "round-bilge hull, L/B 11", "dotted",
        "lowest initial RBH branch; lowest broad-hump shoulder; only curve with high-Fn upturn",
        0.20, 1.00,
        (
            (0.20, 1.00), (0.25, 1.02), (0.30, 1.08), (0.35, 1.12),
            (0.40, 1.22), (0.45, 1.40), (0.475, 1.49), (0.50, 1.51),
            (0.525, 1.47), (0.55, 1.40), (0.60, 1.27), (0.65, 1.12),
            (0.70, 0.98), (0.75, 0.91), (0.80, 0.89), (0.85, 0.89),
            (0.90, 0.91), (0.95, 0.96), (1.00, 1.09),
        ),
    ),
)

X_LEFT = 535
X_RIGHT = 2403
Y_TOP = 314
Y_BOTTOM = 1582


def pixel_x(fn: float) -> float:
    return X_LEFT + (fn - 0.1) / 0.9 * (X_RIGHT - X_LEFT)


def pixel_y(tau: float) -> float:
    return Y_BOTTOM - tau / 2.5 * (Y_BOTTOM - Y_TOP)


def guide_tau(curve: Curve, fn: float) -> float:
    return float(np.interp(
        fn,
        [point[0] for point in curve.guides],
        [point[1] for point in curve.guides],
    ))


def locate_line(gray: np.ndarray, x: float, expected_y: float) -> float:
    xi = round(x)
    centre = round(expected_y)
    best_y = centre
    best_score = -float("inf")
    for y in range(max(0, centre - 18), min(gray.shape[0], centre + 19)):
        patch = gray[max(0, y - 3):y + 4, max(0, xi - 10):xi + 11]
        ink = np.clip(188.0 - patch.astype(float), 0.0, None)
        score = float(ink.sum()) - 8.0 * (y - expected_y) ** 2
        if score > best_score:
            best_score = score
            best_y = y
    strip = gray[max(0, best_y - 4):best_y + 5, max(0, xi - 5):xi + 6]
    weights = np.clip(198.0 - strip.astype(float), 0.0, None).sum(axis=1)
    pixels = np.arange(max(0, best_y - 4), best_y + 5, dtype=float)
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
            y = locate_line(gray, x, pixel_y(guide_tau(curve, fn)))
            rows.append({
                "pass": "A",
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
            draw.ellipse((x - 2, y - 2, x + 2, y + 2), outline=colours[curve.model], width=1)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    overlay.save(OVERLAY)
    print(f"wrote {len(rows)} source-only points to {OUTPUT}")


if __name__ == "__main__":
    main()

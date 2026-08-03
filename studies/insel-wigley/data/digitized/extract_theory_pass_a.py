#!/usr/bin/env python3
"""Source-only first trace of Insel's solid C2 theory curves.

The guide ordinates were read from the scanned figures without consulting a
prediction or another digitization pass. Local pixel scoring then centres each
common-Fn sample on the printed solid line. The ignored overlay is for checking
curve identity against the source, not for comparison with the solver.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RENDER_DIR = ROOT / "tmp/pdfs/insel-followup/theory-highres"
OVERLAY_DIR = ROOT / "tmp/pdfs/insel-followup/theory-pass-a-overlays"
OUTPUT = HERE / "passes/pass_a_theory.tsv"


@dataclass(frozen=True)
class Panel:
    figure: int
    pdf_page: int
    printed_page: int
    separation: float
    x_left: int
    x_right: int
    y_top: int
    y_bottom: int
    fn_last: float
    guides: tuple[tuple[float, float], ...]


PANELS = (
    Panel(359, 368, 358, 0.2, 535, 2403, 314, 1582, 0.80, (
        (0.15, 1.20), (0.175, 1.10), (0.20, 1.35), (0.225, 0.82),
        (0.25, 1.25), (0.275, 0.80), (0.30, 1.50), (0.325, 1.05),
        (0.35, 0.75), (0.375, 1.25), (0.40, 1.78), (0.425, 1.96),
        (0.45, 1.98), (0.475, 1.90), (0.50, 1.76), (0.525, 1.62),
        (0.55, 1.50), (0.60, 1.30), (0.65, 1.14), (0.70, 1.00),
        (0.75, 0.90), (0.80, 0.80),
    )),
    Panel(360, 368, 358, 0.3, 535, 2403, 2053, 3318, 0.95, (
        (0.15, 1.40), (0.175, 1.05), (0.20, 1.27), (0.225, 1.05),
        (0.25, 0.70), (0.275, 1.52), (0.30, 1.05), (0.325, 0.55),
        (0.35, 0.90), (0.375, 1.55), (0.40, 1.72), (0.425, 1.70),
        (0.45, 1.60), (0.475, 1.47), (0.50, 1.35), (0.55, 1.16),
        (0.60, 1.04), (0.65, 0.98), (0.70, 0.95), (0.75, 0.92),
        (0.80, 0.90), (0.85, 0.87), (0.90, 0.87), (0.95, 0.90),
    )),
    Panel(361, 369, 359, 0.4, 535, 2403, 314, 1582, 0.95, (
        (0.15, 1.00), (0.175, 0.99), (0.20, 0.85), (0.225, 1.30),
        (0.25, 0.85), (0.275, 1.30), (0.30, 1.05), (0.325, 0.80),
        (0.35, 0.85), (0.375, 1.40), (0.40, 1.54), (0.425, 1.47),
        (0.45, 1.36), (0.475, 1.27), (0.50, 1.20), (0.55, 1.10),
        (0.60, 1.03), (0.65, 1.00), (0.70, 0.99), (0.75, 0.97),
        (0.80, 0.93), (0.85, 0.90), (0.90, 0.88), (0.95, 0.88),
    )),
    Panel(362, 369, 359, 0.5, 535, 2403, 2053, 3318, 0.95, (
        (0.15, 1.00), (0.175, 1.00), (0.20, 0.98), (0.225, 1.02),
        (0.25, 1.10), (0.275, 0.95), (0.30, 1.15), (0.325, 0.90),
        (0.35, 0.95), (0.375, 1.28), (0.40, 1.36), (0.425, 1.30),
        (0.45, 1.22), (0.475, 1.15), (0.50, 1.10), (0.55, 1.04),
        (0.60, 1.02), (0.65, 1.00), (0.70, 1.00), (0.75, 0.98),
        (0.80, 0.95), (0.85, 0.92), (0.90, 0.91), (0.95, 0.93),
    )),
)


def pixel_x(panel: Panel, fn: float) -> float:
    return panel.x_left + (fn - 0.1) / 0.9 * (panel.x_right - panel.x_left)


def pixel_y(panel: Panel, tau: float) -> float:
    return panel.y_bottom - tau / 2.5 * (panel.y_bottom - panel.y_top)


def guide_tau(panel: Panel, fn: float) -> float:
    xs = np.array([point[0] for point in panel.guides])
    ys = np.array([point[1] for point in panel.guides])
    return float(np.interp(fn, xs, ys))


def locate_line(gray: np.ndarray, x: float, expected_y: float) -> tuple[float, float]:
    xi = round(x)
    centre = round(expected_y)
    candidates = range(max(0, centre - 55), min(gray.shape[0], centre + 56))
    best_y = centre
    best_score = -float("inf")
    for y in candidates:
        patch = gray[max(0, y - 3):y + 4, max(0, xi - 6):xi + 7]
        ink = np.clip(185.0 - patch.astype(float), 0.0, None)
        score = float(ink.sum()) - 11.0 * abs(y - expected_y)
        if score > best_score:
            best_score = score
            best_y = y
    strip = gray[max(0, best_y - 4):best_y + 5, max(0, xi - 3):xi + 4]
    weights = np.clip(200.0 - strip.astype(float), 0.0, None).sum(axis=1)
    offsets = np.arange(max(0, best_y - 4), best_y + 5, dtype=float)
    refined_y = float((weights * offsets).sum() / weights.sum()) if weights.sum() else float(best_y)
    return float(x), refined_y


def main() -> None:
    OVERLAY_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for panel in PANELS:
        source = RENDER_DIR / f"page-{panel.pdf_page}.jpg"
        image = Image.open(source).convert("RGB")
        gray = np.asarray(image.convert("L"))
        overlay = image.copy()
        draw = ImageDraw.Draw(overlay)
        for index in range(round((panel.fn_last - 0.15) / 0.005) + 1):
            fn = round(0.15 + 0.005 * index, 3)
            x, y = locate_line(gray, pixel_x(panel, fn), pixel_y(panel, guide_tau(panel, fn)))
            rows.append({
                "pass": "A",
                "figure": panel.figure,
                "pdf_page": panel.pdf_page,
                "printed_page": panel.printed_page,
                "separation_over_length": f"{panel.separation:.1f}",
                "fn_anchor": f"{fn:.3f}",
                "pixel_x": f"{x:.3f}",
                "pixel_y": f"{y:.3f}",
                "pixel_x_left": panel.x_left,
                "pixel_x_right": panel.x_right,
                "pixel_y_top": panel.y_top,
                "pixel_y_bottom": panel.y_bottom,
            })
            draw.ellipse((x - 3, y - 3, x + 3, y + 3), outline=(220, 0, 0), width=2)
        overlay.save(OVERLAY_DIR / f"figure-{panel.figure}.png")
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} source-only points to {OUTPUT}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent source-only trace of Insel's solid C2 theory curves.

This pass uses a continuity-constrained whole-line pixel path through
independently read source checkpoints, rather than the per-anchor local
readings used by pass A. It reads only the rendered source pages.
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
OVERLAY_DIR = ROOT / "tmp/pdfs/insel-followup/theory-pass-b-overlays"
OUTPUT = HERE / "passes/pass_b_theory.tsv"


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
    checkpoints: tuple[tuple[float, float], ...]


PANELS = (
    Panel(359, 368, 358, 0.2, 536, 2402, 315, 1581, 0.80, (
        (0.15, 1.20), (0.175, 1.10), (0.20, 1.35), (0.225, 0.82),
        (0.25, 1.25), (0.275, 0.80), (0.30, 1.50), (0.325, 1.05),
        (0.35, 0.75), (0.375, 1.25), (0.40, 1.78), (0.425, 1.96),
        (0.45, 1.98), (0.475, 1.90), (0.50, 1.76), (0.55, 1.50),
        (0.60, 1.30), (0.65, 1.14), (0.70, 1.00), (0.75, 0.90),
        (0.80, 0.80),
    )),
    Panel(360, 368, 358, 0.3, 536, 2402, 2054, 3317, 0.95, (
        (0.15, 1.40), (0.175, 1.05), (0.20, 1.27), (0.225, 1.05),
        (0.25, 0.70), (0.275, 1.52), (0.30, 1.05), (0.325, 0.55),
        (0.35, 0.90), (0.375, 1.55), (0.40, 1.72), (0.425, 1.70),
        (0.45, 1.60), (0.475, 1.47), (0.50, 1.35), (0.55, 1.16),
        (0.60, 1.04), (0.65, 0.98), (0.70, 0.95), (0.75, 0.92),
        (0.80, 0.90), (0.85, 0.87), (0.90, 0.87), (0.95, 0.90),
    )),
    Panel(361, 369, 359, 0.4, 536, 2402, 315, 1581, 0.95, (
        (0.15, 1.00), (0.175, 0.99), (0.20, 0.85), (0.225, 1.30),
        (0.25, 0.85), (0.275, 1.30), (0.30, 1.05), (0.325, 0.80),
        (0.35, 0.85), (0.375, 1.40), (0.40, 1.54), (0.425, 1.47),
        (0.45, 1.36), (0.475, 1.27), (0.50, 1.20), (0.55, 1.10),
        (0.60, 1.03), (0.65, 1.00), (0.70, 0.99), (0.75, 0.97),
        (0.80, 0.93), (0.85, 0.90), (0.90, 0.88), (0.95, 0.88),
    )),
    Panel(362, 369, 359, 0.5, 536, 2402, 2054, 3317, 0.95, (
        (0.15, 1.00), (0.175, 1.00), (0.20, 0.98), (0.225, 1.02),
        (0.25, 1.10), (0.275, 0.95), (0.30, 1.15), (0.325, 0.90),
        (0.35, 0.95), (0.375, 1.28), (0.40, 1.36), (0.425, 1.30),
        (0.45, 1.22), (0.475, 1.15), (0.50, 1.10), (0.55, 1.04),
        (0.60, 1.02), (0.65, 1.00), (0.70, 1.00), (0.75, 0.98),
        (0.80, 0.95), (0.85, 0.92), (0.90, 0.91), (0.95, 0.93),
    )),
)


def x_at(panel: Panel, fn: float) -> float:
    return panel.x_left + (fn - 0.1) / 0.9 * (panel.x_right - panel.x_left)


def y_at(panel: Panel, tau: float) -> float:
    return panel.y_bottom - tau / 2.5 * (panel.y_bottom - panel.y_top)


def guide_y(panel: Panel, xs: np.ndarray) -> np.ndarray:
    fn = 0.1 + (xs - panel.x_left) * 0.9 / (panel.x_right - panel.x_left)
    guide_fn = np.array([point[0] for point in panel.checkpoints])
    guide_tau = np.array([point[1] for point in panel.checkpoints])
    return y_at(panel, np.interp(fn, guide_fn, guide_tau))


def local_ink(gray: np.ndarray, xs: np.ndarray, y0: int, y1: int) -> np.ndarray:
    darkness = np.clip(190.0 - gray[y0:y1, :].astype(float), 0.0, None)
    scores = np.zeros((len(xs), y1 - y0), dtype=np.float32)
    for column, x in enumerate(xs):
        patch = darkness[:, max(0, x - 2):x + 3]
        vertical = np.pad(patch.sum(axis=1), (2, 2), mode="edge")
        scores[column] = sum(vertical[offset:offset + y1 - y0] for offset in range(5))
    return scores


def trace(panel: Panel, gray: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    first_x = round(x_at(panel, 0.15))
    last_x = round(x_at(panel, panel.fn_last))
    xs = np.arange(first_x, last_x + 1, dtype=int)
    y0 = round(y_at(panel, 2.15))
    y1 = round(y_at(panel, 0.45)) + 1
    scores = local_ink(gray, xs, y0, y1)
    height = y1 - y0
    relative_y = np.arange(height, dtype=float)
    source_guide = guide_y(panel, xs) - y0
    scores -= 1.0 * (relative_y[None, :] - source_guide[:, None]) ** 2
    back = np.zeros((len(xs), height), dtype=np.int16)
    previous = np.full(height, -1.0e9, dtype=np.float32)
    start = source_guide[0]
    previous[:] = scores[0] - 80.0 * (np.arange(height) - start) ** 2
    shifts = np.arange(-16, 17)
    penalties = 2.0 * shifts.astype(float) ** 2
    for column in range(1, len(xs)):
        choices = np.full((len(shifts), height), -1.0e9, dtype=np.float32)
        for row, (shift, penalty) in enumerate(zip(shifts, penalties, strict=True)):
            if shift < 0:
                choices[row, :shift] = previous[-shift:] - penalty
            elif shift > 0:
                choices[row, shift:] = previous[:-shift] - penalty
            else:
                choices[row] = previous
        selected = choices.argmax(axis=0)
        back[column] = selected.astype(np.int16)
        previous = scores[column] + choices[selected, np.arange(height)]
    path = np.zeros(len(xs), dtype=int)
    path[-1] = int(previous.argmax())
    for column in range(len(xs) - 1, 0, -1):
        path[column - 1] = path[column] - int(shifts[back[column, path[column]]])
    return xs, path + y0


def refine_sample(gray: np.ndarray, x: float, expected_y: float) -> float:
    xi = round(x)
    centre = round(expected_y)
    best_y = centre
    best_score = -float("inf")
    for y in range(max(0, centre - 55), min(gray.shape[0], centre + 56)):
        patch = gray[max(0, y - 4):y + 5, max(0, xi - 5):xi + 6]
        ink = np.clip(180.0 - patch.astype(float), 0.0, None)
        score = float(ink.sum()) - 8.0 * abs(y - expected_y)
        if score > best_score:
            best_score = score
            best_y = y
    return float(best_y)


def main() -> None:
    OVERLAY_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for panel in PANELS:
        source = RENDER_DIR / f"page-{panel.pdf_page}.jpg"
        image = Image.open(source).convert("RGB")
        gray = np.asarray(image.convert("L"))
        xs, ys = trace(panel, gray)
        overlay = image.copy()
        draw = ImageDraw.Draw(overlay)
        for index in range(round((panel.fn_last - 0.15) / 0.005) + 1):
            fn = round(0.15 + 0.005 * index, 3)
            x = x_at(panel, fn)
            path_y = float(np.interp(x, xs, ys))
            checkpoint_y = float(guide_y(panel, np.array([x]))[0])
            expected_y = path_y if abs(path_y - checkpoint_y) <= 20 else checkpoint_y
            y = refine_sample(gray, x, expected_y)
            rows.append({
                "pass": "B",
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
            draw.ellipse((x - 2, y - 2, x + 2, y + 2), outline=(0, 80, 230), width=1)
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

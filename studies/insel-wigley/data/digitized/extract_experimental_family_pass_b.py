#!/usr/bin/env python3
"""Independent source-only principal-hump trace of Figures 347--350."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE_DIR = ROOT / "tmp/pdfs/insel-final/experimental-figures"
OVERLAY_DIR = ROOT / "tmp/pdfs/insel-final/experimental-pass-b"
OUTPUT = HERE / "passes/final_experimental_family_pass_b.tsv"


@dataclass(frozen=True)
class Trace:
    model: str
    hull: str
    line_style: str
    relative_position: str
    checkpoints: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class Panel:
    figure: int
    pdf_page: int
    printed_page: int
    separation: float
    y_top: int
    y_bottom: int
    traces: tuple[Trace, ...]


def trace(model: str, hull: str, style: str, position: str, values: tuple[float, ...]) -> Trace:
    fn = (0.35, 0.375, 0.40, 0.425, 0.45, 0.475, 0.50, 0.525, 0.55)
    return Trace(model, hull, style, position, tuple(zip(fn, values, strict=True)))


PANELS = (
    Panel(347, 362, 352, 0.2, 409, 2104, (
        trace("C2", "Wigley hull", "solid", "solid trough then second-highest close-spacing crest", (1.50, 0.91, 1.13, 1.44, 1.65, 1.75, 1.77, 1.74, 1.68)),
        trace("C3", "round-bilge hull, L/B 7", "long dashed", "uppermost right-shifted principal crest", (1.74, 1.67, 1.58, 1.55, 1.62, 1.75, 1.85, 1.88, 1.86)),
        trace("C4", "round-bilge hull, L/B 9", "dash-dot", "low central shoulder and shallow rise", (1.39, 1.16, 1.20, 1.41, 1.49, 1.51, 1.53, 1.55, 1.57)),
        trace("C5", "round-bilge hull, L/B 11", "dotted", "dotted crest below the two upper curves", (1.35, 1.12, 1.13, 1.39, 1.59, 1.68, 1.70, 1.66, 1.60)),
    )),
    Panel(348, 362, 352, 0.3, 2729, 4425, (
        trace("C2", "Wigley hull", "solid", "solid principal crest below C3", (1.09, 0.73, 1.29, 1.51, 1.54, 1.53, 1.47, 1.37, 1.30)),
        trace("C3", "round-bilge hull, L/B 7", "long dashed", "uppermost broad principal crest", (1.27, 1.17, 1.15, 1.37, 1.55, 1.67, 1.65, 1.57, 1.47)),
        trace("C4", "round-bilge hull, L/B 9", "dash-dot", "left shoulder crossing into the main family", (1.27, 0.71, 1.41, 1.59, 1.55, 1.47, 1.39, 1.34, 1.30)),
        trace("C5", "round-bilge hull, L/B 11", "dotted", "lowest smooth principal crest", (1.14, 0.69, 1.24, 1.41, 1.47, 1.45, 1.39, 1.33, 1.28)),
    )),
    Panel(349, 363, 353, 0.4, 399, 2092, (
        trace("C2", "Wigley hull", "solid", "highest narrow crest", (1.09, 0.76, 1.34, 1.59, 1.39, 1.27, 1.21, 1.17, 1.13)),
        trace("C3", "round-bilge hull, L/B 7", "long dashed", "right-shifted NPL maximum", (1.31, 1.06, 1.12, 1.29, 1.44, 1.49, 1.44, 1.34, 1.26)),
        trace("C4", "round-bilge hull, L/B 9", "dash-dot", "middle broad shoulder", (1.11, 0.77, 1.11, 1.27, 1.34, 1.33, 1.29, 1.23, 1.18)),
        trace("C5", "round-bilge hull, L/B 11", "dotted", "lowest broad shoulder", (1.07, 0.77, 1.07, 1.24, 1.32, 1.31, 1.26, 1.21, 1.16)),
    )),
    Panel(350, 363, 353, 0.5, 2727, 4405, (
        trace("C2", "Wigley hull", "solid", "isolated upper crest", (1.04, 0.91, 1.34, 1.53, 1.27, 1.17, 1.13, 1.09, 1.06)),
        trace("C3", "round-bilge hull, L/B 7", "long dashed", "right-shifted lower crest", (0.99, 0.83, 0.92, 1.07, 1.19, 1.22, 1.19, 1.14, 1.09)),
        trace("C4", "round-bilge hull, L/B 9", "dash-dot", "middle NPL shoulder", (0.99, 0.91, 1.07, 1.19, 1.24, 1.24, 1.21, 1.17, 1.13)),
        trace("C5", "round-bilge hull, L/B 11", "dotted", "smooth lower crest", (1.01, 0.93, 1.09, 1.19, 1.24, 1.25, 1.23, 1.19, 1.14)),
    )),
)

X_LEFT = 731
X_RIGHT = 3228


def pixel_x(fn: float) -> float:
    return X_LEFT + (fn - 0.1) / 0.9 * (X_RIGHT - X_LEFT)


def pixel_y(panel: Panel, tau: float) -> float:
    return panel.y_bottom - tau / 2.5 * (panel.y_bottom - panel.y_top)


def expected_tau(trace_: Trace, fn: float) -> float:
    return float(np.interp(fn, [p[0] for p in trace_.checkpoints], [p[1] for p in trace_.checkpoints]))


def locate_line(gray: np.ndarray, x: float, expected_y: float) -> float:
    xi = round(x)
    centre = round(expected_y)
    best_y = centre
    best_score = -float("inf")
    for y in range(max(0, centre - 22), min(gray.shape[0], centre + 23)):
        patch = gray[max(0, y - 5):y + 6, max(0, xi - 19):xi + 20]
        darkness = np.clip(184.0 - patch.astype(float), 0.0, None)
        response = darkness.sum(axis=1)
        score = float(response.sum() + 0.4 * response.max()) - 6.0 * (y - expected_y) ** 2
        if score > best_score:
            best_score = score
            best_y = y
    strip = gray[max(0, best_y - 5):best_y + 6, max(0, xi - 8):xi + 9]
    weights = np.clip(202.0 - strip.astype(float), 0.0, None).sum(axis=1)
    pixels = np.arange(max(0, best_y - 5), best_y + 6, dtype=float)
    return float((weights * pixels).sum() / weights.sum()) if weights.sum() else float(best_y)


def main() -> None:
    OVERLAY_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str | int]] = []
    colours = {"C2": (220, 0, 0), "C3": (0, 90, 230), "C4": (0, 160, 70), "C5": (190, 0, 190)}
    for panel in PANELS:
        image = Image.open(SOURCE_DIR / f"page-{panel.pdf_page}.jpg").convert("RGB")
        gray = np.asarray(image.convert("L"))
        overlay = image.copy()
        draw = ImageDraw.Draw(overlay)
        for trace_ in panel.traces:
            for index in range(41):
                fn = round(0.35 + 0.005 * index, 3)
                x = pixel_x(fn)
                y = locate_line(gray, x, pixel_y(panel, expected_tau(trace_, fn)))
                rows.append({
                    "pass": "B", "figure": panel.figure, "pdf_page": panel.pdf_page,
                    "printed_page": panel.printed_page,
                    "separation_over_length": f"{panel.separation:.1f}",
                    "model": trace_.model, "hull": trace_.hull,
                    "line_style": trace_.line_style,
                    "relative_position_description": trace_.relative_position,
                    "fn_anchor": f"{fn:.3f}", "pixel_x": f"{x:.3f}",
                    "pixel_y": f"{y:.3f}", "pixel_x_left": X_LEFT,
                    "pixel_x_right": X_RIGHT, "pixel_y_top": panel.y_top,
                    "pixel_y_bottom": panel.y_bottom,
                })
                draw.rectangle((x - 2, y - 2, x + 2, y + 2), outline=colours[trace_.model], width=1)
        overlay.save(OVERLAY_DIR / f"figure-{panel.figure}.png")
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} source-only points to {OUTPUT}")


if __name__ == "__main__":
    main()

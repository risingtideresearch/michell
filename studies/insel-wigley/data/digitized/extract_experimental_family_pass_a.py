#!/usr/bin/env python3
"""First source-only principal-hump trace of Figures 347--350."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE_DIR = ROOT / "tmp/pdfs/insel-final/experimental-figures"
OVERLAY_DIR = ROOT / "tmp/pdfs/insel-final/experimental-pass-a"
OUTPUT = HERE / "passes/final_experimental_family_pass_a.tsv"


@dataclass(frozen=True)
class Curve:
    model: str
    hull: str
    line_style: str
    relative_position: str
    guides: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class Panel:
    figure: int
    pdf_page: int
    printed_page: int
    separation: float
    y_top: int
    y_bottom: int
    curves: tuple[Curve, ...]


def curve(model: str, hull: str, style: str, position: str, values: tuple[float, ...]) -> Curve:
    fn = (0.35, 0.375, 0.40, 0.425, 0.45, 0.475, 0.50, 0.525, 0.55)
    return Curve(model, hull, style, position, tuple(zip(fn, values, strict=True)))


PANELS = (
    Panel(347, 362, 352, 0.2, 407, 2106, (
        curve("C2", "Wigley hull", "solid", "deep trough near Fn 0.38; solid broad crest below dashed C3", (1.52, 0.90, 1.12, 1.45, 1.66, 1.76, 1.78, 1.75, 1.69)),
        curve("C3", "round-bilge hull, L/B 7", "long dashed", "highest broad crest, displaced right of C2", (1.75, 1.68, 1.58, 1.55, 1.63, 1.76, 1.86, 1.89, 1.87)),
        curve("C4", "round-bilge hull, L/B 9", "dash-dot", "lower central shoulder, then middle descending branch", (1.40, 1.15, 1.20, 1.42, 1.50, 1.52, 1.53, 1.55, 1.58)),
        curve("C5", "round-bilge hull, L/B 11", "dotted", "dotted crest below C2 and C3", (1.36, 1.12, 1.12, 1.40, 1.60, 1.70, 1.72, 1.68, 1.62)),
    )),
    Panel(348, 362, 352, 0.3, 2727, 4427, (
        curve("C2", "Wigley hull", "solid", "solid crest below dashed C3", (1.10, 0.72, 1.30, 1.52, 1.55, 1.54, 1.48, 1.38, 1.31)),
        curve("C3", "round-bilge hull, L/B 7", "long dashed", "highest broad crest and rightmost peak", (1.28, 1.18, 1.15, 1.38, 1.56, 1.68, 1.66, 1.58, 1.48)),
        curve("C4", "round-bilge hull, L/B 9", "dash-dot", "sharp central rise, upper-left shoulder", (1.28, 0.70, 1.42, 1.60, 1.56, 1.48, 1.40, 1.34, 1.31)),
        curve("C5", "round-bilge hull, L/B 11", "dotted", "lowest broad crest through the main hump", (1.15, 0.68, 1.25, 1.42, 1.48, 1.46, 1.40, 1.34, 1.29)),
    )),
    Panel(349, 363, 353, 0.4, 397, 2094, (
        curve("C2", "Wigley hull", "solid", "highest narrow crest near Fn 0.41", (1.10, 0.75, 1.35, 1.60, 1.40, 1.28, 1.22, 1.17, 1.13)),
        curve("C3", "round-bilge hull, L/B 7", "long dashed", "right-shifted broad NPL crest", (1.32, 1.05, 1.12, 1.30, 1.45, 1.50, 1.45, 1.35, 1.26)),
        curve("C4", "round-bilge hull, L/B 9", "dash-dot", "lower broad NPL shoulder", (1.12, 0.76, 1.12, 1.28, 1.35, 1.34, 1.30, 1.24, 1.19)),
        curve("C5", "round-bilge hull, L/B 11", "dotted", "lowest smooth NPL crest", (1.08, 0.76, 1.08, 1.25, 1.33, 1.32, 1.27, 1.22, 1.17)),
    )),
    Panel(350, 363, 353, 0.5, 2725, 4407, (
        curve("C2", "Wigley hull", "solid", "isolated highest crest near Fn 0.41", (1.05, 0.90, 1.35, 1.54, 1.28, 1.18, 1.14, 1.10, 1.07)),
        curve("C3", "round-bilge hull, L/B 7", "long dashed", "right-shifted low NPL crest", (1.00, 0.82, 0.92, 1.08, 1.20, 1.23, 1.20, 1.15, 1.10)),
        curve("C4", "round-bilge hull, L/B 9", "dash-dot", "middle NPL shoulder", (1.00, 0.90, 1.08, 1.20, 1.25, 1.25, 1.22, 1.18, 1.14)),
        curve("C5", "round-bilge hull, L/B 11", "dotted", "smooth NPL crest below C2", (1.02, 0.92, 1.10, 1.20, 1.25, 1.26, 1.24, 1.20, 1.15)),
    )),
)

X_LEFT = 729
X_RIGHT = 3230


def pixel_x(fn: float) -> float:
    return X_LEFT + (fn - 0.1) / 0.9 * (X_RIGHT - X_LEFT)


def pixel_y(panel: Panel, tau: float) -> float:
    return panel.y_bottom - tau / 2.5 * (panel.y_bottom - panel.y_top)


def guide_tau(curve_: Curve, fn: float) -> float:
    return float(np.interp(fn, [p[0] for p in curve_.guides], [p[1] for p in curve_.guides]))


def locate_line(gray: np.ndarray, x: float, expected_y: float) -> float:
    xi = round(x)
    centre = round(expected_y)
    best_y = centre
    best_score = -float("inf")
    for y in range(max(0, centre - 22), min(gray.shape[0], centre + 23)):
        patch = gray[max(0, y - 4):y + 5, max(0, xi - 14):xi + 15]
        ink = np.clip(188.0 - patch.astype(float), 0.0, None)
        score = float(ink.sum()) - 7.0 * (y - expected_y) ** 2
        if score > best_score:
            best_score = score
            best_y = y
    strip = gray[max(0, best_y - 4):best_y + 5, max(0, xi - 6):xi + 7]
    weights = np.clip(200.0 - strip.astype(float), 0.0, None).sum(axis=1)
    pixels = np.arange(max(0, best_y - 4), best_y + 5, dtype=float)
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
        for curve_ in panel.curves:
            for index in range(41):
                fn = round(0.35 + 0.005 * index, 3)
                x = pixel_x(fn)
                y = locate_line(gray, x, pixel_y(panel, guide_tau(curve_, fn)))
                rows.append({
                    "pass": "A", "figure": panel.figure, "pdf_page": panel.pdf_page,
                    "printed_page": panel.printed_page,
                    "separation_over_length": f"{panel.separation:.1f}",
                    "model": curve_.model, "hull": curve_.hull,
                    "line_style": curve_.line_style,
                    "relative_position_description": curve_.relative_position,
                    "fn_anchor": f"{fn:.3f}", "pixel_x": f"{x:.3f}",
                    "pixel_y": f"{y:.3f}", "pixel_x_left": X_LEFT,
                    "pixel_x_right": X_RIGHT, "pixel_y_top": panel.y_top,
                    "pixel_y_bottom": panel.y_bottom,
                })
                draw.ellipse((x - 2, y - 2, x + 2, y + 2), outline=colours[curve_.model], width=1)
        overlay.save(OVERLAY_DIR / f"figure-{panel.figure}.png")
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} source-only points to {OUTPUT}")


if __name__ == "__main__":
    main()

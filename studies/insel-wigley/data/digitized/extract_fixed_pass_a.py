#!/usr/bin/env python3
"""Extract blind pass-A marker centers from rendered Insel C2-FX figures.

This is image-component analysis, not OCR.  It deliberately retains calibrated
pixel coordinates and does not snap points to a presumed Froude-number grid.
Every accepted component was checked against the rendered source plot.
"""

from __future__ import annotations

import csv
from collections import deque
from pathlib import Path

import numpy as np
from PIL import Image


STUDY = Path(__file__).resolve().parents[2]
RENDERS = STUDY / "tmp/pdfs/insel-wigley/thesis-c2-highres"
OUTPUT = Path(__file__).resolve().parent / "passes"

# (configuration, CT figure, CWP figure, PDF page, printed page)
SERIES = {
    254: ("monohull", 135, 136, 254, 244),
    255: ("s_l_0_2", 137, 138, 255, 245),
    256: ("s_l_0_3", 139, 140, 256, 246),
    257: ("s_l_0_4", 141, 142, 257, 247),
    258: ("s_l_0_5", 143, 144, 258, 248),
}

# Axis endpoints read independently from the rendered frame lines.
# Each tuple is (x at Fn=.1, x at Fn=1, y at max coefficient,
# y at min coefficient).  Top plots span CT=.003..013; bottom plots span
# CWP=0..010.
CT_AXES = {
    254: (536, 2418, 307, 1576),
    255: (535, 2412, 302, 1568),
    256: (538, 2421, 336, 1604),
    257: (548, 2427, 309, 1577),
    258: (535, 2416, 295, 1563),
}
CWP_AXES = {
    254: (532, 2409, 2047, 3317),
    255: (529, 2403, 2038, 3307),
    256: (534, 2411, 2073, 3341),
    257: (538, 2415, 2046, 3317),
    258: (528, 2406, 2035, 3306),
}

FN_UNCERTAINTY = 0.003
COEFFICIENT_UNCERTAINTY = 0.00005


def coefficient(x: float, y: float, axes: tuple[int, int, int, int],
                minimum: float, maximum: float) -> tuple[float, float]:
    x0, x1, y_top, y_bottom = axes
    fn = 0.1 + (x - x0) * 0.9 / (x1 - x0)
    value = minimum + (y_bottom - y) * (maximum - minimum) / (y_bottom - y_top)
    return fn, value


def white_holes(image: np.ndarray, page: int) -> list[tuple[float, float]]:
    """Locate the enclosed white centers of measured hollow CT squares."""
    x0, x1, y_top, y_bottom = CT_AXES[page]
    # The measured CT curve stays in this band.  Excluding the legend also
    # prevents letter counters from becoming candidates.
    y0 = int(y_bottom - (0.0105 - 0.003) / 0.010 * (y_bottom - y_top))
    y1 = int(y_bottom - (0.0050 - 0.003) / 0.010 * (y_bottom - y_top))
    x_lo, x_hi = x0 + 5, x1 - 5
    y_lo, y_hi = max(y_top + 5, y0), min(y_bottom - 5, y1)
    white = image[y_lo:y_hi, x_lo:x_hi] > 180
    height, width = white.shape
    visited = np.zeros_like(white, dtype=bool)
    holes: list[tuple[float, float, int]] = []

    for row in range(height):
        for col in range(width):
            if not white[row, col] or visited[row, col]:
                continue
            queue = deque([(row, col)])
            visited[row, col] = True
            pixels: list[tuple[int, int]] = []
            touches_edge = False
            while queue:
                current_row, current_col = queue.popleft()
                pixels.append((current_row, current_col))
                touches_edge |= (
                    current_row == 0
                    or current_col == 0
                    or current_row == height - 1
                    or current_col == width - 1
                )
                for delta_row, delta_col in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    next_row = current_row + delta_row
                    next_col = current_col + delta_col
                    if (
                        0 <= next_row < height
                        and 0 <= next_col < width
                        and white[next_row, next_col]
                        and not visited[next_row, next_col]
                    ):
                        visited[next_row, next_col] = True
                        queue.append((next_row, next_col))
            if touches_edge:
                continue
            rows = [pixel[0] for pixel in pixels]
            cols = [pixel[1] for pixel in pixels]
            box_width = max(cols) - min(cols) + 1
            box_height = max(rows) - min(rows) + 1
            area = len(pixels)
            if 3 <= box_width <= 18 and 3 <= box_height <= 20 and 8 <= area <= 220:
                holes.append(
                    (
                        x_lo + sum(cols) / area,
                        y_lo + sum(rows) / area,
                        area,
                    )
                )

    # A fitted curve can split one square's white center into two regions.
    # Merge only regions closer than 12 rendered pixels; adjacent experimental
    # points are farther apart.
    used = [False] * len(holes)
    markers: list[tuple[float, float]] = []
    for index, hole in enumerate(holes):
        if used[index]:
            continue
        group = [hole]
        used[index] = True
        changed = True
        while changed:
            changed = False
            center_x = sum(item[0] for item in group) / len(group)
            center_y = sum(item[1] for item in group) / len(group)
            for candidate_index, candidate in enumerate(holes):
                if used[candidate_index]:
                    continue
                if (candidate[0] - center_x) ** 2 + (candidate[1] - center_y) ** 2 <= 12**2:
                    used[candidate_index] = True
                    group.append(candidate)
                    changed = True
        total_area = sum(item[2] for item in group)
        markers.append(
            (
                sum(item[0] * item[2] for item in group) / total_area,
                sum(item[1] * item[2] for item in group) / total_area,
            )
        )
    return sorted(markers)


def rectangle_sum(integral: np.ndarray, x: int, y: int, radius_x: int,
                  radius_y: int) -> int:
    return int(
        integral[y + radius_y + 1, x + radius_x + 1]
        - integral[y - radius_y, x + radius_x + 1]
        - integral[y + radius_y + 1, x - radius_x]
        + integral[y - radius_y, x - radius_x]
    )


def accept_cwp(page: int, inner: int, outer: int) -> bool:
    # Page-specific cutoffs reflect scan contrast.  The borderline components
    # were accepted or rejected only after visual overlay on the source curve.
    if page == 254:
        return inner >= 95
    if page == 255:
        return inner >= 102
    if page == 256:
        return inner > 97 or outer >= 220
    if page == 257:
        return inner >= 97
    if page == 258:
        return inner >= 96
    raise AssertionError(page)


def solid_markers(image: np.ndarray, page: int) -> list[tuple[float, float]]:
    """Locate dense rectangular measured CWP symbols, excluding thin curves."""
    x0, x1, y_top, y_bottom = CWP_AXES[page]
    black = (image < 120).astype(np.int32)
    integral = np.pad(black, ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    y_min = int(y_bottom - 0.0055 / 0.010 * (y_bottom - y_top))
    candidates: list[tuple[int, int, int, int]] = []
    for y in range(y_min, y_bottom - 3):
        for x in range(x0 + 10, int(x0 + 0.82 * (x1 - x0))):
            inner = rectangle_sum(integral, x, y, 3, 7)
            if inner >= 95:
                outer = rectangle_sum(integral, x, y, 7, 13)
                candidates.append((inner, outer, x, y))
    candidates.sort(reverse=True)

    separated: list[tuple[int, int, int, int]] = []
    for candidate in candidates:
        _, _, x, y = candidate
        if all(
            (x - accepted_x) ** 2 / 18**2 + (y - accepted_y) ** 2 / 22**2 > 1
            for _, _, accepted_x, accepted_y in separated
        ):
            separated.append(candidate)

    markers = [
        (float(x), float(y))
        for inner, outer, x, y in separated
        if accept_cwp(page, inner, outer)
    ]
    return sorted(markers)


def write_pass(page: int, ct_markers: list[tuple[float, float]],
               cwp_markers: list[tuple[float, float]]) -> None:
    configuration, ct_figure, cwp_figure, pdf_page, printed_page = SERIES[page]
    destination = OUTPUT / f"pass_a_fixed_{configuration}.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        handle.write("# source: Mustafa Insel, 1990 PhD thesis\n")
        handle.write(
            f"# source_location: printed page {printed_page}; PDF page {pdf_page}; "
            f"Figures {ct_figure} and {cwp_figure}\n"
        )
        handle.write("# pass: A; image-component-assisted, visually checked; no OCR\n")
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "configuration",
                "attitude",
                "observable",
                "replicate_id",
                "x_px",
                "y_px",
                "fn",
                "coefficient",
                "fn_digitization_uncertainty",
                "coefficient_digitization_uncertainty",
                "source_figure",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        for observable, figure, axes, minimum, maximum, markers in (
            ("ct", ct_figure, CT_AXES[page], 0.003, 0.013, ct_markers),
            ("cwp", cwp_figure, CWP_AXES[page], 0.0, 0.010, cwp_markers),
        ):
            for replicate_id, (x, y) in enumerate(markers, start=1):
                fn, value = coefficient(x, y, axes, minimum, maximum)
                writer.writerow(
                    {
                        "configuration": configuration,
                        "attitude": "fixed",
                        "observable": observable,
                        "replicate_id": replicate_id,
                        "x_px": f"{x:.2f}",
                        "y_px": f"{y:.2f}",
                        "fn": f"{fn:.6f}",
                        "coefficient": f"{value:.7f}",
                        "fn_digitization_uncertainty": f"{FN_UNCERTAINTY:.3f}",
                        "coefficient_digitization_uncertainty": f"{COEFFICIENT_UNCERTAINTY:.5f}",
                        "source_figure": figure,
                    }
                )


def main() -> None:
    for page in SERIES:
        image = np.asarray(Image.open(RENDERS / f"page-{page}.jpg").convert("L"))
        write_pass(page, white_holes(image, page), solid_markers(image, page))


if __name__ == "__main__":
    main()

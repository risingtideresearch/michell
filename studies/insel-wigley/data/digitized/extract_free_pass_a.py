#!/usr/bin/env python3
"""Extract blind pass-A marker centers from rendered free-attitude C2 figures.

This script uses image components only (no OCR), retains unsnapped calibrated
coordinates, and writes the independently selected marker centers used by pass
A.  Page-specific acceptance rules are recorded beside the relevant detector;
each accepted center was checked on a labelled source overlay.
"""

from __future__ import annotations

import csv
from collections import deque
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


STUDY = Path(__file__).resolve().parents[2]
RENDERS = STUDY / "tmp/pdfs/insel-wigley/thesis-c2-highres"
OUTPUT = Path(__file__).resolve().parent / "passes"
OVERLAYS = STUDY / "tmp/digitization-pass-a"

# page: (configuration, resistance figure, attitude figure, CWP figure,
#        printed resistance/attitude page, printed CWP page)
SERIES = {
    268: ("monohull", 161, 162, 163, 258, 259),
    270: ("s_l_0_2", 165, 166, 167, 260, 261),
    272: ("s_l_0_3", 169, 170, 171, 262, 263),
    274: ("s_l_0_4", 173, 174, 175, 264, 265),
    276: ("s_l_0_5", 177, 178, 179, 266, 267),
}

# (x at Fn=.1, x at Fn=1, y at ordinate maximum, y at ordinate minimum)
CT_AXES = {
    268: (530, 2406, 278, 1543),
    270: (539, 2412, 337, 1603),
    272: (542, 2427, 292, 1563),
    274: (514, 2404, 202, 1468),
    276: (518, 2404, 282, 1556),
}
CWP_AXES = {
    269: (521, 2403, 270, 1535),
    271: (523, 2405, 244, 1514),
    273: (521, 2409, 265, 1535),
    275: (523, 2410, 243, 1511),
    277: (525, 2411, 244, 1512),
}
ATTITUDE_AXES = {
    268: (524, 2398, 2080, 3350),
    270: (523, 2393, 2145, 3419),
    272: (536, 2415, 2037, 3317),
    274: (506, 2388, 1942, 3288),
    276: (509, 2391, 2025, 3296),
}

FN_UNCERTAINTY = 0.003
COEFFICIENT_UNCERTAINTY = 0.00005
TRIM_UNCERTAINTY_DEG = 0.03
SINKAGE_UNCERTAINTY = 0.001

# Independent visual guides through the measured marker clouds.  They are
# used only to reject hollow-square/cross-shaped features belonging to the
# other observable or to plot lettering; final values retain the detected
# marker centers rather than values from these guides.
TRIM_GUIDES = {
    268: [(524, 2930), (1050, 2930), (1220, 2740), (1400, 2600), (1750, 2555), (2398, 2440)],
    270: [(523, 2975), (1047, 2980), (1107, 2910), (1162, 2790), (1280, 2560), (1390, 2490), (1625, 2490)],
    272: [(536, 2895), (1064, 2885), (1124, 2780), (1233, 2596), (1290, 2521), (1480, 2520), (1800, 2510), (2415, 2405)],
    274: [(506, 2770), (980, 2780), (1037, 2760), (1095, 2670), (1151, 2576), (1210, 2536), (1320, 2460), (1460, 2420), (1700, 2413), (2388, 2325)],
    276: [(509, 2863), (1043, 2852), (1100, 2770), (1214, 2640), (1330, 2560), (1463, 2555), (1615, 2515), (1735, 2480), (2391, 2405)],
}
SINKAGE_GUIDES = {
    268: [(524, 2920), (800, 2850), (1050, 2710), (1200, 2580), (1350, 2550), (1500, 2600), (1800, 2690), (2398, 2770)],
    270: [(523, 2960), (800, 2880), (1000, 2730), (1150, 2470), (1250, 2390), (1350, 2450), (1500, 2680), (1900, 2830)],
    272: [(536, 2880), (800, 2790), (1000, 2675), (1150, 2450), (1250, 2400), (1350, 2500), (1600, 2700), (1800, 2740), (2100, 2890), (2415, 2990)],
    274: [(506, 2760), (850, 2650), (1000, 2540), (1150, 2350), (1250, 2335), (1380, 2420), (1600, 2530), (1900, 2630), (2388, 2840)],
    276: [(509, 2835), (800, 2760), (1050, 2600), (1200, 2435), (1350, 2470), (1500, 2550), (1750, 2610), (2000, 2700), (2391, 2860)],
}
CWP_VISUAL_MAX = {269: 0.0049, 271: 0.0036, 273: 0.0062, 275: 0.0060, 277: 0.0054}


def calibrate(x: float, y: float, axes: tuple[int, int, int, int],
              minimum: float, maximum: float) -> tuple[float, float]:
    x0, x1, y_top, y_bottom = axes
    fn = 0.1 + (x - x0) * 0.9 / (x1 - x0)
    value = minimum + (y_bottom - y) * (maximum - minimum) / (y_bottom - y_top)
    return fn, value


def enclosed_white_centers(image: np.ndarray, axes: tuple[int, int, int, int],
                           *, exclude_legend: bool) -> list[tuple[float, float]]:
    """Find hollow square centers inside a calibrated plotting frame."""
    x0, x1, y_top, y_bottom = axes
    x_lo, x_hi = x0 + 5, x1 - 5
    y_lo, y_hi = y_top + 5, y_bottom - 5
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
            x = x_lo + sum(cols) / area
            y = y_lo + sum(rows) / area
            if exclude_legend and x < x0 + 760 and y < y_top + 330:
                continue
            if 3 <= box_width <= 19 and 3 <= box_height <= 22 and 8 <= area <= 250:
                holes.append((x, y, area))

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


def dense_candidates(image: np.ndarray, axes: tuple[int, int, int, int],
                     *, threshold: int = 95) -> list[tuple[int, int, int, int, int, int]]:
    """Return NMS centers with vertical/horizontal density diagnostics."""
    x0, x1, y_top, y_bottom = axes
    black = (image < 120).astype(np.int32)
    integral = np.pad(black, ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    candidates: list[tuple[int, int, int, int, int, int]] = []
    for y in range(y_top + 12, y_bottom - 13):
        for x in range(x0 + 10, x1 - 10):
            vertical = rectangle_sum(integral, x, y, 3, 7)
            if vertical >= threshold:
                outer = rectangle_sum(integral, x, y, 7, 13)
                horizontal = rectangle_sum(integral, x, y, 7, 3)
                tall = rectangle_sum(integral, x, y, 2, 13)
                candidates.append((vertical, outer, horizontal, tall, x, y))
    candidates.sort(reverse=True)
    separated: list[tuple[int, int, int, int, int, int]] = []
    for candidate in candidates:
        *_, x, y = candidate
        if all(
            (x - accepted_x) ** 2 / 18**2 + (y - accepted_y) ** 2 / 22**2 > 1
            for *_, accepted_x, accepted_y in separated
        ):
            separated.append(candidate)
    return separated


def crossed_candidates(image: np.ndarray, axes: tuple[int, int, int, int]
                       ) -> list[tuple[int, int, int, int, int, int]]:
    """Locate the two-diagonal measured sinkage symbols."""
    x0, x1, y_top, y_bottom = axes
    black = image < 120
    dilated = np.zeros_like(black)
    for delta_y in (-1, 0, 1):
        for delta_x in (-1, 0, 1):
            dilated |= np.roll(black, (delta_y, delta_x), axis=(0, 1))
    diag_down = np.zeros(image.shape, dtype=np.uint8)
    diag_up = np.zeros(image.shape, dtype=np.uint8)
    for offset in range(-7, 8):
        diag_down += np.roll(dilated, (offset, offset), axis=(0, 1))
        diag_up += np.roll(dilated, (-offset, offset), axis=(0, 1))

    integral = np.pad(black.astype(np.int32), ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    candidates: list[tuple[int, int, int, int, int, int]] = []
    rows, columns = np.nonzero(
        (diag_down >= 13)
        & (diag_up >= 13)
    )
    for y, x in zip(rows, columns, strict=True):
        if not (x0 + 12 <= x < x1 - 12 and y_top + 12 <= y < y_bottom - 12):
            continue
        center = rectangle_sum(integral, int(x), int(y), 4, 4)
        if center >= 28:
            candidates.append(
                (int(diag_down[y, x] + diag_up[y, x]), center, 0, 0,
                 int(x), int(y))
            )
    candidates.sort(reverse=True)
    separated: list[tuple[int, int, int, int, int, int]] = []
    for candidate in candidates:
        *_, x, y = candidate
        if all(
            (x - accepted_x) ** 2 / 18**2 + (y - accepted_y) ** 2 / 18**2 > 1
            for *_, accepted_x, accepted_y in separated
        ):
            separated.append(candidate)
    return separated


def make_overlay(page: int, name: str, candidates: list[tuple[int, int, int, int, int, int]]) -> None:
    image = Image.open(RENDERS / f"page-{page}.jpg").convert("RGB")
    draw = ImageDraw.Draw(image)
    for index, (*_, x, y) in enumerate(candidates, start=1):
        draw.ellipse((x - 12, y - 12, x + 12, y + 12), outline="red", width=2)
        draw.text((x + 12, y - 12), str(index), fill="red")
    OVERLAYS.mkdir(parents=True, exist_ok=True)
    image.save(OVERLAYS / f"page-{page}-{name}-candidates.png")


def make_point_overlay(page: int, name: str,
                       points: list[tuple[float, float]]) -> None:
    image = Image.open(RENDERS / f"page-{page}.jpg").convert("RGB")
    draw = ImageDraw.Draw(image)
    for index, (x, y) in enumerate(points, start=1):
        draw.ellipse((x - 12, y - 12, x + 12, y + 12), outline="blue", width=2)
        draw.text((x + 12, y - 12), str(index), fill="blue")
    OVERLAYS.mkdir(parents=True, exist_ok=True)
    image.save(OVERLAYS / f"page-{page}-{name}-candidates.png")


def guide_y(guide: list[tuple[int, int]], x: float) -> float:
    if x <= guide[0][0]:
        return float(guide[0][1])
    for (x0, y0), (x1, y1) in zip(guide, guide[1:], strict=True):
        if x <= x1:
            fraction = (x - x0) / (x1 - x0)
            return y0 + fraction * (y1 - y0)
    return float(guide[-1][1])


def selected_ct(page: int, points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    selected = []
    for x, y in points:
        _, value = calibrate(x, y, CT_AXES[page], 0.003, 0.013)
        # Cwt triangles and lower-curve features lie below this boundary.  The
        # isolated upper-left feature is the legend's hollow-square example.
        if value >= 0.00525 and not (x < 650 and y < 800):
            selected.append((x, y))
    return selected


def selected_cwp(page: int, candidates: list[tuple[int, int, int, int, int, int]]
                 ) -> list[tuple[float, float]]:
    selected = []
    for vertical, outer, _, tall, x, y in candidates:
        _, value = calibrate(x, y, CWP_AXES[page], 0.0, 0.010)
        if (
            vertical == 105
            and value <= CWP_VISUAL_MAX[page]
            and (outer >= 295 or tall >= 112)
        ):
            selected.append((float(x), float(y)))
    return sorted(selected)


def selected_trim(page: int, points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    return sorted(
        (x, y)
        for x, y in points
        if abs(y - guide_y(TRIM_GUIDES[page], x)) <= 50
    )


def selected_sinkage(page: int,
                      candidates: list[tuple[int, int, int, int, int, int]]
                      ) -> list[tuple[float, float]]:
    return sorted(
        (float(x), float(y))
        for _, center, _, _, x, y in candidates
        if center >= 70 and abs(y - guide_y(SINKAGE_GUIDES[page], x)) <= 42
    )


def write_pass(page: int, ct: list[tuple[float, float]],
               cwp: list[tuple[float, float]], trim: list[tuple[float, float]],
               sinkage: list[tuple[float, float]]) -> None:
    configuration, resistance_figure, attitude_figure, cwp_figure, printed_page, cwp_printed = SERIES[page]
    destination = OUTPUT / f"pass_a_free_{configuration}.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        handle.write("# source: Mustafa Insel, 1990 PhD thesis\n")
        handle.write(
            f"# source_location: printed pages {printed_page} and {cwp_printed}; "
            f"PDF pages {page} and {page + 1}; Figures {resistance_figure}, "
            f"{attitude_figure}, and {cwp_figure}\n"
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
                "value",
                "value_unit",
                "fn_digitization_uncertainty",
                "value_digitization_uncertainty",
                "source_figure",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        series = (
            ("ct", resistance_figure, CT_AXES[page], 0.003, 0.013,
             COEFFICIENT_UNCERTAINTY, "dimensionless", ct),
            ("cwp", cwp_figure, CWP_AXES[page + 1], 0.0, 0.010,
             COEFFICIENT_UNCERTAINTY, "dimensionless", cwp),
            ("trim", attitude_figure, ATTITUDE_AXES[page], -3.0, 6.0,
             TRIM_UNCERTAINTY_DEG, "degree", trim),
            ("sinkage_over_draught", attitude_figure, ATTITUDE_AXES[page],
             -0.1, 0.2, SINKAGE_UNCERTAINTY, "dimensionless", sinkage),
        )
        for observable, figure, axes, minimum, maximum, uncertainty, unit, points in series:
            for replicate_id, (x, y) in enumerate(points, start=1):
                fn, value = calibrate(x, y, axes, minimum, maximum)
                writer.writerow(
                    {
                        "configuration": configuration,
                        "attitude": "free",
                        "observable": observable,
                        "replicate_id": replicate_id,
                        "x_px": f"{x:.2f}",
                        "y_px": f"{y:.2f}",
                        "fn": f"{fn:.6f}",
                        "value": f"{value:.7f}",
                        "value_unit": unit,
                        "fn_digitization_uncertainty": f"{FN_UNCERTAINTY:.3f}",
                        "value_digitization_uncertainty": f"{uncertainty:.5f}",
                        "source_figure": figure,
                    }
                )


def main() -> None:
    for page in SERIES:
        page_image = np.asarray(Image.open(RENDERS / f"page-{page}.jpg").convert("L"))
        ct = selected_ct(
            page,
            enclosed_white_centers(page_image, CT_AXES[page], exclude_legend=True),
        )
        trim = selected_trim(
            page,
            enclosed_white_centers(page_image, ATTITUDE_AXES[page], exclude_legend=True),
        )
        sinkage = selected_sinkage(page, crossed_candidates(page_image, ATTITUDE_AXES[page]))
        cwp_page = page + 1
        image = np.asarray(Image.open(RENDERS / f"page-{cwp_page}.jpg").convert("L"))
        cwp = selected_cwp(cwp_page, dense_candidates(image, CWP_AXES[cwp_page]))
        make_point_overlay(page, "ct-selected", ct)
        make_point_overlay(page, "trim-selected", trim)
        make_point_overlay(page, "sinkage-selected", sinkage)
        make_point_overlay(cwp_page, "cwp-selected", cwp)
        write_pass(page, ct, cwp, trim, sinkage)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Reconcile blinded Insel C2 digitization passes without guessing a grid."""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
PASSES = HERE / "passes"

CONFIGURATIONS = {
    "monohull": "monohull",
    "SL0.2": "s_l_0_2",
    "SL0.3": "s_l_0_3",
    "SL0.4": "s_l_0_4",
    "SL0.5": "s_l_0_5",
    "S/L=0.2": "s_l_0_2",
    "S/L=0.3": "s_l_0_3",
    "S/L=0.4": "s_l_0_4",
    "S/L=0.5": "s_l_0_5",
}

FIGURE_TO_OBSERVABLE = {
    135: "ct", 136: "cwp", 137: "ct", 138: "cwp", 139: "ct",
    140: "cwp", 141: "ct", 142: "cwp", 143: "ct", 144: "cwp",
    161: "ct", 162: "attitude", 163: "cwp", 165: "ct",
    166: "attitude", 167: "cwp", 169: "ct", 170: "attitude",
    171: "cwp", 173: "ct", 174: "attitude", 175: "cwp",
    177: "ct", 178: "attitude", 179: "cwp",
}

BASE_UNCERTAINTY = {
    "ct": (0.002, 0.00008, "dimensionless"),
    "cwp": (0.002, 0.00008, "dimensionless"),
    "trim": (0.002, 0.05, "degree"),
    "sinkage_over_draught": (0.002, 0.002, "dimensionless"),
}


@dataclass
class Point:
    attitude: str
    configuration: str
    observable: str
    figure: int
    render_page: int
    printed_page: int
    x: float
    y: float
    pass_name: str


@dataclass
class Calibration:
    figure: int
    render_page: int
    printed_page: int
    x0: float
    x1: float
    y_top: float
    y_bottom: float


def normalize_configuration(raw: str) -> str:
    for token, normalized in CONFIGURATIONS.items():
        if token in raw:
            return normalized
    raise ValueError(f"unrecognized configuration {raw!r}")


def read_calibration() -> dict[int, Calibration]:
    result = {}
    with (PASSES / "pass_b_calibration.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            figure = int(row["figure"])
            result[figure] = Calibration(
                figure=figure,
                render_page=int(row["render_page"]),
                printed_page=int(row["printed_page"]),
                x0=float(row["pixel_x_left"]),
                x1=float(row["pixel_x_right"]),
                y_top=float(row["pixel_y_top"]),
                y_bottom=float(row["pixel_y_bottom"]),
            )
    return result


def read_pass_a(calibrations: dict[int, Calibration]) -> list[Point]:
    result = []
    for path in sorted(PASSES.glob("pass_a_*.csv")):
        with path.open() as handle:
            rows = csv.DictReader(line for line in handle if not line.startswith("#"))
            for row in rows:
                figure = int(row["source_figure"])
                calibration = calibrations[figure]
                result.append(
                    Point(
                        attitude=row["attitude"],
                        configuration=row["configuration"],
                        observable=row["observable"],
                        figure=figure,
                        render_page=calibration.render_page,
                        printed_page=calibration.printed_page,
                        x=float(row["x_px"]),
                        y=float(row["y_px"]),
                        pass_name="A",
                    )
                )
    return result


def read_tab_pass_b(*, include_fixed: bool, include_free: bool) -> list[Point]:
    result = []
    inputs = (
        ("pass_b_fixed_ct.tsv", "fixed", "ct"),
        ("pass_b_fixed_cwp.tsv", "fixed", "cwp"),
        ("pass_b_free_ct.tsv", "free", "ct"),
        ("pass_b_free_cwp.tsv", "free", "cwp"),
        ("pass_b_free_trim.tsv", "free", "trim"),
        ("pass_b_free_sinkage.tsv", "free", "sinkage_over_draught"),
    )
    for filename, attitude, observable in inputs:
        if (attitude == "fixed" and not include_fixed) or (attitude == "free" and not include_free):
            continue
        with (PASSES / filename).open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                result.append(
                    Point(
                        attitude=attitude,
                        configuration=normalize_configuration(row["configuration"]),
                        observable=observable,
                        figure=int(row["figure"]),
                        render_page=int(row["render_page"]),
                        printed_page=int(row["printed_page"]),
                        x=float(row["pixel_x"]),
                        y=float(row["pixel_y"]),
                        pass_name="B",
                    )
                )
    return result


def read_pass_c() -> list[Point]:
    result = []
    observable_names = {
        "CT": "ct",
        "CWP": "cwp",
        "trim_deg": "trim",
        "sinkage_over_draught": "sinkage_over_draught",
    }
    with (PASSES / "pass_c_free_attitude.csv").open() as handle:
        for row in csv.DictReader(handle):
            result.append(
                Point(
                    attitude="free",
                    configuration=normalize_configuration(row["config"]),
                    observable=observable_names[row["quantity"]],
                    figure=int(row["figure"]),
                    render_page=int(row["pdf_page"]),
                    printed_page=int(row["printed_page"]),
                    x=float(row["x_px"]),
                    y=float(row["y_px"]),
                    pass_name="C",
                )
            )
    return result


def read_independent_pass() -> list[Point]:
    """Return the conservative blind pass used to admit source markers."""
    return read_tab_pass_b(include_fixed=True, include_free=True)


def key(point: Point) -> tuple[str, str, str, int]:
    return point.attitude, point.configuration, point.observable, point.figure


def match_points(pass_a: list[Point], pass_b: list[Point]
                 ) -> tuple[list[tuple[Point | None, Point]], list[Point]]:
    by_key_a: dict[tuple[str, str, str, int], list[Point]] = defaultdict(list)
    for point in pass_a:
        by_key_a[key(point)].append(point)
    used: set[int] = set()
    matches = []
    for point_b in sorted(pass_b, key=lambda point: (key(point), point.x, point.y)):
        candidates = []
        for point_a in by_key_a[key(point_b)]:
            if id(point_a) in used:
                continue
            dx = abs(point_a.x - point_b.x)
            dy = abs(point_a.y - point_b.y)
            if dx <= 18 and dy <= 24:
                score = (dx / 12) ** 2 + (dy / 15) ** 2
                candidates.append((score, point_a))
        if candidates:
            _, point_a = min(candidates, key=lambda item: item[0])
            used.add(id(point_a))
            matches.append((point_a, point_b))
        else:
            matches.append((None, point_b))
    unmatched_a = [point for point in pass_a if id(point) not in used]
    return matches, unmatched_a


def calibrated(point: Point, calibration: Calibration) -> tuple[float, float]:
    fn = 0.1 + (point.x - calibration.x0) * 0.9 / (calibration.x1 - calibration.x0)
    if point.observable in ("ct", "cwp"):
        top, bottom = (0.013, 0.003) if point.observable == "ct" else (0.010, 0.0)
    elif point.observable == "trim":
        top, bottom = 6.0, -3.0
    elif point.observable == "sinkage_over_draught":
        top, bottom = 0.2, -0.1
    else:
        raise AssertionError(point.observable)
    value = bottom + (calibration.y_bottom - point.y) * (top - bottom) / (
        calibration.y_bottom - calibration.y_top
    )
    return fn, value


def write_outputs(matches: list[tuple[Point | None, Point]], unmatched_a: list[Point],
                  calibrations: dict[int, Calibration]) -> None:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    mismatch_rows = []
    for point_a, independent in matches:
        calibration = calibrations[independent.figure]
        independent_fn, independent_value = calibrated(independent, calibration)
        base_fn_uncertainty, base_value_uncertainty, unit = BASE_UNCERTAINTY[independent.observable]
        if point_a is None:
            x = independent.x
            y = independent.y
            fn = independent_fn
            value = independent_value
            fn_uncertainty = base_fn_uncertainty
            value_uncertainty = base_value_uncertainty
            pass_count = "1+source_recheck"
            adjudication = f"pass_a_missing; accepted conservative pass_{independent.pass_name.lower()} glyph"
            mismatch_rows.append(
                {
                    "attitude": independent.attitude,
                    "configuration": independent.configuration,
                    "observable": independent.observable,
                    "figure": independent.figure,
                    "mismatch": "pass_a_missing",
                    "resolution": f"accepted_pass_{independent.pass_name.lower()}_after_source_curve_recheck",
                    "x_px": f"{independent.x:.3f}",
                    "y_px": f"{independent.y:.3f}",
                }
            )
        else:
            fn_a, value_a = calibrated(point_a, calibration)
            x = (point_a.x + independent.x) / 2
            y = (point_a.y + independent.y) / 2
            averaged = Point(**{**independent.__dict__, "x": x, "y": y})
            fn, value = calibrated(averaged, calibration)
            fn_uncertainty = max(base_fn_uncertainty, abs(fn_a - independent_fn) / 2)
            value_uncertainty = max(base_value_uncertainty, abs(value_a - independent_value) / 2)
            pass_count = "2"
            adjudication = "independent pixel-center match; averaged"
        grouped[(independent.attitude, independent.configuration)].append(
            {
                "configuration": independent.configuration,
                "attitude": independent.attitude,
                "observable": independent.observable,
                "replicate_id": "",
                "fn": f"{fn:.7f}",
                "value": f"{value:.8f}",
                "value_unit": unit,
                "fn_digitization_uncertainty": f"{fn_uncertainty:.6f}",
                "value_digitization_uncertainty": f"{value_uncertainty:.8f}",
                "source_figure": str(independent.figure),
                "source_printed_page": str(independent.printed_page),
                "source_pdf_page": str(independent.render_page),
                "pass_count": pass_count,
                "adjudication": adjudication,
                "pass_a_x_px": "" if point_a is None else f"{point_a.x:.3f}",
                "pass_a_y_px": "" if point_a is None else f"{point_a.y:.3f}",
                "independent_pass": independent.pass_name,
                "independent_x_px": f"{independent.x:.3f}",
                "independent_y_px": f"{independent.y:.3f}",
            }
        )

    for point in unmatched_a:
        mismatch_rows.append(
            {
                "attitude": point.attitude,
                "configuration": point.configuration,
                "observable": point.observable,
                "figure": point.figure,
                "mismatch": "independent_pass_missing",
                "resolution": "excluded_conservatively_without_second_glyph_assignment",
                "x_px": f"{point.x:.3f}",
                "y_px": f"{point.y:.3f}",
            }
        )

    fields = [
        "configuration", "attitude", "observable", "replicate_id", "fn", "value",
        "value_unit", "fn_digitization_uncertainty", "value_digitization_uncertainty",
        "source_figure", "source_printed_page", "source_pdf_page", "pass_count",
        "adjudication", "pass_a_x_px", "pass_a_y_px", "independent_pass",
        "independent_x_px", "independent_y_px",
    ]
    for (attitude, configuration), rows in sorted(grouped.items()):
        rows.sort(key=lambda row: (row["observable"], float(row["fn"]), float(row["value"])))
        counters: dict[str, int] = defaultdict(int)
        for row in rows:
            counters[row["observable"]] += 1
            row["replicate_id"] = str(counters[row["observable"]])
        destination = HERE / f"{attitude}_{configuration}.csv"
        with destination.open("w", newline="") as handle:
            figures = sorted({row["source_figure"] for row in rows}, key=int)
            pages = sorted({row["source_printed_page"] for row in rows}, key=int)
            handle.write("# source: Mustafa Insel, 1990 PhD thesis\n")
            handle.write(f"# source_location: printed pages {','.join(pages)}; Figures {','.join(figures)}\n")
            handle.write("# method: blinded pixel passes, conservative glyph reconciliation; unsnapped Fn\n")
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    mismatch_fields = [
        "attitude", "configuration", "observable", "figure", "mismatch",
        "resolution", "x_px", "y_px",
    ]
    with (PASSES / "mismatches.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=mismatch_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            sorted(
                mismatch_rows,
                key=lambda row: (
                    row["attitude"], row["configuration"], row["observable"],
                    int(row["figure"]), float(row["x_px"]), float(row["y_px"]),
                ),
            )
        )


def write_free_audit(method_b: list[Point], pass_c: list[Point]) -> None:
    matches, unmatched_b = match_points(method_b, pass_c)
    counts: dict[tuple[str, str, str, int], dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    for point in method_b:
        counts[key(point)]["pass_b_rows"] += 1
    for method_point, point_c in matches:
        counts[key(point_c)]["pass_c_rows"] += 1
        if method_point is not None:
            counts[key(point_c)]["matched"] += 1
    for point in unmatched_b:
        counts[key(point)]["pass_b_only"] += 1
    rows = []
    for (_, configuration, observable, figure), count in sorted(counts.items()):
        rows.append(
            {
                "configuration": configuration,
                "observable": observable,
                "figure": figure,
                "pass_b_rows": count["pass_b_rows"],
                "pass_c_rows": count["pass_c_rows"],
                "matched": count["matched"],
                "pass_b_only": count["pass_b_only"],
                "pass_c_only": count["pass_c_rows"] - count["matched"],
            }
        )
    fields = [
        "configuration", "observable", "figure", "pass_b_rows", "pass_c_rows",
        "matched", "pass_b_only", "pass_c_only",
    ]
    with (PASSES / "free_pass_b_c_audit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    calibrations = read_calibration()
    pass_a = read_pass_a(calibrations)
    independent = read_independent_pass()
    matches, unmatched_a = match_points(pass_a, independent)
    write_outputs(matches, unmatched_a, calibrations)
    write_free_audit(
        read_tab_pass_b(include_fixed=False, include_free=True),
        read_pass_c(),
    )
    matched = sum(point_a is not None for point_a, _ in matches)
    print(f"independent rows: {len(matches)}; two-pass matches: {matched}; "
          f"independent-only accepted: {len(matches) - matched}; "
          f"A-only excluded: {len(unmatched_a)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Validate the frozen Insel C2 digitization archive and provenance fields."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent

EXPECTED_COUNTS = {
    "fixed_monohull.csv": {"ct": 21, "cwp": 31},
    "fixed_s_l_0_2.csv": {"ct": 25, "cwp": 18},
    "fixed_s_l_0_3.csv": {"ct": 40, "cwp": 34},
    "fixed_s_l_0_4.csv": {"ct": 40, "cwp": 45},
    "fixed_s_l_0_5.csv": {"ct": 30, "cwp": 34},
    "free_monohull.csv": {
        "ct": 22,
        "cwp": 27,
        "trim": 22,
        "sinkage_over_draught": 40,
    },
    "free_s_l_0_2.csv": {
        "ct": 26,
        "cwp": 16,
        "trim": 14,
        "sinkage_over_draught": 35,
    },
    "free_s_l_0_3.csv": {
        "ct": 25,
        "cwp": 14,
        "trim": 22,
        "sinkage_over_draught": 40,
    },
    "free_s_l_0_4.csv": {
        "ct": 44,
        "cwp": 26,
        "trim": 32,
        "sinkage_over_draught": 26,
    },
    "free_s_l_0_5.csv": {
        "ct": 40,
        "cwp": 28,
        "trim": 33,
        "sinkage_over_draught": 24,
    },
}

VALUE_BOUNDS = {
    "ct": (0.003, 0.013),
    "cwp": (0.0, 0.010),
    "trim": (-3.0, 6.0),
    "sinkage_over_draught": (-0.1, 0.2),
}

THEORY_COUNTS = {"359": 121, "360": 153, "361": 153, "362": 151}
THEORY_SEPARATIONS = {"359": 0.2, "360": 0.3, "361": 0.4, "362": 0.5}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def validate_file(path: Path) -> tuple[int, int]:
    rows = read_rows(path)
    counts = Counter(row["observable"] for row in rows)
    assert counts == Counter(EXPECTED_COUNTS[path.name]), (path.name, counts)

    expected_attitude = path.name.split("_", maxsplit=1)[0]
    by_observable: dict[str, list[dict[str, str]]] = defaultdict(list)
    matched = 0
    independent_only = 0
    for row in rows:
        observable = row["observable"]
        assert row["attitude"] == expected_attitude
        assert row["independent_pass"] == "B"
        assert 0.1 <= float(row["fn"]) <= 1.0
        lower, upper = VALUE_BOUNDS[observable]
        assert lower <= float(row["value"]) <= upper
        assert float(row["fn_digitization_uncertainty"]) > 0
        assert float(row["value_digitization_uncertainty"]) > 0
        assert row["independent_x_px"] and row["independent_y_px"]
        if row["pass_count"] == "2":
            assert row["pass_a_x_px"] and row["pass_a_y_px"]
            matched += 1
        else:
            assert row["pass_count"] == "1+source_recheck"
            assert not row["pass_a_x_px"] and not row["pass_a_y_px"]
            independent_only += 1
        by_observable[observable].append(row)

    for observable_rows in by_observable.values():
        assert [int(row["replicate_id"]) for row in observable_rows] == list(
            range(1, len(observable_rows) + 1)
        )
        ordering = [(float(row["fn"]), float(row["value"])) for row in observable_rows]
        assert ordering == sorted(ordering)
    return matched, independent_only


def validate_theory() -> tuple[int, int]:
    rows = read_rows(HERE / "theory_interference.csv")
    counts = Counter(row["figure"] for row in rows)
    assert counts == Counter(THEORY_COUNTS), counts

    by_figure: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        figure = row["figure"]
        assert float(row["separation_over_length"]) == THEORY_SEPARATIONS[figure]
        assert row["source_printed_page"] in {"358", "359"}
        assert int(row["source_pdf_page"]) == int(row["source_printed_page"]) + 10
        assert row["pass_count"] == "2"
        assert 0.15 <= float(row["fn"]) <= 0.95
        assert 0.0 <= float(row["tau"]) <= 2.5
        assert float(row["fn_pass_difference"]) <= 0.015
        assert float(row["tau_pass_difference"]) <= 0.08
        assert float(row["fn_digitization_uncertainty"]) >= 0.003
        assert float(row["tau_digitization_uncertainty"]) >= 0.02
        assert abs(float(row["fn"]) - float(row["fn_anchor"])) < 1e-6
        assert all(row[field] for field in (
            "pass_a_x_px", "pass_a_y_px", "pass_b_x_px", "pass_b_y_px"
        ))
        by_figure[figure].append(row)

    for figure, figure_rows in by_figure.items():
        anchors = [float(row["fn_anchor"]) for row in figure_rows]
        assert anchors == sorted(anchors)
        scoring_anchors = [
            anchor for anchor in anchors
            if 0.20 <= anchor <= 0.80 and round(anchor * 1000) % 10 == 0
        ]
        assert len(scoring_anchors) >= 50, (figure, len(scoring_anchors))
        assert min(anchors) <= 0.35 and max(anchors) >= 0.55

    mismatch_rows = read_rows(HERE / "passes/theory_mismatches.csv")
    assert len(mismatch_rows) == 36
    assert all(
        row["resolution"] == "omitted_after_source_only_reinspection"
        and float(row["tau_pass_difference"]) > 0.08
        for row in mismatch_rows
    )
    return len(rows), len(mismatch_rows)


def main() -> None:
    matched = 0
    independent_only = 0
    for filename in EXPECTED_COUNTS:
        file_matched, file_independent_only = validate_file(HERE / filename)
        matched += file_matched
        independent_only += file_independent_only

    mismatch_rows = read_rows(HERE / "passes/mismatches.csv")
    resolutions = Counter(row["mismatch"] for row in mismatch_rows)
    assert resolutions == {
        "independent_pass_missing": 525,
        "pass_a_missing": 78,
    }
    assert matched == 796
    assert independent_only == 78
    theory_admitted, theory_omitted = validate_theory()
    print(
        "validated 874 admitted markers: "
        f"{matched} two-pass matches, {independent_only} conservative pass-B-only markers; "
        "525 pass-A-only candidates excluded"
    )
    print(
        f"validated {theory_admitted} two-pass theory anchors; "
        f"{theory_omitted} source-rechecked anchors omitted"
    )


if __name__ == "__main__":
    main()

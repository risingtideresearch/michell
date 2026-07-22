"""Round-trip the ``.msw`` sweep-archive reader against a hand-built archive.

Encodes a small archive in Python exactly per the format in
``crates/michell-cli/src/archive.rs`` and checks :func:`pymichell.read_sweep`
recovers the manifest, hull files, metadata, and every per-row payload
(parameters, metrics, GZ curve, spectrum).
"""

import json
import struct

import numpy as np

from pymichell import read_sweep


def _blob(kind, name, data):
    name = name.encode("utf-8")
    return (
        struct.pack("<I", kind)
        + struct.pack("<I", len(name))
        + name
        + struct.pack("<Q", len(data))
        + data
    )


def _build_archive():
    meta = {
        "axis_labels": ["speed_kn", "vaka:mass"],
        "metric_labels": ["rw", "rt"],
        "speeds_ms": [3.0],
        "name": "unit",
    }
    # ROWS blob: 1 row, 2 axes, 2 metrics.
    rows = struct.pack("<III", 1, 2, 2)
    rows += struct.pack("<2d", 6.0, 1500.0)  # params
    rows += struct.pack("<2d", 12.5, 40.0)  # metrics
    gz = [(0.0, 0.0), (0.1, 0.4), (0.2, -0.1)]
    rows += struct.pack("<I", len(gz))
    for h, g in gz:
        rows += struct.pack("<2d", h, g)
    rows += struct.pack("<2d", 1.09, 5.75)  # wavenumber, transverse wavelength
    spec = [(-0.1, 1.0, -2.0, 3.0), (0.0, 4.0, 0.0, 5.0), (0.1, 6.0, 1.0, 7.0)]
    rows += struct.pack("<I", len(spec))
    for th, re, im, d in spec:
        rows += struct.pack("<4d", th, re, im, d)

    body = b"MSWP" + struct.pack("<I", 1)
    body += _blob(1, "study.json", b'{"name":"unit"}')
    body += _blob(2, "vaka.hull", b"michell-hull v1\n")
    body += _blob(3, "meta.json", json.dumps(meta).encode("utf-8"))
    body += _blob(4, "rows", rows)
    return body, meta


def test_read_sweep_round_trip():
    data, meta = _build_archive()
    sw = read_sweep(data)

    assert sw.manifest_name == "study.json"
    assert sw.manifest == '{"name":"unit"}'
    assert sw.hull_files == {"vaka.hull": b"michell-hull v1\n"}
    assert sw.meta == meta
    assert sw.axis_labels == ["speed_kn", "vaka:mass"]

    assert len(sw.rows) == 1
    row = sw.rows[0]
    assert row.params == {"speed_kn": 6.0, "vaka:mass": 1500.0}
    assert row.metrics == {"rw": 12.5, "rt": 40.0}

    assert row.gz_curve.shape == (3, 2)
    np.testing.assert_allclose(row.gz_curve[1], [0.1, 0.4])

    sp = row.spectrum
    assert sp.wavenumber == 1.09
    assert sp.transverse_wavelength == 5.75
    assert sp.theta.shape == (3,)
    np.testing.assert_allclose(sp.theta, [-0.1, 0.0, 0.1])
    np.testing.assert_allclose(sp.amp, [1 - 2j, 4 + 0j, 6 + 1j])
    np.testing.assert_allclose(sp.drw_dtheta, [3.0, 5.0, 7.0])


def test_rejects_bad_magic():
    try:
        read_sweep(b"NOPE" + b"\x00" * 8)
    except ValueError:
        return
    raise AssertionError("expected a ValueError on bad magic")

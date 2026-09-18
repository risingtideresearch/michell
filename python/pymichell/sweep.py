"""Reader for ``michell sweep`` binary archives (``.msw``).

A ``.msw`` file bundles a whole study into one self-contained container: the
JSON manifest, the referenced hull files (verbatim), and, for every output
row, the swept parameter values, the scalar metrics, the full righting-arm
(GZ) curve, and the free-wave spectrum ``A(θ)``. See ``crates/michell-cli/
src/archive.rs`` for the byte-level format; this reader mirrors it.

Because each row carries its spectrum, a stored sweep is enough to regenerate a
wake elevation field or heatmap on demand, at any resolution, without re-running
the study::

    from pymichell import read_sweep
    sweep = read_sweep("study.msw")
    row = sweep.rows[0]
    theta, amp = row.spectrum.theta, row.spectrum.amp   # A(θ) as complex

All integers and floats are little-endian.
"""

from __future__ import annotations

import json
import struct
from dataclasses import dataclass

import numpy as np

_MAGIC = b"MSWP"

# Blob kinds (see archive.rs).
KIND_MANIFEST = 1
KIND_HULLFILE = 2
KIND_META = 3
KIND_ROWS = 4


@dataclass
class Spectrum:
    """The free-wave spectrum stored with a row."""

    wavenumber: float
    transverse_wavelength: float
    theta: np.ndarray  # (n,) propagation angle [rad]
    amp: np.ndarray  # (n,) complex free-wave amplitude A(θ)
    drw_dtheta: np.ndarray  # (n,) resistance density dR_w/dθ


@dataclass
class Row:
    """One output row: parameters, metrics, GZ curve, and spectrum."""

    params: dict  # axis label -> swept value
    metrics: dict  # metric label -> value
    gz_curve: np.ndarray  # (m, 2) columns (heel_rad, gz_m); empty when undefined
    spectrum: Spectrum


@dataclass
class Sweep:
    """A parsed sweep archive."""

    manifest_name: str
    manifest: str  # raw manifest JSON text
    hull_files: dict  # filename -> raw bytes
    meta: dict  # parsed meta.json
    rows: list  # list[Row]

    @property
    def axis_labels(self):
        return self.meta.get("axis_labels", [])

    @property
    def metric_labels(self):
        return self.meta.get("metric_labels", [])


class _Cursor:
    __slots__ = ("buf", "pos")

    def __init__(self, buf, pos=0):
        self.buf = buf
        self.pos = pos

    def take(self, n):
        b = self.buf[self.pos : self.pos + n]
        self.pos += n
        return b

    def u32(self):
        return struct.unpack_from("<I", self.buf, self._adv(4))[0]

    def u64(self):
        return struct.unpack_from("<Q", self.buf, self._adv(8))[0]

    def f64_array(self, n):
        a = np.frombuffer(self.buf, dtype="<f8", count=n, offset=self.pos)
        self.pos += 8 * n
        return a.astype(np.float64)

    def _adv(self, n):
        p = self.pos
        self.pos += n
        return p

    def at_end(self):
        return self.pos >= len(self.buf)


def read_sweep(path) -> Sweep:
    """Parse a ``.msw`` archive from a path or bytes."""
    if isinstance(path, (bytes, bytearray)):
        data = bytes(path)
    else:
        with open(path, "rb") as fh:
            data = fh.read()

    if data[:4] != _MAGIC:
        raise ValueError("not an MSWP archive (bad magic)")
    version = struct.unpack_from("<I", data, 4)[0]
    if version != 1:
        raise ValueError(f"unsupported .msw version {version}")

    manifest_name = "manifest.json"
    manifest = ""
    hull_files: dict = {}
    meta: dict = {}
    rows: list = []

    cur = _Cursor(data, 8)
    while not cur.at_end():
        kind = cur.u32()
        name_len = cur.u32()
        name = cur.take(name_len).decode("utf-8")
        data_len = cur.u64()
        blob = cur.take(data_len)
        if kind == KIND_MANIFEST:
            manifest_name = name
            manifest = blob.decode("utf-8")
        elif kind == KIND_HULLFILE:
            hull_files[name] = bytes(blob)
        elif kind == KIND_META:
            meta = json.loads(blob.decode("utf-8"))
        elif kind == KIND_ROWS:
            rows = _decode_rows(blob, meta)
        # Unknown kinds are ignored for forward compatibility.

    return Sweep(
        manifest_name=manifest_name,
        manifest=manifest,
        hull_files=hull_files,
        meta=meta,
        rows=rows,
    )


def _decode_rows(blob: bytes, meta: dict) -> list:
    cur = _Cursor(blob)
    n_rows = cur.u32()
    n_axes = cur.u32()
    n_metrics = cur.u32()
    axis_labels = meta.get("axis_labels") or [f"axis{i}" for i in range(n_axes)]
    metric_labels = meta.get("metric_labels") or [f"metric{i}" for i in range(n_metrics)]

    rows = []
    for _ in range(n_rows):
        params_vals = cur.f64_array(n_axes)
        metric_vals = cur.f64_array(n_metrics)
        params = dict(zip(axis_labels, params_vals.tolist()))
        metrics = dict(zip(metric_labels, metric_vals.tolist()))

        gz_n = cur.u32()
        gz = cur.f64_array(2 * gz_n).reshape(gz_n, 2) if gz_n else np.empty((0, 2))

        wavenumber = cur.f64_array(1)[0]
        twl = cur.f64_array(1)[0]
        spec_n = cur.u32()
        flat = cur.f64_array(4 * spec_n).reshape(spec_n, 4) if spec_n else np.empty((0, 4))
        spectrum = Spectrum(
            wavenumber=float(wavenumber),
            transverse_wavelength=float(twl),
            theta=flat[:, 0].copy(),
            amp=(flat[:, 1] + 1j * flat[:, 2]),
            drw_dtheta=flat[:, 3].copy(),
        )
        rows.append(Row(params=params, metrics=metrics, gz_curve=gz, spectrum=spectrum))
    return rows

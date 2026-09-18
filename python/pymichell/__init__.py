"""pymichell -- a small, clear reproduction of the Michell wave-amplitude grid.

A pedagogical companion to the Rust ``michell`` crate.  Given a hull as a
B-spline half-breadth control net, it computes the thin-ship free-wave
amplitude ``F = I + iJ``, the far-field spectrum ``A(θ)``, the Kelvin-wake
elevation grid ``ζ(x, y)``, and the wave resistance ``R_w``.

Like the crate, the inner amplitude integral is evaluated **exactly**,
span by span: on each knot span ``∂f/∂x`` is a local polynomial and the
integrals reduce to closed-form moments (:mod:`pymichell.moments`), so
``I`` and ``J`` match the crate to machine precision.  Only the smooth outer
θ-integrals use plain trapezoidal quadrature, for legibility.

Quick start::

    import numpy as np
    from pymichell import wigley, wave_field

    surface = wigley(10.0, 1.0, 0.625)      # L, B, T [m]  (or load_hull("boat.hull"))
    field = wave_field(surface, speed=3.0)   # U [m/s], seawater

    print(field.wave_resistance(), "N")

    x = np.linspace(surface.x_center - 35, surface.x_center + 8, 240)
    y = np.linspace(-18, 18, 200)
    zeta = field.elevation_grid(x, y)        # (len(y), len(x)) elevation grid [m]
"""

from .bspline import BSplineSurface, load_hull, wigley
from .sweep import Row, Spectrum, Sweep, read_sweep
from .wave import (
    Conditions,
    FRESHWATER_DENSITY,
    SEAWATER_DENSITY,
    STANDARD_GRAVITY,
    WaveField,
    wave_field,
)

__all__ = [
    "BSplineSurface",
    "load_hull",
    "wigley",
    "Conditions",
    "WaveField",
    "wave_field",
    "read_sweep",
    "Sweep",
    "Row",
    "Spectrum",
    "SEAWATER_DENSITY",
    "FRESHWATER_DENSITY",
    "STANDARD_GRAVITY",
]

__version__ = "0.1.0"

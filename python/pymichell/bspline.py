"""Clamped, polynomial tensor-product B-spline half-breadth surfaces.

A hull is described by ``y = f(x, z)``, the local half-beam, as a
tensor-product B-spline surface

    f(x, z) = Σ_i Σ_j  P[i, j] · N_i(x) · M_j(z)

with ``P`` the *control net*, ``N_i`` the degree-``p`` basis functions over
the ``x`` knot vector, and ``M_j`` the degree-``q`` basis over ``z``.  This
is exactly the representation the Rust ``michell`` crate loads from a
``.hull`` file; here we reproduce just enough of it to feed the wave
calculation.

Conventions (identical to the crate):

* ``x`` runs along the hull, ``z`` runs vertically **downward** from the
  undisturbed waterline (``z = 0``), domain ``z ∈ [0, T]``.  SI units.
* Knot vectors are *clamped* (end knots repeated ``degree + 1`` times).
* The surface is polynomial (all NURBS weights ``= 1``).

The basis functions are evaluated by the Cox--de Boor recurrence rather
than any of the fast closed-form machinery in the crate: the point of this
package is clarity, not speed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _basis_derivatives(knots: np.ndarray, degree: int, u: np.ndarray, max_deriv: int) -> list:
    """All B-spline basis functions and derivatives up to ``max_deriv``.

    Returns a list ``D`` of length ``max_deriv + 1``; ``D[k]`` is a matrix of
    shape ``(len(u), n_ctrl)`` whose entry ``D[k][r, i]`` is the ``k``-th
    derivative ``N^{(k)}_{i,degree}(u[r])``.  ``n_ctrl = len(knots) - degree - 1``.

    Built straight from the textbook recurrences, raising the degree one
    level at a time (Cox--de Boor for the value, and the standard derivative
    rule -- Piegl & Tiller Eq. 3.10 -- which expresses a degree-``p``
    derivative in terms of degree-``p-1`` derivatives)::

        N_{i,0}(u)      = 1 if t_i <= u < t_{i+1} else 0
        N_{i,p}(u)      = (u - t_i)/(t_{i+p} - t_i)         · N_{i,p-1}(u)
                        + (t_{i+p+1} - u)/(t_{i+p+1} - t_{i+1}) · N_{i+1,p-1}(u)
        N^{(k)}_{i,p}(u) = p [ N^{(k-1)}_{i,p-1}(u)/(t_{i+p} - t_i)
                             - N^{(k-1)}_{i+1,p-1}(u)/(t_{i+p+1} - t_{i+1}) ]

    Points outside the parametric domain are clamped onto it (matching the
    crate's ``eval``); the right domain endpoint is folded into the last
    non-empty span so it evaluates there rather than falling through the
    half-open interval.
    """
    max_deriv = min(max_deriv, degree)
    knots = np.asarray(knots, dtype=float)
    n_ctrl = len(knots) - degree - 1
    u = np.clip(np.atleast_1d(np.asarray(u, dtype=float)), knots[degree], knots[n_ctrl])

    # Degree 0: indicator functions of each half-open knot span (only the
    # 0-th derivative is non-zero at this level).
    basis = ((u[:, None] >= knots[:-1]) & (u[:, None] < knots[1:])).astype(float)
    last = n_ctrl - 1
    while knots[last + 1] <= knots[last]:
        last -= 1
    at_end = u >= knots[n_ctrl]
    basis[at_end, :] = 0.0
    basis[at_end, last] = 1.0
    levels = [basis]  # derivatives 0.. of the current degree

    for p in range(1, degree + 1):
        n_p = len(knots) - p - 1
        new_levels = [np.zeros((len(u), n_p)) for _ in range(min(max_deriv, p) + 1)]
        for i in range(n_p):
            left_den = knots[i + p] - knots[i]
            right_den = knots[i + p + 1] - knots[i + 1]
            for k in range(len(new_levels)):
                if k == 0:  # value recurrence
                    if left_den > 0:
                        new_levels[0][:, i] += (u - knots[i]) / left_den * levels[0][:, i]
                    if right_den > 0:
                        new_levels[0][:, i] += (knots[i + p + 1] - u) / right_den * levels[0][:, i + 1]
                else:  # derivative recurrence, from the level below
                    if left_den > 0:
                        new_levels[k][:, i] += p * levels[k - 1][:, i] / left_den
                    if right_den > 0:
                        new_levels[k][:, i] -= p * levels[k - 1][:, i + 1] / right_den
        levels = new_levels

    # Pad with zero matrices if max_deriv exceeds the degree.
    while len(levels) <= max_deriv:
        levels.append(np.zeros((len(u), n_ctrl)))
    return levels


def _basis_matrix(knots: np.ndarray, degree: int, u: np.ndarray, deriv: int = 0) -> np.ndarray:
    """The ``deriv``-th derivative basis matrix (shape ``(len(u), n_ctrl)``)."""
    return _basis_derivatives(knots, degree, u, deriv)[deriv]


def _nonempty_spans(knots: np.ndarray, degree: int) -> list:
    """Non-empty knot spans as ``(knot_index, start, length)`` triples.

    A clamped knot vector has ``n_ctrl - degree`` candidate spans
    ``[knots[s], knots[s+1])`` for ``s = degree .. n_ctrl - 1``; repeated
    (interior) knots make some of them empty, and those are dropped.
    """
    knots = np.asarray(knots, dtype=float)
    n_ctrl = len(knots) - degree - 1
    spans = []
    for s in range(degree, n_ctrl):
        if knots[s + 1] > knots[s]:
            spans.append((s, float(knots[s]), float(knots[s + 1] - knots[s])))
    return spans


@dataclass
class BSplineSurface:
    """A validated clamped polynomial tensor-product B-spline surface.

    ``control`` has shape ``(n_ctrl_x, n_ctrl_z)`` (``x`` index first,
    ``z`` index second), matching the row-major, ``z``-fastest layout of the
    ``.hull`` file.

    ``waterline`` is the design waterline's depth below ``z = 0``.  For a
    plain wetted-surface hull (the usual case -- ``z = 0`` *is* the
    waterline) it is ``0`` and the whole ``z`` domain is submerged.  For a
    full-band *body* file (which models the hull up above the waterline) it
    is the positive depth at which the water surface cuts the band; the
    submerged region is then ``z ∈ [waterline, draft]``.
    """

    degree_x: int
    degree_z: int
    knots_x: np.ndarray
    knots_z: np.ndarray
    control: np.ndarray
    waterline: float = 0.0

    def __post_init__(self) -> None:
        self.knots_x = np.asarray(self.knots_x, dtype=float)
        self.knots_z = np.asarray(self.knots_z, dtype=float)
        self.control = np.asarray(self.control, dtype=float)
        n_ctrl_x = len(self.knots_x) - self.degree_x - 1
        n_ctrl_z = len(self.knots_z) - self.degree_z - 1
        if self.control.shape != (n_ctrl_x, n_ctrl_z):
            raise ValueError(
                f"control net has shape {self.control.shape}, expected "
                f"({n_ctrl_x}, {n_ctrl_z}) from the knot vectors"
            )

    # -- domain -----------------------------------------------------------
    @property
    def x_domain(self) -> tuple[float, float]:
        return float(self.knots_x[self.degree_x]), float(self.knots_x[len(self.knots_x) - self.degree_x - 1])

    @property
    def z_domain(self) -> tuple[float, float]:
        return float(self.knots_z[self.degree_z]), float(self.knots_z[len(self.knots_z) - self.degree_z - 1])

    @property
    def x_center(self) -> float:
        x0, x1 = self.x_domain
        return 0.5 * (x0 + x1)

    @property
    def length(self) -> float:
        x0, x1 = self.x_domain
        return x1 - x0

    @property
    def draft(self) -> float:
        return self.z_domain[1]

    # -- evaluation -------------------------------------------------------
    def evaluate(self, x, z, dx: int = 0, dz: int = 0) -> np.ndarray:
        """Evaluate ``∂^{dx+dz} f / ∂x^{dx} ∂z^{dz}`` on the grid ``x × z``.

        ``x`` and ``z`` are 1-D arrays of coordinates; the result has shape
        ``(len(x), len(z))``.  With ``dx = dz = 0`` this is the half-beam
        surface itself; ``dx = 1`` gives ``∂f/∂x``, the slope Michell's
        integral needs.
        """
        bx = _basis_matrix(self.knots_x, self.degree_x, np.atleast_1d(x), dx)
        bz = _basis_matrix(self.knots_z, self.degree_z, np.atleast_1d(z), dz)
        # errstate guards a spurious FP-flag warning some BLAS backends raise
        # from matmul (notably NumPy 2.0 on Apple Accelerate); the result is
        # finite and correct.
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            return bx @ self.control @ bz.T

    # -- knot spans (the pieces the closed-form integral works over) ------
    def x_spans(self) -> list:
        """Non-empty knot spans in x as ``(knot_index, start, length)``."""
        return _nonempty_spans(self.knots_x, self.degree_x)

    def z_spans(self) -> list:
        """Non-empty knot spans in z as ``(knot_index, start, length)``."""
        return _nonempty_spans(self.knots_z, self.degree_z)

    def corner_partials(self, span_x: int, span_z: int) -> np.ndarray:
        """Mixed partials at a span's lower-left corner.

        Returns ``D`` of shape ``(degree_x + 1, degree_z + 1)`` with
        ``D[a, b] = ∂^{a+b} f / ∂x^a ∂z^b`` evaluated at the corner
        ``(knots_x[span_x], knots_z[span_z])``.  On the span rectangle the
        surface is exactly polynomial, so together with Taylor's theorem
        these partials *are* the local polynomial:
        ``f = Σ_{a,b} D[a, b]/(a! b!) · (x - x0)^a (z - z0)^b``.

        For a tensor-product surface the mixed partial factorises into the
        x- and z-basis derivatives at the corner::

            D[a, b] = Σ_i Σ_j control[i, j] · N^{(a)}_i(x0) · M^{(b)}_j(z0)
        """
        x0 = self.knots_x[span_x]
        z0 = self.knots_z[span_z]
        dx = _basis_derivatives(self.knots_x, self.degree_x, np.array([x0]), self.degree_x)
        dz = _basis_derivatives(self.knots_z, self.degree_z, np.array([z0]), self.degree_z)
        nx = np.vstack([d[0] for d in dx])  # (degree_x + 1, n_ctrl_x)
        mz = np.vstack([d[0] for d in dz])  # (degree_z + 1, n_ctrl_z)
        return nx @ self.control @ mz.T


def load_hull(path: str) -> BSplineSurface:
    """Parse a ``michell-hull v1`` ``.hull`` file into a :class:`BSplineSurface`.

    The format (SI units, ``#`` comments) is::

        michell-hull v1
        degree-x 2
        degree-z 2
        knots-x -5 -5 -5 5 5 5
        knots-z 0 0 0 0.625 0.625 0.625
        row 0 0 0        # one line per x control index, n_ctrl_z values each
        row 1 1 0
        row 0 0 0

    ``waterline`` and ``centerplane`` lines (written for full-band hull
    bodies) are accepted and ignored -- they are not needed for the
    wave-amplitude calculation.
    """
    degree_x = degree_z = None
    knots_x = knots_z = None
    waterline = 0.0
    rows: list[list[float]] = []
    header_seen = False

    with open(path, "r") as handle:
        for raw in handle:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            key, *rest = line.split()
            if key == "michell-hull":
                header_seen = True
            elif key == "degree-x":
                degree_x = int(rest[0])
            elif key == "degree-z":
                degree_z = int(rest[0])
            elif key == "knots-x":
                knots_x = [float(v) for v in rest]
            elif key == "knots-z":
                knots_z = [float(v) for v in rest]
            elif key == "row":
                rows.append([float(v) for v in rest])
            elif key == "waterline":
                waterline = float(rest[0])  # design WL depth below the band top
            elif key == "centerplane":
                continue  # transverse position; irrelevant to a single-hull wave field
            else:
                raise ValueError(f"{path}: unrecognised line: {line!r}")

    if not header_seen:
        raise ValueError(f"{path}: missing 'michell-hull' header")
    if degree_x is None or degree_z is None or knots_x is None or knots_z is None:
        raise ValueError(f"{path}: missing one of degree-x, degree-z, knots-x, knots-z")

    control = np.array(rows, dtype=float)
    return BSplineSurface(
        degree_x, degree_z, np.array(knots_x), np.array(knots_z), control, waterline
    )


def wigley(length: float = 10.0, beam: float = 1.0, draft: float = 0.625) -> BSplineSurface:
    """The standard Wigley parabolic hull as an exact biquadratic surface.

    ``f(x, z) = (B/2)(1 - (2x/L)^2)(1 - (z/T)^2)`` on ``x ∈ [-L/2, L/2]``,
    ``z ∈ [0, T]``.  Handy as a self-contained reference with a known
    analytic wave resistance.
    """
    a = length / 2.0
    knots_x = [-a, -a, -a, a, a, a]
    knots_z = [0.0, 0.0, 0.0, draft, draft, draft]
    gx = np.array([0.0, 2.0, 0.0])  # quadratic Bezier weights for 1 - (2x/L)^2
    hz = np.array([1.0, 1.0, 0.0])  # ... and for 1 - (z/T)^2
    control = (beam / 2.0) * np.outer(gx, hz)
    return BSplineSurface(2, 2, np.array(knots_x), np.array(knots_z), control)

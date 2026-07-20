"""Closed-form moment integrals against oscillatory and exponential kernels.

These are the primitives that make the Michell inner integrals *exact* for a
piecewise-polynomial (B-spline) hull.  On every knot span the integrand is a
polynomial times ``e^{ikx}`` in x and a polynomial times ``e^{-κz}`` in z, and
both families of moments below have closed forms -- so ``I`` and ``J`` carry
no quadrature error at all.

Two evaluation branches (a direct port of the Rust ``moments.rs``):

* **Series** for small ``|k|h`` -- the entire Taylor series, which avoids the
  catastrophic cancellation the recurrence would suffer near ``k = 0``.
* **Recurrence** for large ``|k|h`` -- a stable upward recurrence in the
  polynomial power, exact up to round-off.

Each function is vectorized over the kernel rate (an array of ``k`` or
``κ``, one per propagation angle ``λ = sec θ``); ``h`` is one span length.
"""

from __future__ import annotations

import numpy as np

# Switch between the power series (small argument) and the closed-form
# recurrence (large argument), matching the crate.
SERIES_THRESHOLD = 4.0
# Series term count: at the threshold |k|h = 4 the terms 4^m/m! have fallen
# below 1e-24 by m ~ 40, so 64 reaches machine precision with room to spare.
_SERIES_TERMS = 64


def osc_moments(k: np.ndarray, h: float, a_max: int) -> np.ndarray:
    """Oscillatory moments ``M_a = ∫_0^h t^a e^{i k t} dt`` for ``a = 0..a_max``.

    ``k`` is an array (one wavenumber per angle); returns a complex array of
    shape ``(len(k), a_max + 1)`` with ``M[:, a] = M_a``.

    * Series: ``M_a = h^{a+1} Σ_m (ikh)^m / (m! (a + m + 1))``.
    * Recurrence: ``M_a = (h^a e^{ikh} - a M_{a-1}) / (ik)``.
    """
    k = np.atleast_1d(np.asarray(k, dtype=float))
    moments = np.zeros((len(k), a_max + 1), dtype=complex)
    small = np.abs(k * h) <= SERIES_THRESHOLD

    if small.any():
        ikh = 1j * k[small] * h
        for a in range(a_max + 1):
            term = np.ones(np.count_nonzero(small), dtype=complex)  # (ikh)^m / m!
            series = np.zeros_like(term)
            for m in range(_SERIES_TERMS):
                series += term / (a + m + 1)
                term = term * ikh / (m + 1)
            moments[small, a] = series * h ** (a + 1)

    if (~small).any():
        kb = k[~small]
        e = np.exp(1j * kb * h)
        inv_ik = 1.0 / (1j * kb)
        prev = (e - 1.0) * inv_ik
        moments[~small, 0] = prev
        h_pow = h
        for a in range(1, a_max + 1):
            prev = (e * h_pow - a * prev) * inv_ik
            moments[~small, a] = prev
            h_pow *= h

    return moments


def exp_moments(kappa: np.ndarray, h: float, b_max: int) -> np.ndarray:
    """Exponential moments ``N_b = ∫_0^h t^b e^{-κ t} dt`` for ``b = 0..b_max``.

    ``κ ≥ 0`` is an array (one decay rate per angle); returns a real array of
    shape ``(len(κ), b_max + 1)`` with ``N[:, b] = N_b``.

    * Series: ``N_b = h^{b+1} Σ_m (-κh)^m / (m! (b + m + 1))``.
    * Recurrence: ``N_b = (b N_{b-1} - h^b e^{-κh}) / κ``.
    """
    kappa = np.atleast_1d(np.asarray(kappa, dtype=float))
    moments = np.zeros((len(kappa), b_max + 1), dtype=float)
    x = kappa * h
    small = x <= SERIES_THRESHOLD

    if small.any():
        neg_x = -x[small]
        for b in range(b_max + 1):
            term = np.ones(np.count_nonzero(small), dtype=float)  # (-κh)^m / m!
            series = np.zeros_like(term)
            for m in range(_SERIES_TERMS):
                series += term / (b + m + 1)
                term = term * neg_x / (m + 1)
            moments[small, b] = series * h ** (b + 1)

    if (~small).any():
        kb = kappa[~small]
        e = np.exp(-x[~small])
        prev = (1.0 - e) / kb
        moments[~small, 0] = prev
        h_pow = h
        for b in range(1, b_max + 1):
            prev = (b * prev - h_pow * e) / kb
            moments[~small, b] = prev
            h_pow *= h

    return moments

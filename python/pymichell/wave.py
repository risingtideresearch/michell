"""Michell thin-ship wave amplitude, free-wave spectrum, and Kelvin wake.

This is the heart of the package: from a B-spline half-breadth surface it
builds the far-field wave pattern a slender hull leaves behind, and the wave
resistance that pattern carries away.

The physics (Tuck 1989; Tuck, Scullen & Lazauskas 2001/2002; Dambrine,
Pierre & Rousseaux 2016), with ``ν = g/U²`` and the half-beam ``f(x, z)``:

1. **Amplitude function.**  Each propagation angle ``θ ∈ (-π/2, π/2)`` off
   the ship's track carries one free plane wave.  Its complex weight comes
   from the Michell inner integral evaluated at ``λ = sec θ``::

       F(λ) = I(λ) + i J(λ)
            = ∬ (∂f/∂x) · exp(-ν λ² z) · exp(i ν λ (x - x_c)) dx dz

   (phases referenced to the hull mid-length ``x_c``; ``|F|`` does not
   depend on that choice).  Because the hull is piecewise polynomial this
   integral is evaluated **exactly, span by span**: on each knot span
   ``∂f/∂x`` is a local polynomial (:meth:`BSplineSurface.corner_partials`)
   and the two 1-D integrals reduce to the closed-form moments in
   :mod:`pymichell.moments`.  ``I`` and ``J`` carry no quadrature error.

2. **Free-wave amplitude density.**

       A(θ) = -(2ν/π) sec³θ · conj(F(sec θ))     [m/rad]

   whose magnitude is pinned by the deep-water free-wave resistance
   identity below and whose phase follows Tuck, Scullen & Lazauskas.

3. **Wave resistance.**  Two equivalent forms, both computed here as a
   check on each other:

       R_w = (4 ρ g²)/(π U²) ∫₀^{π/2} |F(sec θ)|² sec³θ dθ          (Michell)
           = ½ π ρ U²        ∫_{-π/2}^{π/2} |A(θ)|² cos³θ dθ        (spectrum)

4. **Kelvin wake.**  The far-field wave elevation is the superposition of
   all those plane waves::

       ζ(x, y) = Re ∫_{-π/2}^{π/2} A(θ) · exp(i ν sec θ ((x - x_c) + y tan θ)) dθ

   The ship advances toward ``+x``; the wake trails toward ``-x`` and this
   free-wave reconstruction is physical only *astern* of the hull.

The inner integral ``F`` (step 1) is exact; only the smooth *outer*
integrals over the propagation angle ``θ`` -- resistance (step 3) and the
wake reconstruction (step 4) -- are done by simple trapezoidal quadrature,
chosen for legibility over the crate's adaptive, oscillation-aware panels.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .bspline import BSplineSurface
from .moments import exp_moments, osc_moments

STANDARD_GRAVITY = 9.80665
SEAWATER_DENSITY = 1025.9  # kg/m³, salt water at 15 °C (ITTC)
FRESHWATER_DENSITY = 999.1


@dataclass
class Conditions:
    """Steady-ahead operating conditions in calm, deep water."""

    speed: float  # U [m/s]
    density: float = SEAWATER_DENSITY  # ρ [kg/m³]
    gravity: float = STANDARD_GRAVITY  # g [m/s²]

    @property
    def nu(self) -> float:
        """Fundamental wavenumber ν = g/U² [1/m]."""
        return self.gravity / (self.speed * self.speed)

    def froude_number(self, length: float) -> float:
        return self.speed / np.sqrt(self.gravity * length)


# np.trapz was renamed to np.trapezoid in NumPy 2.0; support both.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def _trapezoid_weights(nodes: np.ndarray) -> np.ndarray:
    """Composite-trapezoid weights for the given (uniform) sample nodes."""
    weights = np.full(len(nodes), nodes[1] - nodes[0])
    weights[0] *= 0.5
    weights[-1] *= 0.5
    return weights


class WaveField:
    """The free-wave field of one thin hull at a fixed speed.

    On construction it extracts, once, the exact local polynomial of
    ``∂f/∂x`` on every knot span (via corner Taylor data).  The inner
    integral is then *closed form*: on each span the integrand is a
    polynomial times ``e^{ikx}`` in x and a polynomial times ``e^{-κz}`` in
    z, and both reduce to the moment integrals in :mod:`pymichell.moments`.
    ``I`` and ``J`` therefore carry no quadrature error -- only the smooth
    outer θ-integrals (resistance, wake) are done numerically.
    """

    def __init__(self, surface: BSplineSurface, conditions: Conditions) -> None:
        self.surface = surface
        self.cond = conditions
        self.nu = conditions.nu
        self.x_center = surface.x_center
        self.waterline = surface.waterline
        self.p = surface.degree_x
        self.q = surface.degree_z

        # x-spans (integrated in full) and z-spans (clipped to the submerged
        # region z >= waterline; spans entirely above water are dropped).
        self.x_spans = [(start, length) for _, start, length in surface.x_spans()]
        self.z_spans = []  # (start, length, clip) below the waterline
        z_span_ids = []
        for idx, (z_knot, z_start, z_len) in enumerate(surface.z_spans()):
            if z_start + z_len <= self.waterline:
                continue  # this span is entirely above the water surface
            self.z_spans.append((z_start, z_len, max(0.0, self.waterline - z_start)))
            z_span_ids.append((idx, z_knot))

        # Local polynomial coefficients of ∂f/∂x on each (x-span, z-span):
        # coeff[sx][sz] has shape (p, q+1) with the coefficient of
        # (x - x0)^a (z - z0)^b.  From f = Σ D[a][b]/(a! b!) X^a Z^b, the
        # coefficient of X^a Z^b in ∂f/∂x is D[a+1][b] / (a! b!).
        fact = np.array([math.factorial(i) for i in range(max(self.p, self.q) + 1)], float)
        inv_fact = 1.0 / np.outer(fact[: self.p], fact[: self.q + 1])  # (p, q+1)
        self.coeff = []  # per x-span: array (n_zspan, p, q+1)
        for x_knot, _, _ in surface.x_spans():
            per_x = []
            for _, z_knot in z_span_ids:
                d = surface.corner_partials(x_knot, z_knot)  # (p+1, q+1)
                per_x.append(d[1 : self.p + 1, :] * inv_fact)
            self.coeff.append(np.array(per_x))

    # -- amplitude function ----------------------------------------------
    def inner_integral(self, lam) -> np.ndarray:
        """Michell inner integral ``F(λ) = I + i J`` for ``λ = sec θ ≥ 1``.

        Exact, span by span::

            F(λ) = Σ_{x-spans} e^{i k_x (x0 - x_c)}
                   Σ_a M_a(k_x, h_x) Σ_{z-spans} Σ_b coeff[a, b] · Z_b

        where ``M_a = ∫_0^{h_x} X^a e^{i k_x X} dX`` and ``Z_b`` is the
        z-moment of ``(z - z0)^b e^{-κ (z - z_wl)}`` over the submerged part
        of the span -- both from :mod:`pymichell.moments`.  ``k_x = ν λ``,
        ``κ = ν λ²``.  Accepts a scalar or array of ``λ`` (vectorized over
        the angles) and returns complex values of matching length.
        """
        lam = np.atleast_1d(np.asarray(lam, dtype=float))
        kx = self.nu * lam  # (L,)  oscillation rate in x
        kappa = self.nu * lam * lam  # (L,)  decay rate with depth z

        # z-moments per submerged z-span, shape (L, n_zspan, q+1).
        #   Z_b = e^{-κ (z0 - z_wl)} · ( N_b(κ, h) - N_b(κ, clip) )
        # where the clip subtracts the part of the span above the waterline.
        z_moments = []
        for z_start, z_len, clip in self.z_spans:
            n_full = exp_moments(kappa, z_len, self.q)  # (L, q+1)
            n_clip = exp_moments(kappa, clip, self.q) if clip > 0 else 0.0
            shift = np.exp(-kappa * (z_start - self.waterline))[:, None]
            z_moments.append(shift * (n_full - n_clip))
        z_stack = np.stack(z_moments, axis=1)  # (L, n_zspan, q+1)

        # x-spans: closed-form x-moments, then assemble.
        f = np.zeros(len(lam), dtype=complex)
        for (x_start, x_len), coeff in zip(self.x_spans, self.coeff):
            mx = osc_moments(kx, x_len, self.p - 1)  # (L, p)
            phase = np.exp(1j * kx * (x_start - self.x_center))  # (L,)
            # g[l, a] = Σ_{z-span} Σ_b coeff[z-span, a, b] · Z_b[l, z-span]
            g = np.einsum("lzb,zab->la", z_stack, coeff)  # (L, p)
            f += phase * np.sum(mx * g, axis=1)
        return f

    def amplitude(self, theta) -> np.ndarray:
        """Complex free-wave amplitude density ``A(θ)`` [m/rad]."""
        theta = np.atleast_1d(np.asarray(theta, dtype=float))
        sec = 1.0 / np.cos(theta)
        f = self.inner_integral(sec)
        return -(2.0 * self.nu / np.pi) * sec**3 * np.conj(f)

    def transverse_wavelength(self) -> float:
        """Wavelength of the transverse (θ = 0) waves, 2π U²/g [m]."""
        return 2.0 * np.pi / self.nu

    # -- resistance -------------------------------------------------------
    def wave_resistance(self, n_theta: int = 4000, sec_max: float = 15.0) -> float:
        """Michell wave resistance ``R_w`` [N].

        Integrates ``|F(sec θ)|² sec³θ`` over ``θ ∈ [0, π/2)`` (short waves
        beyond ``sec θ = sec_max`` are dropped -- viscosity destroys them in
        reality).  Uses the Michell form; :meth:`resistance_from_spectrum`
        computes the algebraically identical spectrum form as a check.
        """
        theta_max = np.arccos(1.0 / sec_max)
        theta = np.linspace(0.0, theta_max, n_theta)
        sec = 1.0 / np.cos(theta)
        integrand = np.abs(self.inner_integral(sec)) ** 2 * sec**3
        integral = _trapezoid(integrand, theta)
        coeff = 4.0 * self.cond.density * self.cond.gravity**2 / (np.pi * self.cond.speed**2)
        return coeff * integral

    def resistance_from_spectrum(self, n_theta: int = 4000, sec_max: float = 15.0) -> float:
        """The same ``R_w`` via the free-wave spectrum identity.

        ``R_w = ½ π ρ U² ∫ |A(θ)|² cos³θ dθ``.  Agreeing with
        :meth:`wave_resistance` demonstrates that the amplitude ``A`` is
        correctly normalised.
        """
        theta_max = np.arccos(1.0 / sec_max)
        theta = np.linspace(-theta_max, theta_max, 2 * n_theta)
        integrand = np.abs(self.amplitude(theta)) ** 2 * np.cos(theta) ** 3
        integral = _trapezoid(integrand, theta)
        return 0.5 * np.pi * self.cond.density * self.cond.speed**2 * integral

    # -- the wave amplitude grid -----------------------------------------
    def elevation_grid(
        self,
        x: np.ndarray,
        y: np.ndarray,
        n_theta: int = 1200,
        sec_max: float = 10.0,
    ) -> np.ndarray:
        """Kelvin-wake surface elevation ``ζ(x, y)`` [m] on a grid.

        ``x`` and ``y`` are 1-D coordinate arrays (metres, in the hull
        frame; the hull sits around ``x = x_center``, ``y = 0``).  Returns an
        array of shape ``(len(y), len(x))`` -- row-major over ``y``, matching
        the crate's ``WaveGrid`` -- suitable for ``imshow``/``pcolormesh``.

        Each propagation angle contributes one plane wave; we sum them by
        the trapezoidal rule over ``θ``.  ``sec_max`` caps the shortest
        (most steeply diverging) waves included.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        gx, gy = np.meshgrid(x - self.x_center, y)  # both (len(y), len(x))

        theta_max = np.arccos(1.0 / sec_max)
        theta = np.linspace(-theta_max, theta_max, n_theta)
        w_theta = _trapezoid_weights(theta)
        amp = self.amplitude(theta)
        sec = 1.0 / np.cos(theta)
        tan = np.tan(theta)

        zeta = np.zeros_like(gx, dtype=complex)
        for a, w, s, t in zip(amp, w_theta, sec, tan):
            kx = self.nu * s
            ky = self.nu * s * t
            zeta += w * a * np.exp(1j * (kx * gx + ky * gy))
        return zeta.real


def wave_field(surface: BSplineSurface, speed: float, seawater: bool = True, **kwargs) -> WaveField:
    """Convenience constructor: a :class:`WaveField` at ``speed`` [m/s]."""
    density = SEAWATER_DENSITY if seawater else FRESHWATER_DENSITY
    return WaveField(surface, Conditions(speed=speed, density=density), **kwargs)

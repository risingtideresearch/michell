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
   depend on that choice).  The crate evaluates this in closed form span by
   span; here we integrate it numerically on a grid -- the same integral,
   written plainly.

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

Everything below deliberately uses simple, dense numerical quadrature
(trapezoidal rules over uniform grids).  It is not how you would evaluate
this quickly -- it is how you would explain it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .bspline import BSplineSurface

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

    On construction it samples ``∂f/∂x`` once on a dense ``(x, z)`` grid;
    every amplitude and elevation query is then a weighted sum over that
    cached slope field.  Raise ``n_x``/``n_z`` for accuracy on wiggly hulls
    or high speeds (the ``x`` grid must resolve the wave phase ``ν λ x``).
    """

    def __init__(
        self,
        surface: BSplineSurface,
        conditions: Conditions,
        n_x: int = 400,
        n_z: int = 80,
    ) -> None:
        self.surface = surface
        self.cond = conditions
        self.nu = conditions.nu
        self.x_center = surface.x_center
        self.waterline = surface.waterline

        # Cache the slope field ∂f/∂x on a uniform grid over the *submerged*
        # region z ∈ [waterline, draft], with the trapezoid weights that
        # integrate over it.  For a plain wetted-surface hull waterline = 0.
        x0, x1 = surface.x_domain
        _, z_keel = surface.z_domain
        self.x_nodes = np.linspace(x0, x1, n_x)
        self.z_nodes = np.linspace(surface.waterline, z_keel, n_z)
        self.wx = _trapezoid_weights(self.x_nodes)
        self.wz = _trapezoid_weights(self.z_nodes)
        self.fx = surface.evaluate(self.x_nodes, self.z_nodes, dx=1, dz=0)  # (n_x, n_z)

    # -- amplitude function ----------------------------------------------
    def inner_integral(self, lam) -> np.ndarray:
        """Michell inner integral ``F(λ) = I + i J`` for ``λ = sec θ ≥ 1``.

        Accepts a scalar or array of ``λ`` and returns complex values of the
        same shape.  This is the double integral of the docstring, done as
        two nested weighted sums over the cached slope grid.
        """
        lam = np.atleast_1d(np.asarray(lam, dtype=float))
        kx = self.nu * lam  # oscillation rate in x
        kappa = self.nu * lam * lam  # decay rate with depth z

        # ∫ (∂f/∂x) e^{-κ (z - z_wl)} dz  for every x station and every λ,
        # where z - z_wl is the depth below the waterline.
        below = self.z_nodes[None, :] - self.waterline
        depth = np.exp(-kappa[:, None] * below) * self.wz[None, :]  # (L, n_z)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            g = depth @ self.fx.T  # (L, n_x)  (errstate: see BSplineSurface.evaluate)

        # ∫ g(x) e^{i k_x (x - x_c)} dx.
        phase = np.exp(1j * kx[:, None] * (self.x_nodes[None, :] - self.x_center))  # (L, n_x)
        return np.sum(self.wx[None, :] * g * phase, axis=1)

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

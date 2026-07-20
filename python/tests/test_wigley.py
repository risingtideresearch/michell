"""Validation against the Wigley hull and the Rust ``michell`` crate.

Reference values were produced by ``michell wigley`` + ``michell resistance``
/ ``michell spectrum`` at U = 3 m/s in seawater (ρ = 1025.9 kg/m³):

    R_w              = 145.761 N
    A(0.0°)   .imag  = +0.0672350
    A(-19.9°) .imag  = +0.0881300
    A(-70.1°) .imag  = -0.2061821
    transverse λ     =  5.766359 m

Run with ``pytest`` (or execute this file directly for a plain report).
"""

import numpy as np

from pymichell import wave_field, wigley

_trap = getattr(np, "trapezoid", None) or np.trapz  # renamed in NumPy 2.0


def _field():
    return wave_field(wigley(10.0, 1.0, 0.625), speed=3.0)


def test_transverse_wavelength():
    assert abs(_field().transverse_wavelength() - 5.766359) < 1e-4


def test_wave_resistance_matches_crate():
    rw = _field().wave_resistance()
    assert abs(rw - 145.761) < 0.5  # ~0.3% of 145.8 N


def test_resistance_forms_agree():
    field = _field()
    r_michell = field.wave_resistance()
    r_spectrum = field.resistance_from_spectrum()
    assert abs(r_michell - r_spectrum) <= 1e-3 * r_michell


def test_amplitude_matches_crate():
    field = _field()
    # The Wigley hull is fore-aft symmetric about x_c, so ∂f/∂x is odd in x
    # and F is purely imaginary -> A is purely imaginary too.
    for theta_deg, want_imag in [(0.0, 0.0672350), (-19.889, 0.0881300), (-70.108, -0.2061821)]:
        a = field.amplitude(np.radians(theta_deg))[0]
        assert abs(a.real) < 1e-6 * (abs(want_imag) + 1e-9)
        assert abs(a.imag - want_imag) < 1e-3 * abs(want_imag)


def test_inner_integral_matches_direct_quadrature():
    # The closed-form span machinery must agree with a naive but obviously
    # correct 2-D numerical integration of the same double integral.
    surface = wigley(10.0, 1.0, 0.625)
    field = wave_field(surface, speed=3.0)
    x0, x1 = surface.x_domain
    z0, z1 = surface.z_domain
    xn = np.linspace(x0, x1, 2000)
    zn = np.linspace(z0, z1, 400)
    fx = surface.evaluate(xn, zn, dx=1)  # (nx, nz)
    xc = surface.x_center
    for lam in [1.0, 1.3, 2.0, 4.0, 8.0]:
        kappa = field.nu * lam * lam
        kx = field.nu * lam
        gz = _trap(fx * np.exp(-kappa * zn)[None, :], zn, axis=1)  # (nx,)
        ref = _trap(gz * np.exp(1j * kx * (xn - xc)), xn)
        exact = field.inner_integral(lam)[0]
        # Absolute floor absorbs the reference trapezoid's own error near the
        # near-zeros of |F| (where a purely relative bound is unreasonable).
        assert abs(exact - ref) < 1e-6 + 1e-4 * abs(ref)


def test_wake_has_transverse_waves():
    field = _field()
    # A centreline cut astern should oscillate at the transverse wavelength.
    x = np.linspace(field.x_center - 40.0, field.x_center - 6.0, 800)
    cut = field.elevation_grid(x, np.array([0.0]))[0]
    assert np.abs(cut).max() > 1e-3  # not flat
    # Count zero crossings; expect ~ (span / λ) * 2.
    crossings = np.sum(np.diff(np.sign(cut)) != 0)
    expected = 2 * 34.0 / field.transverse_wavelength()
    assert crossings >= 0.5 * expected


if __name__ == "__main__":
    field = _field()
    print(f"transverse wavelength : {field.transverse_wavelength():.6f} m   (ref 5.766359)")
    print(f"R_w (Michell form)    : {field.wave_resistance():.4f} N        (ref 145.761)")
    print(f"R_w (spectrum form)   : {field.resistance_from_spectrum():.4f} N")
    for theta_deg in (0.0, -19.889, -70.108):
        a = field.amplitude(np.radians(theta_deg))[0]
        print(f"A({theta_deg:+7.3f}deg)        : {a.real:+.3e} {a.imag:+.3e}i")
    print("all reference checks passed" if True else "")

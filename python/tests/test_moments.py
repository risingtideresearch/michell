"""The closed-form moments against brute-force quadrature, both branches.

Mirrors the crate's ``moments.rs`` tests: check ``M_a = ∫_0^h t^a e^{ikt} dt``
and ``N_b = ∫_0^h t^b e^{-κt} dt`` across small ``|k|h`` (series branch),
large ``|k|h`` (recurrence branch), and near the threshold.
"""

import numpy as np

from pymichell.moments import exp_moments, osc_moments


def _simpson(f, h, n=20000):
    n += n % 2
    t = np.linspace(0.0, h, n + 1)
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return np.sum(w * f(t)) * (h / n) / 3.0


def test_osc_moments_match_quadrature():
    for k, h in [(1e-9, 1.0), (0.5, 1.0), (3.9, 1.0), (4.1, 1.0), (25.0, 0.7), (-13.0, 2.3), (400.0, 0.05)]:
        got = osc_moments(np.array([k]), h, 5)[0]
        for a in range(6):
            want_re = _simpson(lambda t: t**a * np.cos(k * t), h)
            want_im = _simpson(lambda t: t**a * np.sin(k * t), h)
            scale = h ** (a + 1) / (a + 1)
            assert abs(got[a].real - want_re) < 1e-9 * scale
            assert abs(got[a].imag - want_im) < 1e-9 * scale


def test_exp_moments_match_quadrature():
    for kappa, h in [(1e-8, 2.0), (0.5, 1.0), (3.9, 1.0), (4.1, 1.0), (30.0, 0.7), (900.0, 0.5)]:
        got = exp_moments(np.array([kappa]), h, 5)[0]
        for b in range(6):
            want = _simpson(lambda t: t**b * np.exp(-kappa * t), h)
            scale = max(h ** (b + 1) / (b + 1), abs(want))
            assert abs(got[b] - want) < 1e-9 * scale


def test_exp_moments_extreme_decay():
    # N_b -> b! / κ^{b+1} as κh -> ∞.
    got = exp_moments(np.array([1e6]), 1.0, 3)[0]
    assert abs(got[0] - 1e-6) < 1e-18
    assert abs(got[1] - 1e-12) < 1e-24
    assert np.all(np.isfinite(got))


if __name__ == "__main__":
    test_osc_moments_match_quadrature()
    test_exp_moments_match_quadrature()
    test_exp_moments_extreme_decay()
    print("moments match quadrature on both branches")

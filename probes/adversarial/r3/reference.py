#!/usr/bin/env python3
"""Independent high-precision endpoint/pair reference for stress geometries.

This program does not call the michell crate for its reference values.  It
evaluates B-spline bases recursively with mpmath, reconstructs each polynomial
span by interpolation, performs integration by parts independently, and sums
all endpoint pairs.  Waterline kernels use high-precision adaptive quadrature;
submerged kernels use independently generated, phase-resolved real-axis
Gauss-Legendre panels.
"""

from __future__ import annotations

import argparse
import cmath
import math
import subprocess
from dataclasses import dataclass
from functools import lru_cache

import mpmath as mp
import numpy as np


@dataclass(frozen=True)
class Case:
    name: str
    px: int
    pz: int
    kx: tuple[float, ...]
    kz: tuple[float, ...]
    control: tuple[float, ...]
    froude: float


def tensor(gx, gz, scale=0.5):
    return tuple(scale * x * z for x in gx for z in gz)


def wigley(name, draft, froude):
    return Case(
        name,
        2,
        2,
        (-5.0, -5.0, -5.0, 5.0, 5.0, 5.0),
        (0.0, 0.0, 0.0, draft, draft, draft),
        tensor((0.0, 2.0, 0.0), (1.0, 1.0, 0.0)),
        froude,
    )


def cases():
    high = Case(
        "degree5x4_deep_fn005",
        5,
        4,
        (-5.0,) * 6 + (5.0,) * 6,
        (0.0,) * 5 + (2.0,) * 5,
        tensor((0.0, 0.45, 1.35, 1.55, 0.55, 0.0), (1.0, 0.95, 0.72, 0.33, 0.0)),
        0.05,
    )
    very_high = Case(
        "degree8x8_deep_fn005",
        8,
        8,
        (-5.0,) * 9 + (5.0,) * 9,
        (0.0,) * 9 + (2.0,) * 9,
        tensor(
            (0.0, 0.20, 0.75, 1.35, 1.60, 1.42, 0.88, 0.26, 0.0),
            (1.0, 0.99, 0.94, 0.82, 0.65, 0.44, 0.23, 0.08, 0.0),
        ),
        0.05,
    )
    extreme_px = 16
    extreme_pz = 12
    extreme = Case(
        "degree16x12_deep_fn005",
        extreme_px,
        extreme_pz,
        (-5.0,) * (extreme_px + 1) + (5.0,) * (extreme_px + 1),
        (0.0,) * (extreme_pz + 1) + (2.0,) * (extreme_pz + 1),
        tensor(
            tuple(math.sin(math.pi * index / extreme_px) for index in range(extreme_px + 1)),
            tuple(1.0 - index / extreme_pz for index in range(extreme_pz + 1)),
        ),
        0.05,
    )
    more_extreme_px = 24
    more_extreme_pz = 16
    more_extreme = Case(
        "degree24x16_deep_fn005",
        more_extreme_px,
        more_extreme_pz,
        (-5.0,) * (more_extreme_px + 1) + (5.0,) * (more_extreme_px + 1),
        (0.0,) * (more_extreme_pz + 1) + (2.0,) * (more_extreme_pz + 1),
        tensor(
            tuple(math.sin(math.pi * index / more_extreme_px) for index in range(more_extreme_px + 1)),
            tuple(1.0 - index / more_extreme_pz for index in range(more_extreme_pz + 1)),
        ),
        0.05,
    )
    chine = Case(
        "degree2_chine_multispan_fn003",
        2,
        2,
        (0.0, 0.0, 0.0, 0.8, 0.8, 2.0, 2.0, 2.0),
        (0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0),
        (
            0.0,
            0.0,
            0.0,
            0.5,
            0.4,
            0.1,
            0.7,
            0.5,
            0.1,
            0.6,
            0.45,
            0.05,
            0.3,
            0.2,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ),
        0.03,
    )
    close_control = tensor(
        (0.0, 0.55, 1.0, 1.1, 1.30, 1.20, 1.1, 0.80, 0.40, 0.0),
        (1.0, 0.9, 0.45, 0.0),
    )
    return [
        wigley("wigley_fn008", 0.625, 0.08),
        wigley("wigley_fn005", 0.625, 0.05),
        wigley("wigley_gate_exact_fn020", 0.625, 0.20),
        wigley("wigley_gate_reject_fn02001", 0.625, 0.2001),
        wigley("wigley_deep_fn008", 2.5, 0.08),
        wigley("wigley_shallow_fn005", 0.05, 0.05),
        high,
        very_high,
        extreme,
        more_extreme,
        chine,
        *[
            Case(
                name,
                3,
                3,
                (-5.0,) * 4 + (-delta,) * 3 + (delta,) * 3 + (5.0,) * 4,
                (0.0,) * 4 + (2.0,) * 4,
                close_control,
                froude,
            )
            for name, delta, froude in (
                ("close_chines_gate_pass", 0.01, 0.0088),
                ("close_chines_gate_reject", 0.01, 0.0091),
                ("close_chines_delta_1e-4", 1e-4, math.sqrt(1e-4 / 130.0)),
                ("close_chines_delta_1e-7", 1e-7, math.sqrt(1e-7 / 130.0)),
                ("close_chines_delta_1e-10", 1e-10, math.sqrt(1e-10 / 130.0)),
            )
        ],
    ]


@lru_cache(maxsize=None)
def basis(i, p, u, knots):
    if p == 0:
        return mp.mpf(1) if knots[i] <= u < knots[i + 1] else mp.mpf(0)
    left_den = knots[i + p] - knots[i]
    right_den = knots[i + p + 1] - knots[i + 1]
    left = (u - knots[i]) * basis(i, p - 1, u, knots) / left_den if left_den else 0
    right = (
        (knots[i + p + 1] - u) * basis(i + 1, p - 1, u, knots) / right_den
        if right_den
        else 0
    )
    return left + right


@lru_cache(maxsize=None)
def basis_derivative(i, p, u, knots):
    if p == 0:
        return mp.mpf(0)
    left_den = knots[i + p] - knots[i]
    right_den = knots[i + p + 1] - knots[i + 1]
    left = p * basis(i, p - 1, u, knots) / left_den if left_den else 0
    right = p * basis(i + 1, p - 1, u, knots) / right_den if right_den else 0
    return left - right


def surface_fx(case, x, z):
    nx = len(case.kx) - case.px - 1
    nz = len(case.kz) - case.pz - 1
    total = mp.mpf(0)
    kx = tuple(mp.mpf(v) for v in case.kx)
    kz = tuple(mp.mpf(v) for v in case.kz)
    for i in range(nx):
        dx = basis_derivative(i, case.px, x, kx)
        if not dx:
            continue
        for j in range(nz):
            bz = basis(j, case.pz, z, kz)
            total += mp.mpf(case.control[i * nz + j]) * dx * bz
    return total


def nonzero_spans(knots, degree):
    nctrl = len(knots) - degree - 1
    return [(mp.mpf(knots[i]), mp.mpf(knots[i + 1])) for i in range(degree, nctrl) if knots[i + 1] > knots[i]]


def falling(n, r):
    return mp.factorial(n) / mp.factorial(n - r)


def reconstruct_terms(case):
    length = mp.mpf(case.kx[-1]) - mp.mpf(case.kx[0])
    nu = 1 / (mp.mpf(case.froude) ** 2 * length)
    combined = {}
    for xa, xb in nonzero_spans(case.kx, case.px):
        for za, zb in nonzero_spans(case.kz, case.pz):
            nxp = case.px
            nzp = case.pz + 1
            xs = [xa + (xb - xa) * (mp.mpf(i) + mp.mpf("0.5")) / nxp for i in range(nxp)]
            zs = [za + (zb - za) * (mp.mpf(j) + mp.mpf("0.5")) / nzp for j in range(nzp)]
            vx = mp.matrix([[(x - xa) ** a for a in range(nxp)] for x in xs])
            vz = mp.matrix([[(z - za) ** b for b in range(nzp)] for z in zs])
            samples = mp.matrix([[surface_fx(case, x, z) for z in zs] for x in xs])
            coeff = (vx ** -1) * samples * (vz ** -1).T
            dx = xb - xa
            dz = zb - za
            for a in range(nxp):
                for b in range(nzp):
                    polynomial = coeff[a, b]
                    if abs(polynomial) < mp.mpf("1e-50"):
                        continue
                    for r in range(a + 1):
                        derivative_x = falling(a, r)
                        x_ends = [(xb, (-1) ** r * derivative_x * dx ** (a - r))]
                        if r == a:
                            x_ends.append((xa, -((-1) ** r) * derivative_x))
                        for u in range(b + 1):
                            derivative_z = falling(b, u)
                            z_ends = [(zb, -derivative_z * dz ** (b - u))]
                            if u == b:
                                z_ends.append((za, derivative_z))
                            power = r + 2 * u + 3
                            scale = polynomial / nu ** (r + u + 2) / (1j) ** (r + 1)
                            for xe, xf in x_ends:
                                for ze, zf in z_ends:
                                    key = (xe, ze, power)
                                    combined[key] = combined.get(key, 0j) + scale * xf * zf
    return nu, [(x, z, n, c) for (x, z, n), c in combined.items() if abs(c) > mp.mpf("1e-42")]


def zero_kernel(s):
    return mp.sqrt(mp.pi) * mp.gamma(mp.mpf(s) / 2) / (2 * mp.gamma(mp.mpf(s + 1) / 2))


@lru_cache(maxsize=None)
def water_kernel(s, omega_text):
    omega = mp.mpf(omega_text)
    if not omega:
        return zero_kernel(s)
    if omega < 0:
        return mp.conj(water_kernel(s, mp.nstr(-omega, 50)))
    phase = 2 * mp.e ** (1j * (omega + mp.pi / 4)) / mp.sqrt(omega)

    def integrand(y):
        return mp.e ** (-y * y) / (
            (1 + 1j * y * y / omega) ** s * mp.sqrt(2 + 1j * y * y / omega)
        )

    return phase * mp.quad(integrand, [0, 1, 2, 4, 8, 12])


GL_X, GL_W = np.polynomial.legendre.leggauss(32)


@lru_cache(maxsize=None)
def submerged_kernel(s, omega_float, a_float, phase_divisor):
    omega = float(omega_float)
    a = float(a_float)
    if a > 60.0:
        return 0j
    # Cut when decay relative to lambda=1 is exp(-92), then resolve each
    # oscillatory phase interval by at least phase_divisor panels per pi.
    lambda_max = math.sqrt(1.0 + 92.0 / a)
    u_max = lambda_max - 1.0
    du = min(0.04, math.pi / (phase_divisor * max(abs(omega), 1.0)))
    panels = max(1, math.ceil(u_max / du))
    total = 0j
    for panel in range(panels):
        u0 = panel * u_max / panels
        u1 = (panel + 1) * u_max / panels
        t0, t1 = math.sqrt(u0), math.sqrt(u1)
        half, mid = 0.5 * (t1 - t0), 0.5 * (t1 + t0)
        t = mid + half * GL_X
        lam = 1.0 + t * t
        values = (
            2.0
            * np.exp(1j * omega * lam - a * lam * lam)
            / np.power(lam, s)
            / np.sqrt(2.0 + t * t)
        )
        total += half * np.dot(GL_W, values)
    return complex(total)


def reference(case, phase_divisor):
    nu, terms = reconstruct_terms(case)
    total = mp.mpc(0)
    for i, (xi, zi, ni, ci) in enumerate(terms):
        for j in range(i, len(terms)):
            xj, zj, nj, cj = terms[j]
            s = ni + nj - 2
            omega = nu * (xi - xj)
            a = nu * (zi + zj)
            if not a:
                kernel = water_kernel(s, mp.nstr(omega, 50))
            else:
                kernel = submerged_kernel(s, float(omega), float(a), phase_divisor)
            multiplicity = 1 if i == j else 2
            total += multiplicity * mp.re(ci * mp.conj(cj) * kernel)
    length = mp.mpf(case.kx[-1]) - mp.mpf(case.kx[0])
    speed_sq = mp.mpf(case.froude) ** 2 * mp.mpf("9.80665") * length
    physical = 4 * mp.mpf("999.1") * mp.mpf("9.80665") ** 2 / (mp.pi * speed_sq)
    return physical * total, len(terms), max(n for _, _, n, _ in terms)


def run_production(binary):
    lines = subprocess.check_output([binary], text=True).splitlines()
    return {line.split("|", 1)[0]: line.split("|") for line in lines if "|" in line}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--phase-divisor", type=int, default=8)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--case")
    args = parser.parse_args()
    mp.mp.dps = args.dps
    production = run_production(args.binary)
    print("name|reference|full_terms|max_power|reduced_status|reduced_rel_error|estimate|covers|hybrid_method|hybrid_rel_error")
    for case in cases():
        if args.case and case.name != args.case:
            continue
        ref, full_terms, max_power = reference(case, args.phase_divisor)
        fields = production[case.name]
        if fields[2] == "OK":
            reduced = mp.mpf(fields[3])
            estimate = mp.mpf(fields[4])
            reduced_error = abs(reduced - ref) / abs(ref)
            method = fields[11]
            hybrid = mp.mpf(fields[13])
        else:
            reduced_error = mp.nan
            estimate = mp.nan
            method = fields[4]
            hybrid = mp.mpf(fields[6])
        hybrid_error = abs(hybrid - ref) / abs(ref)
        covers = "NA" if mp.isnan(reduced_error) else str(bool(reduced_error <= estimate))
        print(
            f"{case.name}|{mp.nstr(ref, 17)}|{full_terms}|{max_power}|{fields[2]}|"
            f"{mp.nstr(reduced_error, 8)}|{mp.nstr(estimate, 8)}|{covers}|{method}|{mp.nstr(hybrid_error, 8)}"
        )


if __name__ == "__main__":
    main()

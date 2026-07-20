"""Compute and (optionally) plot the Kelvin wake of a hull.

Usage::

    python examples/demo.py                 # Wigley reference hull, U = 3 m/s
    python examples/demo.py boat.hull 8     # a .hull file at 8 m/s

Prints the wave resistance and, if matplotlib is installed, saves a heatmap
of the wave-amplitude grid ζ(x, y) to ``wake.png``.
"""

import sys

import numpy as np

from pymichell import load_hull, wave_field, wigley


def main() -> None:
    args = sys.argv[1:]
    if args and not args[0].replace(".", "").isdigit():
        surface = load_hull(args[0])
        speed = float(args[1]) if len(args) > 1 else 3.0
        label = args[0]
    else:
        surface = wigley(10.0, 1.0, 0.625)
        speed = float(args[0]) if args else 3.0
        label = "Wigley L=10 B=1 T=0.625"

    field = wave_field(surface, speed=speed)
    length = surface.length

    print(f"hull                  : {label}")
    print(f"speed                 : {speed:.3f} m/s   (Fn = {field.cond.froude_number(length):.3f})")
    print(f"transverse wavelength : {field.transverse_wavelength():.3f} m")
    print(f"wave resistance R_w   : {field.wave_resistance():.3f} N")
    print(f"  (spectrum identity) : {field.resistance_from_spectrum():.3f} N")

    # The wave-amplitude grid: a box trailing astern of the hull (toward -x).
    xc = surface.x_center
    x = np.linspace(xc - 5.0 * length, xc + 1.0 * length, 320)
    y = np.linspace(-2.0 * length, 2.0 * length, 240)
    zeta = field.elevation_grid(x, y)
    print(f"wave grid             : {zeta.shape[1]} x {zeta.shape[0]}, peak |ζ| = {np.abs(zeta).max():.4f} m")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("(install matplotlib to render wake.png)")
        return

    amp = np.abs(zeta).max()
    plt.figure(figsize=(9, 6))
    plt.pcolormesh(x, y, zeta, cmap="RdBu_r", vmin=-amp, vmax=amp, shading="auto")
    plt.colorbar(label="surface elevation ζ [m]")
    plt.axis("equal")
    plt.xlabel("x [m]  (ship advances toward +x)")
    plt.ylabel("y [m]")
    plt.title(f"Kelvin wake -- {label}, U = {speed:g} m/s")
    plt.tight_layout()
    plt.savefig("wake.png", dpi=130)
    print("wrote wake.png")


if __name__ == "__main__":
    main()

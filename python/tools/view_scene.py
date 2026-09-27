#!/usr/bin/env python3
"""Show a `michell field` scene in polyscope.

    michell field e12.igs --waterline -0.95 --froude 0.3 -o scene.json
    python python/tools/view_scene.py scene.json [--wave-scale 20] [--screenshot out.png]

Every object in the scene becomes a polyscope structure (meshes as surface
meshes, station curves as curve networks) with each of its quantities as a
vertex scalar you can switch between in the UI. The free surface's height can
be exaggerated with --wave-scale (its colour always shows the true ζ).

Needs Python >= 3.10 with `polyscope` and `numpy` (e.g. `uv run --with
polyscope --with numpy python python/tools/view_scene.py scene.json`).
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import polyscope as ps


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scene", help="a scene written by `michell field`")
    ap.add_argument("--wave-scale", type=float, default=1.0,
                    help="vertical exaggeration of the free surface (colour stays true)")
    ap.add_argument("--screenshot", help="render once to this PNG and exit")
    ap.add_argument("--hide", action="append", default=[],
                    help="hide objects whose name contains this, or is exactly it after a "
                         "leading '=' (repeatable)")
    ap.add_argument("--camera", nargs=6, type=float, metavar=("EX", "EY", "EZ", "TX", "TY", "TZ"),
                    help="screenshot eye and target points")
    args = ap.parse_args()

    with open(args.scene) as f:
        scene = json.load(f)
    if scene.get("michell") != "scene":
        raise SystemExit(f"{args.scene}: not a michell scene")
    meta = scene.get("meta", {})

    ps.set_program_name("michell scene")
    ps.set_up_dir("z_up")
    ps.set_front_dir("neg_y_front")
    ps.set_ground_plane_mode("none")
    ps.init()

    for obj in scene["objects"]:
        v = np.asarray(obj["vertices"], dtype=float)
        name = obj["name"]
        if obj["kind"] == "mesh":
            if "free surface" in name and args.wave_scale != 1.0:
                v = v.copy()
                v[:, 2] *= args.wave_scale
            s = ps.register_surface_mesh(name, v, np.asarray(obj["faces"], dtype=int))
            if name.endswith(" hull") or "CAD" in name:
                s.set_transparency(0.35)
                s.set_color((0.85, 0.78, 0.66))
            if "free surface" in name or name.endswith(" pressure"):
                # Unshaded, so the colour reads as the value.
                s.set_material("flat")
            if "closure" in name:
                s.set_color((0.91, 0.35, 0.05))
        elif obj["kind"] == "curves":
            s = ps.register_curve_network(name, v, np.asarray(obj["edges"], dtype=int), radius=0.0006)
        else:
            continue
        for i, q in enumerate(obj.get("quantities", [])):
            vals = np.asarray([np.nan if x is None else x for x in q["values"]], dtype=float)
            symmetric = np.nanmin(vals) < 0 < np.nanmax(vals)
            field = "free surface" in name or name.endswith(" pressure")
            kw = dict(
                defined_on="nodes" if obj["kind"] == "curves" else "vertices",
                cmap="coolwarm" if symmetric else "viridis",
                enabled=(i == 0 and obj["kind"] != "mesh") or field,
            )
            if symmetric:
                # A robust scale: a stagnation spike should not wash out the rest.
                peak = float(np.nanpercentile(np.abs(vals), 99.0))
                kw["vminmax"] = (-peak, peak)
            s.add_scalar_quantity(q["name"], vals, **kw)
        if any(h[1:] == name if h.startswith("=") else h in name for h in args.hide):
            s.set_enabled(False)
        if obj.get("note"):
            print(f"{name}: {obj['note']}")

    print(f"U = {meta.get('speed', float('nan')):.3f} m/s, Fn {meta.get('froude', float('nan')):.3f}, "
          f"transverse wavelength {meta.get('transverse_wavelength', float('nan')):.3f} m, closure {meta.get('closure')}")
    if args.screenshot:
        if args.camera:
            ps.look_at(tuple(args.camera[:3]), tuple(args.camera[3:]))
        else:
            ps.look_at((0.0, -30.0, 18.0), (-6.0, 0.0, 0.0))
        ps.screenshot(args.screenshot, transparent_bg=False)
    else:
        ps.show()


if __name__ == "__main__":
    main()

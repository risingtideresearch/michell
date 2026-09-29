// The 3D view: a hull as uploaded and as cut into sections, or a case's
// flow — the pressure on each hull and the free surface around them.
//
// `createViewer(host)` fills `host` (positioned, sized by the page) with the
// canvas and its overlays: the view buttons, a readout under the pointer,
// the colour legends. Hull frame: x along, y transverse, z DOWN from the
// waterline. Scene frame: X = x, Y = -z (up), Z = y.

import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { b64 } from "/static/common.js";

// Diverging colour for t in [−1, 1]: blue (negative), white, red (positive).
const RAMP = [[0.129, 0.4, 0.675], [0.573, 0.773, 0.871], [0.969, 0.969, 0.969], [0.957, 0.647, 0.51], [0.698, 0.094, 0.169]];
function diverging(t) {
  const u = (Math.max(-1, Math.min(1, t)) + 1) * 2, i = Math.min(3, Math.floor(u)), f = u - i;
  return RAMP[i].map((c, k) => c + f * (RAMP[i + 1][k] - c));
}
// A robust symmetric range: the 99th percentile of |v|.
function range99(vals) {
  const a = Float64Array.from(vals, Math.abs).sort();
  return a.length ? Math.max(a[Math.floor(0.99 * (a.length - 1))], 1e-12) : 1;
}
// φ(s) = 1 − 3s² + 2s³: 1 at the transom, 0 at the hollow's end, flat at both.
export const phi = (s) => 1 - 3 * s * s + 2 * s * s * s;

// A tessellation from the server (file frame, z up) as scene positions.
function meshGeometry(m) {
  const v = new Float32Array(b64(m.vertices)), t = new Uint32Array(b64(m.triangles));
  const pos = new Float32Array(v.length);
  for (let i = 0; i < v.length; i += 3) { pos[i] = v[i]; pos[i + 1] = v[i + 2]; pos[i + 2] = v[i + 1]; }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  geo.setIndex(new THREE.BufferAttribute(t, 1));
  geo.computeVertexNormals();
  geo.computeBoundingBox();
  return geo;
}

export function createViewer(host) {
  host.classList.add("viewer");
  host.insertAdjacentHTML("beforeend", `
    <div class="empty">Nothing to show yet</div>
    <div class="toolbar"><div class="group views">
      <button data-view="iso" class="on">3D</button>
      <button data-view="profile">Profile</button>
      <button data-view="plan">Plan</button>
      <button data-view="body">Body</button>
      <button data-view="wake" hidden>Wake</button>
    </div></div>
    <div class="readout"></div>
    <div class="legend cp"><div>Hull pressure C<sub>p</sub></div><div class="bar"></div><div class="ends"><span class="lo"></span><span>0</span><span class="hi"></span></div></div>
    <div class="legend zeta"><div>Free-surface elevation ζ</div><div class="bar"></div><div class="ends"><span class="lo"></span><span>0</span><span class="hi"></span></div></div>`);
  const q = (sel) => host.querySelector(sel);
  const views = q(".views"), readout = q(".readout"), empty = q(".empty");

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(devicePixelRatio);
  renderer.localClippingEnabled = true;
  // Keeps what is above the water (scene y ≥ 0).
  const aboveWater = new THREE.Plane(new THREE.Vector3(0, 1, 0), 0);
  host.prepend(renderer.domElement);
  const scene = new THREE.Scene();
  const camera = new THREE.OrthographicCamera(-1, 1, 1, -1, -1000, 1000);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.addEventListener("change", render);
  scene.add(new THREE.HemisphereLight(0xffffff, 0x445566, 1.6));
  const sun = new THREE.DirectionalLight(0xffffff, 1.6); sun.position.set(0.4, 1, 0.7); scene.add(sun);
  const fillLight = new THREE.DirectionalLight(0xffffff, 0.6); fillLight.position.set(-0.5, -0.3, -0.8); scene.add(fillLight);

  const model = new THREE.Group();
  scene.add(model);
  const NAMES = ["surface", "stations", "keel", "rays", "closure", "water", "geom", "hull", "pressure", "waves", "seaway"];
  const layers = Object.fromEntries(NAMES.map((n) => [n, new THREE.Group()]));
  for (const g of Object.values(layers)) model.add(g);
  // What is shown: the sections (a hull), or the flow (a case). The
  // sections view is in the file's frame (z up): the whole hull above the
  // waterline, the cut below it — the cut's layers, relative to the
  // waterline they were cut at, raised by it. The flow view is relative to
  // the solved waterline.
  const SHOWN = { cut: ["surface", "stations", "closure", "geom", "water"], flow: ["hull", "pressure", "waves"], sea: ["hull", "seaway"] };
  const CUT_FRAME = ["surface", "stations", "keel", "rays", "closure"];
  let mode = "cut", wlCut = 0, geomBounds = null;
  function show(which) {
    mode = which;
    for (const [name, g] of Object.entries(layers)) {
      g.visible = SHOWN[which].includes(name);
      if (CUT_FRAME.includes(name)) g.position.y = which === "cut" ? wlCut : 0;
    }
  }
  show("cut");
  // Keeps what is above the chosen waterline (the geometry's topsides).
  const aboveWl = new THREE.Plane(new THREE.Vector3(0, 1, 0), 0);
  let bounds = null, hullBounds = null, wakeBounds = null, view = "iso", surfaces = [], probes = [];

  const clear = (g) => { for (const c of [...g.children]) { g.remove(c); c.traverse((o) => { o.geometry?.dispose(); o.material?.dispose(); }); } };
  const P = (h, x, z, y) => [x + h.placement.x, -z, h.placement.y + y];
  const seg = (arr, color, opacity = 1, depthTest = true) => {
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.Float32BufferAttribute(arr, 3));
    return new THREE.LineSegments(g, new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity, depthTest }));
  };

  // The whole of an upload, topsides and all, with the water at `wl`.
  function geometry(meshes, wl) {
    clear(layers.geom); clear(layers.water);
    const gb = new THREE.Box3();
    for (const m of meshes) {
      const geo = meshGeometry(m);
      gb.union(geo.boundingBox);
      layers.geom.add(new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ color: 0xe9e5db, roughness: 0.7, side: THREE.DoubleSide, clippingPlanes: [aboveWl] })));
    }
    geomBounds = gb;
    const size = gb.getSize(new THREE.Vector3()), mid = gb.getCenter(new THREE.Vector3());
    const wgeo = new THREE.PlaneGeometry(size.x * 1.25, Math.max(size.z * 1.6, size.x * 0.25));
    wgeo.rotateX(-Math.PI / 2);
    const water = new THREE.Mesh(wgeo, new THREE.MeshBasicMaterial({ color: 0x3d8bd9, transparent: true, opacity: 0.22, side: THREE.DoubleSide, depthWrite: false }));
    water.position.set(mid.x, 0, mid.z);
    layers.water.add(water);
    setWater(wl);
    empty.style.display = "none";
    show("cut");
    setView(view);
  }
  function setWater(wl) {
    aboveWl.constant = -wl;
    layers.water.position.y = wl;
    render();
  }

  // A structured (nx × nz) mesh on both sides of the centreplane.
  function sheet(h, nx, nz, pt, color, opacity) {
    const group = new THREE.Group();
    for (const side of [1, -1]) {
      const pos = new Float32Array(nx * nz * 3);
      for (let i = 0; i < nx; i++) for (let j = 0; j < nz; j++) {
        const [x, y, z] = pt(i, j);
        pos.set(P(h, x, z, side * y), 3 * (i * nz + j));
      }
      const index = [];
      for (let i = 0; i + 1 < nx; i++) for (let j = 0; j + 1 < nz; j++) {
        const a = i * nz + j, b = a + nz, c = b + 1, d = a + 1;
        if (side > 0) index.push(a, b, c, a, c, d); else index.push(a, c, b, a, d, c);
      }
      const geo = new THREE.BufferGeometry();
      geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
      geo.setIndex(index);
      geo.computeVertexNormals();
      const mat = new THREE.MeshStandardMaterial({ color, roughness: 0.65, side: THREE.DoubleSide,
        transparent: opacity < 1, opacity, depthWrite: opacity >= 1 });
      const mesh = new THREE.Mesh(geo, mat);
      mesh.userData = { h };
      group.add(mesh);
    }
    return group;
  }

  // A cut's hulls, drawn as the physics uses them, at the waterline `wl`
  // (file frame) they were cut at.
  let cutHulls = [];
  function load(hs, wl = 0) {
    wlCut = wl;
    cutHulls = hs;
    for (const [name, g] of Object.entries(layers)) if (name !== "geom" && name !== "water") clear(g);
    surfaces = [];
    const box = new THREE.Box3();
    for (const h of hs) {
      // The see-through surface between stations: orientation only.
      const m = h.mesh;
      const s = sheet(h, m.nx, m.nz, (i, j) => { const k = i * m.nz + j; return [m.x[k], m.y[k], m.z[k]]; }, 0xd8b0c8, 0.3);
      layers.surface.add(s);
      surfaces.push(...s.children);
      // Stations as the physics uses them: each at its quadrature nodes.
      const pts = [], dwl = [];
      for (const side of [1, -1]) {
        let prev = null;
        for (const [x, o] of h.stations) {
          for (let k = 0; k + 1 < o.length; k++) pts.push(...P(h, x, o[k][1], side * o[k][0]), ...P(h, x, o[k + 1][1], side * o[k + 1][0]));
          if (o.length && o[0][1] === 0) { const p = P(h, x, 0, side * o[0][0]); if (prev) dwl.push(...prev, ...p); prev = p; } else prev = null;
        }
      }
      layers.stations.add(seg(pts, 0x5b2a4a, 0.85));
      layers.stations.add(seg(dwl, 0x1f6feb));
      // Keel: each section's lowest point, on the centreplane.
      const kp = h.keel.flatMap(([x, d]) => P(h, x, d, 0));
      const kg = new THREE.BufferGeometry(); kg.setAttribute("position", new THREE.Float32BufferAttribute(kp, 3));
      layers.keel.add(new THREE.Line(kg, new THREE.LineBasicMaterial({ color: 0xd62728, depthTest: false })));
      // CAD ray hits the stations were interpolated from.
      if (h.rays) {
        const rp = [];
        for (const [x, o] of h.rays) for (const [y, z] of o) rp.push(...P(h, x, z, y));
        const rg = new THREE.BufferGeometry(); rg.setAttribute("position", new THREE.Float32BufferAttribute(rp, 3));
        layers.rays.add(new THREE.Points(rg, new THREE.PointsMaterial({ color: 0xe8590c, size: 3, sizeAttenuation: false, depthTest: false })));
      }
      const xs = h.stations.map((s) => s[0]);
      box.expandByPoint(new THREE.Vector3(Math.min(...xs) + h.placement.x, -h.draft, h.placement.y - h.beam / 2));
      box.expandByPoint(new THREE.Vector3(Math.max(...xs) + h.placement.x, 0, h.placement.y + h.beam / 2));
    }
    // Room aft for a closure's hollow.
    if (hs.some((h) => h.transom)) box.min.x -= 0.1 * (box.max.x - box.min.x);
    bounds = box;
    hullBounds = box.clone();
    if (hs.length) empty.style.display = "none";
    // The cut's layers sit at the waterline they were cut at.
    show(mode);
    setView(view);
  }

  // The virtual appendage behind each transom: the transom section with its
  // half-beam scaled by φ(s) over the hollow, x = x_T − s·L_v, where
  // `hollowLength(depth)` gives L_v (null: no hollow).
  function closure(hollowLength) {
    clear(layers.closure);
    for (const h of cutHulls) {
      const t = h.transom;
      if (!t || !hollowLength) continue;
      const lv = hollowLength(t.depth);
      if (lv == null || lv <= 0 || !t.outline.length) continue;
      const ns = 25, o = t.outline;
      layers.closure.add(sheet(h, ns, o.length, (i, j) => {
        const s = i / (ns - 1);
        return [t.x - s * lv, o[j][0] * phi(s), o[j][1]];
      }, 0xe8590c, 0.55));
      const lines = [];
      for (const s of [0.25, 0.5, 0.75]) for (const side of [1, -1]) {
        for (let j = 0; j + 1 < o.length; j++) {
          lines.push(...P(h, t.x - s * lv, o[j][1], side * o[j][0] * phi(s)), ...P(h, t.x - s * lv, o[j + 1][1], side * o[j + 1][0] * phi(s)));
        }
      }
      layers.closure.add(seg(lines, 0xb54708, 0.8));
    }
    render();
  }

  // Coloured mesh from positions (scene frame) and per-vertex values.
  function colored(pos, index, vals, vmax, userData, opacity = 1, lit = true) {
    const col = new Float32Array(vals.length * 3);
    vals.forEach((v, i) => col.set(diverging(v / vmax), 3 * i));
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    geo.setAttribute("color", new THREE.BufferAttribute(col, 3));
    geo.setIndex(Array.isArray(index) ? index : new THREE.BufferAttribute(index, 1));
    geo.computeVertexNormals();
    // Unlit when the colour is the whole story: at true scale a surface's
    // shading would draw only its tiniest ripples, as streaks.
    const Mat = lit ? THREE.MeshStandardMaterial : THREE.MeshBasicMaterial;
    const mesh = new THREE.Mesh(geo, new Mat({ vertexColors: true, ...(lit ? { roughness: 0.8 } : {}), side: THREE.DoubleSide,
      transparent: opacity < 1, opacity, depthWrite: opacity >= 1 }));
    mesh.userData = { vals, ...userData };
    return mesh;
  }

  // A case's flow (the stored answer), or none: the pressure on each hull,
  // and the free surface as a height field.
  function flow(d) {
    clear(layers.pressure); clear(layers.waves); clear(layers.hull);
    probes = [];
    const has = !!d;
    q(".legend.cp").style.display = q(".legend.zeta").style.display = has ? "block" : "none";
    show(has ? "flow" : "cut");
    views.querySelector("button[data-view=wake]").hidden = !has;
    if (!has) { wakeBounds = null; bounds = hullBounds; if (view === "wake") pick("iso"); else setView(view, true); return; }
    const cpMax = range99(d.hulls.flatMap((p) => p.cp));
    q(".legend.cp .lo").textContent = `−${cpMax.toPrecision(2)}`;
    q(".legend.cp .hi").textContent = `+${cpMax.toPrecision(2)}`;
    d.hulls.forEach((p, k) => {
      const nx = p.x.length, nz = p.depth.length;
      for (const side of [1, -1]) {
        const pos = new Float32Array(nx * nz * 3);
        for (let i = 0; i < nx; i++) for (let j = 0; j < nz; j++) {
          const n = i * nz + j;
          pos.set([p.x[i], -p.depth[j], p.y + side * p.half_beam[n]], 3 * n);
        }
        const index = [];
        for (let i = 0; i + 1 < nx; i++) for (let j = 0; j + 1 < nz; j++) {
          const a = i * nz + j, b = a + nz;
          // Rows below the keel collapse onto the centreplane: not hull.
          if (!(p.half_beam[a] || p.half_beam[b] || p.half_beam[a + 1] || p.half_beam[b + 1])) continue;
          index.push(a, b, b + 1, a, b + 1, a + 1);
        }
        const m = colored(pos, index, p.cp, cpMax, { what: "cp", hull: k });
        layers.pressure.add(m); probes.push(m);
      }
    });
    // The whole hull at the attitude, topsides included, cut at the water
    // while the pressure is shown below it.
    const hb = new THREE.Box3();
    for (const p of d.hulls) {
      const geo = meshGeometry(p.mesh);
      hb.union(geo.boundingBox);
      layers.hull.add(new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ color: 0xe9e5db, roughness: 0.7, side: THREE.DoubleSide, clippingPlanes: [aboveWater] })));
    }

    const g = d.surface;
    const zeta = typeof g.zeta === "string" ? new Float32Array(b64(g.zeta)) : g.zeta;
    const zmax = range99(zeta);
    const mm = (1000 * zmax) >= 100 ? (1000 * zmax).toFixed(0) : (1000 * zmax).toPrecision(2);
    q(".legend.zeta .lo").textContent = `−${mm} mm`;
    q(".legend.zeta .hi").textContent = `+${mm} mm`;
    // The waterplanes at this attitude, from each pressure mesh's top row
    // (its depth-0 half-beams): inside them the surface is flattened to calm
    // water and left to the hull to cover, rather than cut out cell by cell.
    const wl = d.hulls.map((p) => {
      const nz = p.depth.length;
      return { x: p.x, y: p.y, hb: p.x.map((_, i) => p.half_beam[i * nz]) };
    });
    const insideWl = (x, y) => wl.some((w) => {
      const n = w.x.length;
      if (x < w.x[0] || x > w.x[n - 1]) return false;
      let i = 1;
      while (i < n - 1 && w.x[i] < x) i++;
      const t = (x - w.x[i - 1]) / (w.x[i] - w.x[i - 1] || 1);
      return Math.abs(y - w.y) < w.hb[i - 1] + t * (w.hb[i] - w.hb[i - 1]);
    });
    const pos = new Float32Array(g.nx * g.ny * 3), vals = new Float32Array(g.nx * g.ny);
    for (let iy = 0; iy < g.ny; iy++) for (let ix = 0; ix < g.nx; ix++) {
      const n = iy * g.nx + ix;
      const x = g.x0 + (g.x1 - g.x0) * ix / (g.nx - 1), y = g.y0 + (g.y1 - g.y0) * iy / (g.ny - 1);
      const z = insideWl(x, y) ? 0 : zeta[n];
      vals[n] = z;
      pos.set([x, z, y], 3 * n);
    }
    const index = new Uint32Array(6 * (g.nx - 1) * (g.ny - 1));
    let m = 0;
    for (let iy = 0; iy + 1 < g.ny; iy++) for (let ix = 0; ix + 1 < g.nx; ix++) {
      const a = iy * g.nx + ix, b = a + 1, c = a + g.nx + 1, e = a + g.nx;
      index.set([a, c, b, a, e, c], m); m += 6;
    }
    // See-through, so the hull's pressure shows beneath it.
    const waveMesh = colored(pos, index, vals, zmax, { what: "zeta" }, 0.72, false);
    waveMesh.renderOrder = 1;
    layers.waves.add(waveMesh); probes.push(waveMesh);
    // The Wake view frames the whole pattern; the others stay on the hull.
    const first = !wakeBounds;
    bounds = (hullBounds ? hullBounds.clone() : new THREE.Box3()).union(hb);
    hullBounds = hullBounds || hb.clone();
    wakeBounds = bounds.clone().union(new THREE.Box3(new THREE.Vector3(g.x0, -0.01, g.y0), new THREE.Vector3(g.x1, 0.01, g.y1)));
    empty.style.display = "none";
    if (first) pick("wake"); else render();
  }

  let fitW = 1, fitH = 1;
  function setView(name, keepDir = false) {
    view = name;
    if (!bounds && !geomBounds) return render();
    const box = (mode === "cut" && geomBounds ? geomBounds
      : name === "wake" && wakeBounds ? wakeBounds : bounds).clone();
    const mid = box.getCenter(new THREE.Vector3()), size = box.getSize(new THREE.Vector3());
    const dirs = { wake: new THREE.Vector3(0.9, 0.7, 0.9), iso: new THREE.Vector3(0.8, 0.55, 1.0), profile: new THREE.Vector3(0, 0, 1), plan: new THREE.Vector3(0, 1, 0), body: new THREE.Vector3(1, 0, 0) };
    const dir = keepDir ? camera.position.clone().sub(controls.target).normalize() : dirs[name].clone().normalize();
    const r = size.length();
    controls.target.copy(mid);
    camera.position.copy(mid).addScaledVector(dir, r * 2);
    camera.up.set(0, 1, 0);
    if (name === "plan" && !keepDir) camera.up.set(0, 0, -1);
    camera.near = -r * 10; camera.far = r * 10;
    camera.lookAt(mid);
    camera.updateMatrixWorld();
    const inv = camera.matrixWorldInverse, seen = new THREE.Box3();
    for (let i = 0; i < 8; i++) seen.expandByPoint(new THREE.Vector3(i & 1 ? box.max.x : box.min.x, i & 2 ? box.max.y : box.min.y, i & 4 ? box.max.z : box.min.z).applyMatrix4(inv));
    const ss = seen.getSize(new THREE.Vector3());
    fitW = ss.x * 1.12; fitH = ss.y * 1.12;
    camera.zoom = 1;
    resize();
    controls.update();
  }

  function resize() {
    const w = host.clientWidth, h = host.clientHeight;
    if (!w || !h) return;
    renderer.setSize(w, h);
    const aspect = w / h;
    let hw = fitW / 2, hh = fitH / 2;
    if (hw / hh > aspect) hh = hw / aspect; else hw = hh * aspect;
    Object.assign(camera, { left: -hw, right: hw, top: hh, bottom: -hh });
    camera.updateProjectionMatrix();
    render();
  }
  new ResizeObserver(resize).observe(host);

  let pending = false;
  function render() {
    if (pending) return;
    pending = true;
    requestAnimationFrame(() => { pending = false; if (controls.update()) render(); renderer.render(scene, camera); });
  }

  const ray = new THREE.Raycaster(), ndc = new THREE.Vector2();
  renderer.domElement.addEventListener("pointermove", (e) => {
    const targets = [...(layers.surface.visible ? surfaces : []),
      ...probes.filter((m) => layers[m.userData.what === "cp" ? "pressure" : "waves"].visible)];
    if (!targets.length) { readout.style.display = "none"; return; }
    const r = renderer.domElement.getBoundingClientRect();
    ndc.set(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    const hit = ray.intersectObjects(targets, false)[0];
    if (!hit) { readout.style.display = "none"; return; }
    const ud = hit.object.userData;
    if (ud.what) {
      // Barycentric interpolation of the mesh's values at the hit.
      const f = hit.face, bc = new THREE.Vector3();
      const Pa = hit.object.geometry.attributes.position, A = new THREE.Vector3(), B = new THREE.Vector3(), C = new THREE.Vector3();
      A.fromBufferAttribute(Pa, f.a); B.fromBufferAttribute(Pa, f.b); C.fromBufferAttribute(Pa, f.c);
      const local = hit.object.worldToLocal(hit.point.clone());
      THREE.Triangle.getBarycoord(local, A, B, C, bc);
      const v = bc.x * ud.vals[f.a] + bc.y * ud.vals[f.b] + bc.z * ud.vals[f.c];
      readout.textContent = ud.what === "cp"
        ? `x ${local.x.toFixed(3)} m   depth ${(-local.y).toFixed(3)} m   Cp ${v.toFixed(4)}`
        : `x ${local.x.toFixed(2)} m   y ${local.z.toFixed(2)} m   ζ ${(1000 * v).toFixed(2)} mm`;
      readout.style.display = "block";
      return;
    }
    const p = model.worldToLocal(hit.point.clone()), { h } = hit.object.userData;
    readout.textContent = `x ${(p.x - h.placement.x).toFixed(3)} m   z ${(-p.y).toFixed(3)} m   half-beam ${Math.abs(p.z - h.placement.y).toFixed(4)} m`;
    readout.style.display = "block";
  });
  renderer.domElement.addEventListener("pointerleave", () => { readout.style.display = "none"; });

  function pick(name) {
    views.querySelectorAll("button").forEach((x) => x.classList.toggle("on", x.dataset.view === name));
    setView(name);
  }
  for (const b of views.querySelectorAll("button")) b.onclick = () => pick(b.dataset.view);

  // A study in waves: the whole hull at its attitude (z up from the water),
  // over which `seaway` moves it.
  function hullsAt(meshes) {
    clear(layers.hull);
    const hb = new THREE.Box3();
    for (const m of meshes) {
      const geo = meshGeometry(m);
      hb.union(geo.boundingBox);
      layers.hull.add(new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ color: 0xe9e5db, roughness: 0.7, side: THREE.DoubleSide })));
    }
    bounds = hullBounds = hb;
    show("sea");
    empty.style.display = "none";
    setView(view);
  }

  // A regular wave and the hull moving in it: the incident elevation
  //   ζ = A cos(k((x − x_G) cos β + y sin β) − ω_e t)
  // in the frame moving with the boat, and each motion A·Re[η̂ e^{−iω_e t}]
  // about G — sway to port, heave up, roll port side up, pitch bow up, yaw
  // bow to port — from the study's complex responses (per unit amplitude).
  // `spec`: { k, heading [rad], omega_e, amp [m], g: [x, y, z], eta:
  // { sway, heave, roll, pitch, yaw } as [re, im], extent: [x0, x1, y0, y1],
  // rate }, or null to stop. Scene axes are (x, z up, y): SWAP swaps them.
  const SWAP = new THREE.Matrix4().set(1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1);
  let sea = null;
  function seaway(spec) {
    if (sea) { cancelAnimationFrame(sea.raf); sea = null; }
    clear(layers.seaway);
    layers.hull.matrixAutoUpdate = true;
    layers.hull.position.set(0, 0, 0); layers.hull.rotation.set(0, 0, 0); layers.hull.updateMatrix();
    if (!spec) { render(); return; }
    const [x0, x1, y0, y1] = spec.extent;
    // About 16 points a wavelength, within a budget.
    const lam = 2 * Math.PI / spec.k, step = Math.max(lam / 16, Math.sqrt((x1 - x0) * (y1 - y0) / 60000));
    const nx = Math.max(2, Math.ceil((x1 - x0) / step) + 1), ny = Math.max(2, Math.ceil((y1 - y0) / step) + 1);
    const pos = new Float32Array(nx * ny * 3), phase = new Float32Array(nx * ny);
    const cb = Math.cos(spec.heading), sb = Math.sin(spec.heading);
    for (let iy = 0; iy < ny; iy++) for (let ix = 0; ix < nx; ix++) {
      const n = iy * nx + ix, x = x0 + (x1 - x0) * ix / (nx - 1), y = y0 + (y1 - y0) * iy / (ny - 1);
      pos.set([x, 0, y], 3 * n);
      phase[n] = spec.k * ((x - spec.g[0]) * cb + y * sb);
    }
    const index = new Uint32Array(6 * (nx - 1) * (ny - 1));
    let m = 0;
    for (let iy = 0; iy + 1 < ny; iy++) for (let ix = 0; ix + 1 < nx; ix++) {
      const a = iy * nx + ix, b = a + 1, c = a + nx + 1, e = a + nx;
      index.set([a, c, b, a, e, c], m); m += 6;
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    // Shaded by elevation, troughs dark and crests light: at a real wave's
    // steepness the lighting alone barely shows it.
    const col = new Float32Array(nx * ny * 3);
    geo.setAttribute("color", new THREE.BufferAttribute(col, 3));
    geo.setIndex(new THREE.BufferAttribute(index, 1));
    const water = new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.4, metalness: 0.05,
      transparent: true, opacity: 0.7, side: THREE.DoubleSide, depthWrite: false }));
    water.renderOrder = 1;
    layers.seaway.add(water);
    // Walls down the patch's edges, so side-on views show the wave's
    // profile: each edge point of the grid, and below it one at the walls'
    // depth.
    const edge = [];
    for (let ix = 0; ix < nx; ix++) edge.push(ix);
    for (let iy = 1; iy < ny; iy++) edge.push(iy * nx + nx - 1);
    for (let ix = nx - 2; ix >= 0; ix--) edge.push((ny - 1) * nx + ix);
    for (let iy = ny - 2; iy >= 0; iy--) edge.push(iy * nx);
    const depth = Math.min(0.35 * lam, 0.3 * (x1 - x0));
    const wpos = new Float32Array(edge.length * 6), widx = [];
    edge.forEach((n, i) => {
      wpos.set([pos[3 * n], 0, pos[3 * n + 2], pos[3 * n], -depth, pos[3 * n + 2]], 6 * i);
      if (i + 1 < edge.length) { const a = 2 * i; widx.push(a, a + 1, a + 2, a + 1, a + 3, a + 2); }
    });
    const wgeo = new THREE.BufferGeometry();
    wgeo.setAttribute("position", new THREE.BufferAttribute(wpos, 3));
    wgeo.setIndex(widx);
    const walls = new THREE.Mesh(wgeo, new THREE.MeshBasicMaterial({ color: 0x2f74c0, transparent: true, opacity: 0.28,
      side: THREE.DoubleSide, depthWrite: false }));
    walls.renderOrder = 1;
    layers.seaway.add(walls);
    sea = { spec, t: 0, last: null, raf: 0 };
    const DARK = [0.09, 0.3, 0.55], LIGHT = [0.62, 0.8, 0.95];
    const re = (z, c, s) => z[0] * c + z[1] * s; // Re[(a + ib) e^{−iωt}]
    const R = new THREE.Matrix4(), T = new THREE.Matrix4(), tmp = new THREE.Matrix4();
    const tick = (now) => {
      if (!sea || sea.spec !== spec) return;
      if (sea.last != null) sea.t += Math.min(0.1, (now - sea.last) / 1000) * spec.rate;
      sea.last = now;
      const wt = spec.omega_e * sea.t, c = Math.cos(wt), s = Math.sin(wt), A = spec.amp;
      const attr = geo.attributes.position, ca = geo.attributes.color.array;
      for (let n = 0; n < phase.length; n++) {
        const c1 = Math.cos(phase[n] - wt), u = 0.5 + 0.5 * c1;
        attr.array[3 * n + 1] = A * c1;
        for (let j = 0; j < 3; j++) ca[3 * n + j] = DARK[j] + u * (LIGHT[j] - DARK[j]);
      }
      attr.needsUpdate = true;
      geo.attributes.color.needsUpdate = true;
      geo.computeVertexNormals();
      const wa = wgeo.attributes.position;
      edge.forEach((n, i) => { wa.array[6 * i + 1] = attr.array[3 * n + 1]; });
      wa.needsUpdate = true;
      const e = spec.eta, [gx, gy, gz] = spec.g;
      const d = { sway: A * re(e.sway, c, s), heave: A * re(e.heave, c, s), roll: A * re(e.roll, c, s),
        pitch: A * re(e.pitch, c, s), yaw: A * re(e.yaw, c, s) };
      // Hull frame (x fwd, y port, z up): about G, then into the scene.
      R.makeRotationZ(d.yaw).multiply(tmp.makeRotationY(-d.pitch)).multiply(tmp.makeRotationX(d.roll));
      T.makeTranslation(gx, gy + d.sway, gz + d.heave).multiply(R).multiply(tmp.makeTranslation(-gx, -gy, -gz));
      layers.hull.matrixAutoUpdate = false;
      layers.hull.matrix.copy(SWAP).multiply(T).multiply(SWAP);
      layers.hull.matrixWorldNeedsUpdate = true;
      controls.update();
      renderer.render(scene, camera);
      sea.raf = requestAnimationFrame(tick);
    };
    sea.raf = requestAnimationFrame(tick);
  }

  return {
    geometry, setWater, load, closure, flow, pick, hullsAt, seaway,
    // "cut" (the hull and its sections), "flow" (a calm-water study) or "sea"
    // (a study in waves).
    mode: (m) => { show(m); setView(view, true); },
    // Say why there is nothing to show.
    empty: (msg) => { empty.textContent = msg; empty.style.display = ""; },
  };
}

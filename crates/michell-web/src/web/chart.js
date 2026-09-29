// Line and point charts in SVG, one y-axis each: series in the palette's
// fixed order (by entity, never by rank), a crosshair that snaps to the
// nearest x with every series' value in one tooltip, a legend for two or
// more series with direct end labels for up to four, and a table view.
//
//   lineChart(host, {
//     title, x: { label, fmt }, y: { label, fmt, zero },
//     series: [{ name, points: [[x, y, meta?], ...], dashed?, color? }],
//     onPick: (meta) => …,   // a point clicked
//     marks: "auto" | "always" | "never",
//   })
//
// Series names and labels are set with textContent, never as HTML.

const SVGNS = "http://www.w3.org/2000/svg";
const el = (tag, attrs = {}, parent) => {
  const e = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
};
const h = (tag, cls, parent, text) => {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text != null) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
};

// The categorical slot a series takes: its given colour, else by its index.
export const slot = (i) => `var(--series-${(i % 8) + 1})`;

// About five round ticks spanning [a, b].
export function ticks(a, b, n = 5) {
  const span = b - a || Math.abs(a) || 1, raw = span / n, mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw);
  const out = [];
  for (let t = Math.ceil(a / step - 1e-9) * step; t <= b + 1e-9 * span; t += step) out.push(+t.toPrecision(12));
  return out;
}

const plain = (v) => (Math.abs(v) >= 1000 || (Math.abs(v) < 0.01 && v !== 0) ? v.toPrecision(3) : +v.toPrecision(4)).toString();

export function lineChart(host, spec) {
  host.innerHTML = "";
  host.classList.add("chart");
  // Drawn at the host's own width, so text stays its size; redrawn when
  // that changes.
  let W = 640;
  const H = spec.height || 260, L = 56, R = spec.series.length >= 2 && spec.series.length <= 4 ? 110 : 16, T = 10, B = 40;
  const xf = spec.x?.fmt || plain, yf = spec.y?.fmt || plain;
  const series = spec.series.map((s, i) => ({ ...s, color: s.color || slot(i), points: s.points.filter((p) => Number.isFinite(p[0]) && Number.isFinite(p[1])) }));
  const all = series.flatMap((s) => s.points);

  const bar = h("div", "bar", host);
  if (spec.title) h("span", "title", bar, spec.title);
  h("span", "grow", bar);
  const toggle = h("button", "", bar, "Table");
  const body = h("div", "", host);
  const legend = h("div", "legend-row", host);
  const tip = h("div", "tip", host);
  let showTable = false;
  toggle.onclick = () => { showTable = !showTable; toggle.textContent = showTable ? "Chart" : "Table"; draw(); };

  function table() {
    body.innerHTML = "";
    const wrap = h("div", "tablewrap", body), t = h("table", "", wrap), head = h("tr", "", h("thead", "", t));
    h("th", "", head, spec.x?.label || "x");
    for (const s of series) h("th", "num", head, s.name);
    const xs = [...new Set(all.map((p) => p[0]))].sort((a, b) => a - b);
    const tb = h("tbody", "", t);
    for (const x of xs) {
      const tr = h("tr", "", tb);
      h("td", "num", tr, xf(x));
      for (const s of series) {
        const p = s.points.find((q) => q[0] === x);
        h("td", "num", tr, p ? yf(p[1]) : "");
      }
    }
  }

  function draw() {
    legend.innerHTML = "";
    if (!all.length) { body.innerHTML = ""; h("div", "empty", body, spec.empty || "Nothing to plot yet"); return; }
    if (showTable) return table();
    body.innerHTML = "";
    W = Math.max(320, Math.round(host.clientWidth || 640));
    let [x0, x1] = [Math.min(...all.map((p) => p[0])), Math.max(...all.map((p) => p[0]))];
    let [y0, y1] = [Math.min(...all.map((p) => p[1])), Math.max(...all.map((p) => p[1]))];
    if (spec.y?.zero !== false) { y0 = Math.min(0, y0); y1 = Math.max(0, y1); }
    if (x1 - x0 < 1e-12) { x0 -= 0.5 * (Math.abs(x0) || 1); x1 += 0.5 * (Math.abs(x1) || 1); }
    if (y1 - y0 < 1e-12) { y0 -= 0.5 * (Math.abs(y0) || 1); y1 += 0.5 * (Math.abs(y1) || 1); }
    const pad = 0.06 * (y1 - y0); if (y0 < 0 || spec.y?.zero === false) y0 -= pad; y1 += pad;
    const X = (x) => L + (W - L - R) * (x - x0) / (x1 - x0), Y = (y) => H - B - (H - T - B) * (y - y0) / (y1 - y0);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img" }, body);
    if (spec.title) el("title", {}, svg).textContent = spec.title;
    for (const t of ticks(y0, y1)) {
      el("line", { class: "grid", x1: L, x2: W - R, y1: Y(t), y2: Y(t) }, svg);
      el("text", { x: L - 6, y: Y(t) + 4, "text-anchor": "end" }, svg).textContent = yf(t);
    }
    if (y0 < 0 && y1 > 0) el("line", { class: "axis", x1: L, x2: W - R, y1: Y(0), y2: Y(0) }, svg);
    el("line", { class: "axis", x1: L, x2: W - R, y1: H - B, y2: H - B }, svg);
    for (const t of ticks(x0, x1, 6)) el("text", { x: X(t), y: H - B + 15, "text-anchor": "middle" }, svg).textContent = xf(t);
    if (spec.x?.label) el("text", { x: (L + W - R) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = spec.x.label;
    if (spec.y?.label) el("text", { x: 12, y: (T + H - B) / 2, "text-anchor": "middle", transform: `rotate(-90 12 ${(T + H - B) / 2})` }, svg).textContent = spec.y.label;

    const marks = spec.marks || "auto";
    series.forEach((s) => {
      const ps = [...s.points].sort((a, b) => a[0] - b[0]);
      if (ps.length > 1 && !s.scatter) el("polyline", { fill: "none", stroke: s.color, "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round", "stroke-dasharray": s.dashed ? "5 4" : "none", points: ps.map((p) => `${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join(" ") }, svg);
      const dots = marks === "always" || s.scatter || (marks === "auto" && ps.length <= 30);
      if (dots) for (const p of ps) {
        const c = el("circle", { class: "dot", cx: X(p[0]).toFixed(1), cy: Y(p[1]).toFixed(1), r: 4, fill: s.color }, svg);
        if (spec.onPick && p[2] != null) {
          c.style.cursor = "pointer";
          // A hit area bigger than the mark.
          const hit = el("circle", { cx: X(p[0]).toFixed(1), cy: Y(p[1]).toFixed(1), r: 12, fill: "transparent" }, svg);
          hit.style.cursor = "pointer";
          hit.addEventListener("click", () => spec.onPick(p[2]));
        }
      }
      // Direct label at the line's end, for up to four series.
      if (series.length >= 2 && series.length <= 4 && ps.length) {
        const last = ps[ps.length - 1];
        el("text", { class: "direct", x: X(last[0]) + 8, y: Y(last[1]) + 4 }, svg).textContent = s.name;
      }
    });
    if (series.length >= 2) {
      for (const s of series) {
        const item = h("span", "", legend);
        const key = h("i", "", item);
        key.style.borderTopColor = s.color;
        if (s.dashed) key.style.borderTopStyle = "dashed";
        item.appendChild(document.createTextNode(s.name));
      }
    }

    // The crosshair: the nearest x, every series' value there.
    const cross = el("line", { class: "cross", y1: T, y2: H - B, visibility: "hidden" }, svg);
    const xs = [...new Set(all.map((p) => p[0]))].sort((a, b) => a - b);
    const overlay = el("rect", { x: L, y: T, width: W - L - R, height: H - T - B, fill: "transparent" }, svg);
    overlay.style.pointerEvents = "all";
    // The point hit areas stay on top of the overlay.
    for (const hit of svg.querySelectorAll("circle[fill=transparent]")) svg.appendChild(hit);
    const move = (ev) => {
      const r = svg.getBoundingClientRect(), sx = (ev.clientX - r.left) * W / r.width;
      const xv = x0 + (sx - L) / (W - L - R) * (x1 - x0);
      const near = xs.reduce((a, b) => Math.abs(b - xv) < Math.abs(a - xv) ? b : a, xs[0]);
      cross.setAttribute("x1", X(near)); cross.setAttribute("x2", X(near)); cross.setAttribute("visibility", "visible");
      tip.innerHTML = "";
      h("div", "x", tip, `${spec.x?.label || "x"}: ${xf(near)}`);
      for (const s of series) {
        const p = s.points.find((q) => Math.abs(q[0] - near) < 1e-9 * (Math.abs(near) + 1));
        if (!p) continue;
        const row = h("div", "row", tip);
        const key = h("i", "", row); key.style.borderTopColor = s.color;
        h("b", "", row, yf(p[1]));
        if (series.length > 1 || s.name) h("span", "", row, s.name);
      }
      tip.style.display = "block";
      const hr = host.getBoundingClientRect(), left = ev.clientX - hr.left + 14;
      tip.style.left = `${Math.min(left, hr.width - tip.offsetWidth - 4)}px`;
      tip.style.top = `${ev.clientY - hr.top + 12}px`;
    };
    for (const t of [overlay, ...svg.querySelectorAll("circle[fill=transparent]")]) {
      t.addEventListener("pointermove", move);
      t.addEventListener("pointerleave", () => { tip.style.display = "none"; cross.setAttribute("visibility", "hidden"); });
    }
  }
  draw();
  let lastW = host.clientWidth, pending = null;
  new ResizeObserver(() => {
    if (Math.abs(host.clientWidth - lastW) < 4) return;
    lastW = host.clientWidth;
    clearTimeout(pending);
    pending = setTimeout(() => { if (!showTable) draw(); }, 80);
  }).observe(host);
  return { redraw: draw };
}

// What every page shares: small helpers, the API client, and the top bar.

export const $ = (id) => document.getElementById(id);
export const G = 9.81;
export const KN = 0.514444;
// The solver's fluid: salt water at 15 °C (michell::Fluid::SEAWATER_15C).
export const RHO = 1025.9;

export const fmt = (v, d = 3) => v == null || !Number.isFinite(Number(v)) ? "—" : Number(v).toFixed(d);
export const fmtBytes = (n) => n < 1024 ? `${n} B` : n < 1 << 20 ? `${(n / 1024).toFixed(0)} kB` : `${(n / (1 << 20)).toFixed(1)} MB`;
export const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

// A unix time as a short local date and time; today's as the time alone.
export function fmtTime(t) {
  if (t == null) return "—";
  const d = new Date(t * 1000), now = new Date();
  const time = d.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
  return d.toDateString() === now.toDateString() ? time
    : `${d.toLocaleDateString([], { day: "numeric", month: "short", ...(d.getFullYear() === now.getFullYear() ? {} : { year: "numeric" }) })} ${time}`;
}
export function fmtDuration(s) {
  if (s == null) return "—";
  return s < 60 ? `${s.toFixed(s < 10 ? 1 : 0)} s` : s < 3600 ? `${Math.floor(s / 60)} min ${Math.round(s % 60)} s` : `${(s / 3600).toFixed(1)} h`;
}

// Base64 to an ArrayBuffer (the server's little-endian typed arrays).
export function b64(s) {
  const bin = atob(s), out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out.buffer;
}

// The name this browser's uploads and requests are labelled with.
export function whoami() {
  try { return localStorage.getItem("michell.by") || ""; } catch { return ""; }
}
function setWhoami(v) {
  try { localStorage.setItem("michell.by", v); } catch { /* private window */ }
}

// A JSON API call: the parsed answer, or an Error with the server's message.
// With a `json` or raw `body` it is a POST unless `method` says otherwise.
export async function api(path, { method, body, json } = {}) {
  const opts = {};
  if (json !== undefined) { opts.body = JSON.stringify(json); opts.headers = { "Content-Type": "application/json" }; }
  else if (body !== undefined) opts.body = body;
  opts.method = method || (opts.body !== undefined ? "POST" : "GET");
  const r = await fetch(path, opts);
  const text = await r.text();
  let data;
  try { data = text ? JSON.parse(text) : null; } catch { data = { error: text || r.statusText }; }
  if (!r.ok) throw new Error((data && data.error) || `${r.status} ${r.statusText}`);
  return data;
}

// The id in a path like /hulls/12 (the last number in it).
export function pathId() {
  const m = location.pathname.match(/\/(\d+)\/?$/);
  return m ? Number(m[1]) : null;
}

// A hull's summary in words: principal dimensions and displacement.
export function hullDims(summary) {
  const hs = summary?.hulls || [];
  if (!hs.length) return "";
  const h = hs[0];
  const vol = hs.reduce((s, h) => s + (h.displaced_volume || 0), 0);
  const n = hs.length > 1 ? `${hs.length} hulls · ` : "";
  return `${n}L ${fmt(h.length, 2)} · B ${fmt(h.beam, 2)} · T ${fmt(h.draft, 3)} m · Δ ${fmt(RHO * vol, 0)} kg`;
}

// The reference length speeds are made Froude numbers on: the longest hull.
export const lRef = (summary) => Math.max(0, ...(summary?.hulls || []).map((h) => h.length || 0));
// A Froude number on length `l` as knots.
export const knots = (fn, l) => fn * Math.sqrt(G * l) / KN;

// A case's parameters in words; `l` (the hull's reference length) adds the
// speed in knots.
export function caseLabel(p, l) {
  const speed = `Fn ${fmt(p.froude, 3)}${l ? ` (${fmt(knots(p.froude, l), 2)} kn)` : ""}`;
  const how = { sinking: "", scale: " scaled xyz", scale_yz: " scaled yz" }[p.mass_by] ?? "";
  const load = p.mass == null && p.lcg == null ? "design load"
    : `${p.mass == null ? "design mass" : `${fmt(p.mass, 1)} kg${how}`}${p.lcg == null ? "" : ` at x ${fmt(p.lcg, 3)}`}`;
  const layout = p.spans ? `cat spans ${p.spans.map((s) => fmt(s, 2)).join(", ")} m (held)`
    : p.span != null ? `cat span ${fmt(p.span, 2)} m` : "";
  const c = p.closure || {};
  const closure = c.type === "off" ? "no closure" : c.type === "fixed" ? `hollow ${fmt(c.length, 2)} m` : "";
  return [speed, layout, load, p.dynamic === false ? "design attitude" : "", closure].filter(Boolean).join(" · ");
}

const TABS = [
  ["hulls", "/hulls", "Hulls"],
  ["new-cases", "/cases/new", "New cases"],
  ["queue", "/queue", "Queue"],
  ["results", "/results", "Results"],
];

// The top bar: the pages, the queue's length, and who is asking.
export function nav(on) {
  const bar = document.createElement("nav");
  bar.className = "top";
  bar.innerHTML = `<a class="brand" href="/hulls">michell</a>`
    + TABS.map(([k, href, label]) => `<a class="tab${k === on ? " on" : ""}" href="${href}" data-tab="${k}">${label}${k === "queue" ? `<span class="badge" id="qbadge" hidden></span>` : ""}</a>`).join("")
    + `<label class="who">you <input id="whoami" placeholder="your name" autocomplete="name"></label>`;
  document.body.prepend(bar);
  const who = bar.querySelector("#whoami");
  who.value = whoami();
  who.addEventListener("change", () => setWhoami(who.value.trim()));
  // The queue's length, kept up to date while the page is open.
  const badge = bar.querySelector("#qbadge");
  async function poll() {
    try {
      const q = await api("/api/queue");
      const n = q.cases.length;
      badge.hidden = !n;
      badge.textContent = n;
    } catch { /* offline: leave it */ }
  }
  poll();
  setInterval(() => { if (!document.hidden) poll(); }, 5000);
}

// A status line: a message, or an error in red.
export function status(msg, err = false) {
  const s = $("status");
  if (!s) return;
  s.textContent = msg;
  s.className = err ? "err" : "";
}

// Each hull's hydrostatics, as the cut reports them.
export function hullCards(hulls, notes = []) {
  return notes.map((n) => `<div class="note">${esc(n)}</div>`).join("") + hulls.map((h, i) => {
    const title = hulls.length > 1 ? `Hull ${i + 1}` : "Hull";
    const where = h.placement.y || h.placement.x ? `y ${fmt(h.placement.y)} m` : "";
    const t = h.transom;
    return `<div class="hull">
      <h3>${title}<span class="tag">${esc(where)}</span></h3>
      <dl>
        <dt>Length</dt><dd>${fmt(h.length)} m</dd>
        <dt>Beam (WL)</dt><dd>${fmt(h.beam)} m</dd>
        <dt>Draft</dt><dd>${fmt(h.draft)} m</dd>
        <dt>Displaced volume</dt><dd>${fmt(h.displaced_volume, 4)} m³</dd>
        <dt>Displacement</dt><dd>${fmt(RHO * h.displaced_volume, 1)} kg</dd>
        <dt>Wetted surface</dt><dd>${fmt(h.wetted_surface)} m²</dd>
        <dt>Waterplane area</dt><dd>${fmt(h.waterplane_area)} m²</dd>
        <dt>LCB</dt><dd>x ${fmt(h.lcb_x)} m</dd>
        <dt>Stations</dt><dd>${h.stations.length}</dd>
        ${t ? `<dt>Transom</dt><dd>x ${fmt(t.x)} m, ${(100 * t.area_ratio).toFixed(1)}% of max section</dd>` : ""}
      </dl>
      <details><summary>Sectioning diagnostics</summary><pre>${esc(h.diagnostics.join("\n"))}</pre></details>
    </div>`;
  }).join("");
}

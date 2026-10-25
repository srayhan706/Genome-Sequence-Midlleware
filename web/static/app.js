"use strict";

// ---------------------------------------------------------------- helpers
const $ = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];
const nf = new Intl.NumberFormat();
const fmtInt = (n) => nf.format(Math.round(n));
const fmtMs = (ms) => (ms >= 1000 ? (ms / 1000).toFixed(2) + " s" : ms >= 10 ? ms.toFixed(0) + " ms" : ms.toFixed(2) + " ms");
const fmtBp = (n) => (n >= 1e6 ? (n / 1e6).toFixed(n >= 1e7 ? 0 : 1) + " Mbp" : n >= 1e3 ? (n / 1e3).toFixed(n >= 1e4 ? 0 : 1) + " kbp" : n + " bp");
function fmtCompact(n) {
  if (n >= 1e9) return (n / 1e9).toFixed(1) + "B";
  if (n >= 1e6) return (n / 1e6).toFixed(1) + "M";
  if (n >= 1e4) return (n / 1e3).toFixed(0) + "K";
  if (n >= 1e3) return (n / 1e3).toFixed(1) + "K";
  return String(Math.round(n));
}
function fmtBytes(b) {
  if (b >= 1 << 30) return (b / (1 << 30)).toFixed(1) + " GB";
  if (b >= 1 << 20) return (b / (1 << 20)).toFixed(1) + " MB";
  if (b >= 1 << 10) return (b / (1 << 10)).toFixed(1) + " KB";
  return b + " B";
}
function fmtUptime(s) {
  if (s < 60) return Math.floor(s) + "s";
  if (s < 3600) return Math.floor(s / 60) + "m " + Math.floor(s % 60) + "s";
  return Math.floor(s / 3600) + "h " + Math.floor((s % 3600) / 60) + "m";
}
const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const cleanSeq = (s) => s.replace(/\s+/g, "").toUpperCase();
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

const store = {
  get(k) { try { return localStorage.getItem(k); } catch { return null; } },
  set(k, v) { try { localStorage.setItem(k, v); } catch { /* storage unavailable */ } },
};

// Worker identity: categorical slots in fixed rank order; past 8 workers the
// colour folds to neutral (the "Rank N" label always carries identity).
const workerColor = (rank) => (rank >= 1 && rank <= 8 ? `var(--s${rank})` : "var(--s-other)");
const workerName = (rank) => (rank === 0 ? "Coordinator" : `Rank ${rank}`);
const workerTag = (rank) => `<span class="wk"><span class="swatch" style="background:${workerColor(rank)}"></span>${workerName(rank)}</span>`;

async function api(path, body) {
  const opts = body === undefined ? {} : { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) };
  let res;
  try {
    res = await fetch(path, opts);
  } catch {
    throw new Error("Gateway not reachable — is web/gateway.py running?");
  }
  const data = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
  return data;
}

// ---------------------------------------------------------------- DNA generators
const BASES = "ACGT";
function randomDna(n) {
  const buf = new Uint8Array(n);
  crypto.getRandomValues(buf);
  let s = "";
  for (let i = 0; i < n; i++) s += BASES[buf[i] & 3];
  return s;
}
// Copy with ~rate substitutions and occasional short indels, so alignments
// against the original are long and interesting.
function mutate(seq, rate = 0.08) {
  let out = "";
  for (let i = 0; i < seq.length; i++) {
    const r = Math.random();
    if (r < rate * 0.7) out += BASES[(BASES.indexOf(seq[i]) + 1 + Math.floor(Math.random() * 3)) % 4];
    else if (r < rate * 0.85) continue; // deletion
    else if (r < rate) out += seq[i] + BASES[Math.floor(Math.random() * 4)]; // insertion
    else out += seq[i];
  }
  return out;
}
// Simulated read chunk: random background with telomeric repeats (TTAGGG)
// mixed in, so the k-mer histogram has a visible signal.
function sequencerChunk(n) {
  let s = "";
  while (s.length < n) s += Math.random() < 0.25 ? "TTAGGG".repeat(1 + Math.floor(Math.random() * 4)) : randomDna(24);
  return s.slice(0, n);
}
function parseSeqList(text) {
  const t = text.trim();
  if (!t) return [];
  if (t.startsWith(">")) {
    const out = [];
    let cur = null;
    for (const line of t.split(/\r?\n/)) {
      if (line.startsWith(">")) { if (cur !== null) out.push(cur); cur = ""; }
      else if (cur !== null) cur += line.trim();
    }
    if (cur !== null) out.push(cur);
    return out.filter(Boolean);
  }
  return t.split(/\r?\n/).map((l) => l.trim()).filter(Boolean);
}

let sampleCache = null;
async function sample() {
  if (!sampleCache) sampleCache = (await api("/api/sample")).records;
  const by = Object.fromEntries(sampleCache.map((r) => [r.name, r.sequence]));
  return { records: sampleCache, by };
}

// ---------------------------------------------------------------- theme
(function initTheme() {
  const saved = store.get("theme");
  if (saved) document.documentElement.dataset.theme = saved;
  $("#theme").addEventListener("click", () => {
    const cur = document.documentElement.dataset.theme
      || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
    const next = cur === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = next;
    store.set("theme", next);
  });
})();

// ---------------------------------------------------------------- tabs
function selectTab(name) {
  $$(".tab").forEach((t) => t.setAttribute("aria-selected", String(t.dataset.tab === name)));
  $$(".panel").forEach((p) => (p.hidden = p.dataset.panel !== name));
  store.set("tab", name);
}
$$(".tab").forEach((t) => t.addEventListener("click", () => selectTab(t.dataset.tab)));
if (store.get("tab") && $(`.tab[data-tab="${store.get("tab")}"]`)) selectTab(store.get("tab"));

// ---------------------------------------------------------------- tooltip
const tip = $("#tooltip");
function showTip(ev, html) {
  tip.innerHTML = html;
  tip.classList.add("show");
  const pad = 14, w = tip.offsetWidth, h = tip.offsetHeight;
  let x = ev.clientX + pad, y = ev.clientY + pad;
  if (x + w > innerWidth - 8) x = ev.clientX - w - pad;
  if (y + h > innerHeight - 8) y = ev.clientY - h - pad;
  tip.style.left = Math.max(8, x) + "px";
  tip.style.top = Math.max(8, y) + "px";
}
const hideTip = () => tip.classList.remove("show");
const ttRows = (title, rows) => `<div class="tt-title">${title}</div>` + rows.map(([k, v]) => `<div class="tt-row"><span>${k}</span><b>${v}</b></div>`).join("");

// ---------------------------------------------------------------- charts (SVG)
const SVGNS = "http://www.w3.org/2000/svg";
function el(tag, attrs = {}, parent) {
  const e = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}
// Round tick step (1/2/5 × 10^n) giving at most ~`intervals` intervals;
// integer data gets integer ticks.
function niceScale(v, intervals = 5, integer = false) {
  v = Math.max(v, integer ? 1 : 1e-9);
  const raw = v / intervals;
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  let step = [1, 2, 5, 10].map((m) => m * p).find((s) => s >= raw);
  if (integer) step = Math.max(1, Math.round(step));
  const n = Math.ceil(v / step - 1e-9);
  return { max: n * step, ticks: Array.from({ length: n + 1 }, (_, i) => i * step) };
}
// Horizontal bar: square at the baseline (x0), 4px rounded data end.
function hbarPath(x0, y, w, h, r = 4) {
  r = Math.min(r, w, h / 2);
  return `M${x0},${y}H${x0 + w - r}Q${x0 + w},${y} ${x0 + w},${y + r}V${y + h - r}Q${x0 + w},${y + h} ${x0 + w - r},${y + h}H${x0}Z`;
}

// data: [{label, value, color?, tip?}] — one series, sorted by caller.
function barChart(container, data, { valueLabel = "Count", labelWidth = 90, mono = true } = {}) {
  container.innerHTML = "";
  const W = Math.max(320, container.clientWidth || 600);
  const rowH = 22, barH = 14, top = 6, bottom = 24, right = 56;
  const H = top + data.length * rowH + bottom;
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": `${valueLabel} bar chart` }, container);
  const x0 = labelWidth, plotW = W - labelWidth - right;
  const { max, ticks } = niceScale(Math.max(...data.map((d) => d.value)), 4, true);
  for (const t of ticks) {
    const x = x0 + (t / max) * plotW;
    el("line", { x1: x, x2: x, y1: top, y2: H - bottom, class: t === 0 ? "baseline" : "grid" }, svg);
    el("text", { x, y: H - 8, "text-anchor": "middle", class: "axis-label" }, svg).textContent = fmtCompact(t);
  }
  data.forEach((d, i) => {
    const y = top + i * rowH + (rowH - barH) / 2;
    const w = Math.max(1, (d.value / max) * plotW);
    const lbl = el("text", { x: x0 - 8, y: y + barH / 2 + 4, "text-anchor": "end" }, svg);
    lbl.textContent = d.label;
    if (mono) lbl.style.fontFamily = "var(--mono)";
    const bar = el("path", { d: hbarPath(x0, y, w, barH), fill: d.color || "var(--accent)", class: "mark" }, svg);
    el("text", { x: x0 + w + 6, y: y + barH / 2 + 4 }, svg).textContent = fmtInt(d.value);
    const hit = el("rect", { x: 0, y: top + i * rowH, width: W, height: rowH, class: "hit" }, svg);
    hit.addEventListener("mousemove", (ev) => { bar.classList.add("hover"); showTip(ev, d.tip || ttRows(esc(d.label), [[valueLabel, fmtInt(d.value)]])); });
    hit.addEventListener("mouseleave", () => { bar.classList.remove("hover"); hideTip(); });
  });
}

// Gantt-style timeline: one lane per rank, one bar per request covering the
// time the worker spent computing it. The reply returns right after the
// worker finishes, so compute ≈ [end - compute_ms, end]; time before that was
// spent queueing for an idle worker (shown in the tooltip).
function timeline(container, jobs, { title = "Request" } = {}) {
  container.innerHTML = "";
  const ranks = [...new Set(jobs.map((j) => j.worker))].sort((a, b) => a - b);
  const W = Math.max(320, container.clientWidth || 600);
  const laneH = 26, barH = 16, left = 96, right = 16, top = 6, bottom = 28;
  const H = top + ranks.length * laneH + bottom;
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Per-worker request timeline" }, container);
  const { max: tMax, ticks } = niceScale(Math.max(...jobs.map((j) => j.end_ms)), 5);
  const plotW = W - left - right;
  const X = (t) => left + (t / tMax) * plotW;
  for (const t of ticks) {
    const x = X(t);
    el("line", { x1: x, x2: x, y1: top, y2: H - bottom, class: t === 0 ? "baseline" : "grid" }, svg);
    el("text", { x, y: H - 10, "text-anchor": "middle", class: "axis-label" }, svg).textContent = t === 0 ? "0" : fmtMs(t);
  }
  ranks.forEach((rank, li) => {
    const y = top + li * laneH;
    if (li > 0) el("line", { x1: left, x2: W - right, y1: y, y2: y, class: "grid" }, svg);
    el("rect", { x: 8, y: y + laneH / 2 - 5, width: 10, height: 10, rx: 3, fill: workerColor(rank) }, svg);
    el("text", { x: 24, y: y + laneH / 2 + 4 }, svg).textContent = workerName(rank);
  });
  for (const j of jobs) {
    const li = ranks.indexOf(j.worker);
    const y = top + li * laneH + (laneH - barH) / 2;
    const cStart = Math.max(j.start_ms, j.end_ms - j.compute_ms);
    const x = X(cStart), w = Math.max(3, X(j.end_ms) - x);
    const bar = el("rect", { x: x + 1, y, width: Math.max(1, w - 2), height: barH, rx: 3, fill: workerColor(j.worker), class: "mark" }, svg);
    const hit = el("rect", { x, y: top + li * laneH, width: w, height: laneH, class: "hit" }, svg);
    hit.addEventListener("mousemove", (ev) => {
      bar.classList.add("hover");
      showTip(ev, ttRows(`${title} #${j.i + 1}`, [
        ["Worker", workerName(j.worker)], ["Sent at", fmtMs(j.start_ms)],
        ["Waited for a worker", fmtMs(Math.max(0, cStart - j.start_ms))],
        ["Compute", fmtMs(j.compute_ms)], ["Finished at", fmtMs(j.end_ms)],
      ]));
    });
    hit.addEventListener("mouseleave", () => { bar.classList.remove("hover"); hideTip(); });
  }
}

function jobsPerWorker(items) {
  const m = new Map();
  for (const it of items) m.set(it.worker, (m.get(it.worker) || 0) + 1);
  return [...m.entries()].sort((a, b) => a[0] - b[0]).map(([rank, n]) => ({
    label: workerName(rank), value: n, color: workerColor(rank),
    tip: ttRows(workerName(rank), [["Jobs", fmtInt(n)]]),
  }));
}

const chips = (pairs) => `<div class="chips">${pairs.map(([k, v]) => `<span class="chip">${k} <b>${v}</b></span>`).join("")}</div>`;
const tiles = (items) => `<div class="tiles" style="margin-top:0">${items.map(([label, value, hint]) => `<div class="tile"><div class="label">${label}</div><div class="value">${value}</div>${hint ? `<div class="hint">${hint}</div>` : ""}</div>`).join("")}</div>`;
const showError = (box, err) => { box.innerHTML = `<div class="error"><b>Request failed.</b> ${esc(err.message)}</div>`; };
function jobChips(r, rtt) {
  return chips([
    ["Computed by", r.cache_hit ? "cache (rank 0)" : workerTag(r.worker)],
    ["Compute", r.cache_hit ? "—" : fmtMs(r.compute_ms)],
    ["Middleware", fmtMs(r.total_ms)],
    ["Round trip", fmtMs(rtt)],
    ["Cache", r.cache_hit ? "hit" : "miss"],
  ]);
}

// ---------------------------------------------------------------- activity log
const logBody = $("#log");
function logRow({ op, input, r, rtt, error }) {
  if (logBody.querySelector("td[colspan]")) logBody.innerHTML = "";
  const tr = document.createElement("tr");
  const time = new Date().toLocaleTimeString();
  if (error) {
    tr.innerHTML = `<td>${time}</td><td>${op}</td><td>${esc(input)}</td><td colspan="5" class="err-text">${esc(error)}</td>`;
  } else {
    tr.innerHTML = `<td>${time}</td><td>${op}</td><td>${esc(input)}</td>
      <td>${r.cache_hit ? '<span class="muted">cache</span>' : workerTag(r.worker)}</td>
      <td>${r.cache_hit ? '<span class="ok-text">hit</span>' : "miss"}</td>
      <td class="num">${r.cache_hit ? "—" : r.compute_ms.toFixed(2)}</td>
      <td class="num">${r.total_ms.toFixed(2)}</td>
      <td class="num">${rtt.toFixed(1)}</td>`;
  }
  logBody.prepend(tr);
  while (logBody.children.length > 200) logBody.lastChild.remove();
}
async function timed(op, input, path, body) {
  const t0 = performance.now();
  try {
    const r = await api(path, body);
    const rtt = performance.now() - t0;
    if (r.worker !== undefined) logRow({ op, input, r, rtt });
    return { r, rtt };
  } catch (e) {
    logRow({ op, input, error: e.message });
    throw e;
  }
}

function busy(form, on) { $$("button", form).forEach((b) => (b.disabled = on)); }
function bindForm(id, handler) {
  const form = $(id);
  form.addEventListener("submit", async (ev) => {
    ev.preventDefault();
    busy(form, true);
    try { await handler(new FormData(form), form); } finally { busy(form, false); }
  });
}

// ---------------------------------------------------------------- alignment
function renderAlignment(a, b) {
  const W = 60;
  let html = "";
  for (let s = 0; s < a.length; s += W) {
    let la = "", mid = "", lb = "";
    for (let i = s; i < Math.min(a.length, s + W); i++) {
      const x = a[i], y = b[i];
      if (x === "-" || y === "-") {
        la += x === "-" ? '<span class="gap">-</span>' : x;
        lb += y === "-" ? '<span class="gap">-</span>' : y;
        mid += " ";
      } else if (x === y) { la += x; lb += y; mid += "|"; }
      else { la += `<span class="mm">${x}</span>`; lb += `<span class="mm">${y}</span>`; mid += "·"; }
    }
    html += `<div class="blk"><span class="pos">${s + 1}</span>${la}\n<span class="pos"></span><span class="mid">${mid}</span>\n<span class="pos"></span>${lb}</div>`;
  }
  return html || '<span class="muted">No local alignment with a positive score.</span>';
}
function alignStats(a, b) {
  let m = 0, mm = 0, g = 0;
  for (let i = 0; i < a.length; i++) {
    if (a[i] === "-" || b[i] === "-") g++; else if (a[i] === b[i]) m++; else mm++;
  }
  return { m, mm, g, len: a.length };
}
bindForm("#f-align", async (fd) => {
  const box = $("#r-align");
  const a = cleanSeq(fd.get("a")), b = cleanSeq(fd.get("b"));
  try {
    const { r, rtt } = await timed("align", `${fmtInt(a.length)} × ${fmtInt(b.length)} bp`, "/api/align", {
      a, b, match: +fd.get("match"), mismatch: +fd.get("mismatch"), gap: +fd.get("gap"),
    });
    const st = alignStats(r.aligned_a, r.aligned_b);
    box.innerHTML = `
      <div class="row" style="align-items:center">
        <div><div class="muted" style="font-size:12px">Score</div><div class="big">${fmtInt(r.score)}</div></div>
        <div style="flex:1">${jobChips(r, rtt)}</div>
      </div>
      ${chips([["Aligned length", fmtInt(st.len)], ["Identity", st.len ? ((100 * st.m) / st.len).toFixed(1) + "%" : "—"], ["Matches", fmtInt(st.m)], ["Mismatches", fmtInt(st.mm)], ["Gaps", fmtInt(st.g)]])}
      <div class="legend"><span><i></i>match ( | )</span><span><i style="background:var(--critical-wash)"></i>mismatch ( · )</span><span><i style="background:var(--surface-2)"></i>gap ( - )</span></div>
      <div class="aln">${renderAlignment(r.aligned_a, r.aligned_b)}</div>`;
  } catch (e) { showError(box, e); }
});

// ---------------------------------------------------------------- pattern
function renderMatches(text, positions, m, limit = 20000) {
  const shown = text.slice(0, limit);
  let html = "", cur = 0, i = 0;
  // Merge overlapping matches into highlighted runs.
  while (i < positions.length && positions[i] < shown.length) {
    let s = positions[i], e = s + m;
    while (i + 1 < positions.length && positions[i + 1] <= e) e = Math.max(e, positions[++i] + m);
    i++;
    e = Math.min(e, shown.length);
    html += esc(shown.slice(cur, s)) + `<mark>${esc(shown.slice(s, e))}</mark>`;
    cur = e;
  }
  html += esc(shown.slice(cur));
  if (text.length > limit) html += `<span class="muted"> … (${fmtBp(text.length - limit)} more not shown)</span>`;
  return html;
}
bindForm("#f-pattern", async (fd) => {
  const box = $("#r-pattern");
  const text = cleanSeq(fd.get("text")), pattern = cleanSeq(fd.get("pattern"));
  try {
    const { r, rtt } = await timed("pattern", `${fmtBp(text.length)} / "${pattern.slice(0, 12)}"`, "/api/pattern", { text, pattern });
    box.innerHTML = `
      <div class="row" style="align-items:center">
        <div><div class="muted" style="font-size:12px">Matches</div><div class="big">${fmtInt(r.count)}</div></div>
        <div style="flex:1">${jobChips(r, rtt)}</div>
      </div>
      <div class="seqview">${renderMatches(text, r.positions, pattern.length)}</div>
      <div><h3>Positions (0-based)${r.truncated ? ` — first ${fmtInt(r.positions.length)}` : ""}</h3>
      <div class="positions">${r.positions.length ? r.positions.join(", ") : "none"}</div></div>`;
  } catch (e) { showError(box, e); }
});

// ---------------------------------------------------------------- k-mer
bindForm("#f-kmer", async (fd) => {
  const box = $("#r-kmer");
  const sequence = cleanSeq(fd.get("sequence")), k = +fd.get("k");
  try {
    const { r, rtt } = await timed("kmer", `${fmtBp(sequence.length)}, k=${k}`, "/api/kmer", { sequence, k, top: 30 });
    box.innerHTML = `
      ${tiles([["Total k-mers", fmtInt(r.total)], ["Distinct", fmtInt(r.distinct), `of ${k <= 15 ? fmtInt(4 ** k) : "4^" + k} possible`], ["Most frequent", r.top[0] ? `<span class="mono" style="font-size:18px">${r.top[0].kmer}</span>` : "—", r.top[0] ? fmtInt(r.top[0].count) + "×" : ""]])}
      ${jobChips(r, rtt)}
      <div><h3>Top ${r.top.length} ${k}-mers</h3><div class="chart" id="kmer-chart"></div></div>`;
    barChart($("#kmer-chart"), r.top.map((t) => ({ label: t.kmer, value: t.count })), { labelWidth: Math.min(260, 24 + k * 8) });
  } catch (e) { showError(box, e); }
});

// ---------------------------------------------------------------- batch
bindForm("#f-batch", async (fd) => {
  const box = $("#r-batch");
  const reference = cleanSeq(fd.get("reference"));
  const sequences = parseSeqList(fd.get("sequences")).map(cleanSeq);
  const t0 = performance.now();
  try {
    const r = await api("/api/batch", { reference, sequences });
    const rtt = performance.now() - t0;
    r.results.forEach((x, i) => logRow({ op: "batch align", input: `#${i + 1}: ${fmtInt(reference.length)} × ${fmtInt(sequences[i].length)} bp`, r: x, rtt }));
    const sum = r.results.reduce((s, x) => s + x.compute_ms, 0);
    const used = new Set(r.results.filter((x) => !x.cache_hit).map((x) => x.worker)).size;
    box.innerHTML = `
      ${tiles([
        ["Alignments", fmtInt(r.results.length)],
        ["Wall time", fmtMs(r.wall_ms), "whole batch, in the middleware"],
        ["Σ compute time", fmtMs(sum), "sum over all jobs"],
        ["Effective parallelism", r.wall_ms > 0 ? (sum / r.wall_ms).toFixed(2) + "×" : "—", `${used} worker${used === 1 ? "" : "s"} used`],
      ])}
      <div class="two">
        <div><h3>Jobs per worker</h3><div class="chart" id="batch-chart"></div></div>
        <div class="table-wrap" style="max-height:320px;overflow:auto"><table>
          <thead><tr><th>#</th><th class="num">Score</th><th>Computed by</th><th class="num">Compute ms</th><th>Aligned (start)</th></tr></thead>
          <tbody>${r.results.map((x, i) => `<tr><td>${i + 1}</td><td class="num">${fmtInt(x.score)}</td><td>${x.cache_hit ? '<span class="muted">cache</span>' : workerTag(x.worker)}</td><td class="num">${x.cache_hit ? "—" : x.compute_ms.toFixed(2)}</td><td class="mono">${esc(x.aligned_a.slice(0, 24))}${x.aligned_a.length > 24 ? "…" : ""}</td></tr>`).join("")}</tbody>
        </table></div>
      </div>`;
    barChart($("#batch-chart"), jobsPerWorker(r.results), { valueLabel: "Jobs", labelWidth: 96, mono: false });
  } catch (e) { logRow({ op: "batch align", input: `${sequences.length} sequences`, error: e.message }); showError(box, e); }
});

// ---------------------------------------------------------------- streaming
// Not bindForm: the submit button doubles as "Stop" and must stay enabled.
let streaming = false;
$("#f-stream").addEventListener("submit", async (ev) => {
  ev.preventDefault();
  const btn = $("#stream-btn");
  if (streaming) { streaming = false; return; }
  const fd = new FormData(ev.target);
  const chunks = +fd.get("chunks"), size = +fd.get("size"), k = +fd.get("k"), delay = +fd.get("delay");
  const box = $("#r-stream");
  const agg = new Map(), perWorker = new Map();
  let done = 0, bases = 0, kmers = 0, latency = 0, failed = null;
  streaming = true;
  btn.textContent = "Stop stream";
  box.innerHTML = `<div id="s-tiles"></div><div class="two"><div><h3 id="s-top-h">Top ${k}-mers so far</h3><div class="chart" id="s-chart"></div></div><div><h3>Chunks per worker</h3><div class="chart" id="s-workers"></div></div></div>`;
  const draw = () => {
    $("#s-tiles").innerHTML = tiles([
      ["Chunks", `${fmtInt(done)} / ${fmtInt(chunks)}`, streaming ? "streaming…" : failed ? "stopped on error" : "finished"],
      ["Bases sequenced", fmtBp(bases)],
      ["K-mers counted", fmtCompact(kmers)],
      ["Distinct k-mers", fmtInt(agg.size), k > 6 ? "from each chunk's top 5,000" : ""],
      ["Avg round trip", done ? fmtMs(latency / done) : "—", "per chunk"],
    ]);
    const top = [...agg.entries()].sort((a, b) => b[1] - a[1]).slice(0, 15);
    barChart($("#s-chart"), top.map(([kmer, n]) => ({ label: kmer, value: n })), { labelWidth: Math.min(260, 24 + k * 8) });
    barChart($("#s-workers"), [...perWorker.entries()].sort((a, b) => a[0] - b[0]).map(([rank, n]) => ({ label: workerName(rank), value: n, color: workerColor(rank) })), { valueLabel: "Chunks", labelWidth: 96, mono: false });
  };
  draw();
  try {
    for (let i = 0; i < chunks && streaming; i++) {
      const seq = sequencerChunk(size);
      const { r, rtt } = await timed("stream kmer", `chunk ${i + 1}: ${fmtBp(size)}, k=${k}`, "/api/kmer", { sequence: seq, k, top: 5000 });
      for (const t of r.top) agg.set(t.kmer, (agg.get(t.kmer) || 0) + t.count);
      const w = r.cache_hit ? 0 : r.worker;
      perWorker.set(w, (perWorker.get(w) || 0) + 1);
      done++; bases += size; kmers += r.total; latency += rtt;
      draw();
      if (delay) await sleep(delay);
    }
  } catch (e) { failed = e; }
  streaming = false;
  btn.textContent = "Start stream";
  draw();
  if (failed) box.insertAdjacentHTML("afterbegin", `<div class="error"><b>Stream stopped.</b> ${esc(failed.message)}</div>`);
});

// ---------------------------------------------------------------- load test
bindForm("#f-load", async (fd) => {
  const box = $("#r-load");
  const body = { op: fd.get("op"), requests: +fd.get("requests"), concurrency: +fd.get("concurrency"), length: +fd.get("length") };
  box.innerHTML = `<div class="muted">Running ${fmtInt(body.requests)} requests with ${body.concurrency} concurrent clients…</div>`;
  try {
    const r = await api("/api/loadtest", body);
    for (const j of r.jobs) logRow({ op: `load ${r.op}`, input: `#${j.i + 1}: ${fmtBp(r.length)}`, r: j, rtt: j.end_ms - j.start_ms });
    const sum = r.jobs.reduce((s, j) => s + j.compute_ms, 0);
    const avg = r.jobs.reduce((s, j) => s + (j.end_ms - j.start_ms), 0) / r.jobs.length;
    box.innerHTML = `
      ${tiles([
        ["Throughput", r.throughput.toFixed(1) + "/s", `${fmtInt(r.requests)} requests`],
        ["Wall time", fmtMs(r.wall_ms)],
        ["Effective parallelism", (sum / r.wall_ms).toFixed(2) + "×", "Σ compute ÷ wall time"],
        ["Avg round trip", fmtMs(avg), `with ${r.concurrency} concurrent clients`],
      ])}
      <div><h3>Timeline — each bar is one request's compute time, on the worker that ran it</h3><div class="chart" id="load-tl"></div></div>
      <div class="two">
        <div><h3>Requests per worker</h3><div class="chart" id="load-dist"></div></div>
        <div class="muted" style="font-size:12.5px;align-self:center">
          With more concurrent clients than workers, requests queue in the coordinator until a worker is idle (hover a bar to see the wait).
          Try concurrency 1 (only one row busy at any moment) versus concurrency ≥ the number of workers (every row busy at once).
        </div>
      </div>`;
    timeline($("#load-tl"), r.jobs);
    barChart($("#load-dist"), jobsPerWorker(r.jobs), { valueLabel: "Requests", labelWidth: 96, mono: false });
  } catch (e) { showError(box, e); }
});

// ---------------------------------------------------------------- sample/random fillers
const fillers = {
  async "align-sample"() {
    const { by, records } = await sample();
    const f = $("#f-align");
    f.a.value = by.seq3 || records[0]?.sequence || "";
    f.b.value = by.ref || records[1]?.sequence || "";
  },
  "align-random"() {
    const a = randomDna(1500);
    const f = $("#f-align");
    f.a.value = a;
    f.b.value = randomDna(200) + mutate(a.slice(300, 1300)) + randomDna(200);
  },
  "pattern-random"() {
    const f = $("#f-pattern");
    f.text.value = randomDna(2_000_000);
    f.pattern.value = "GATTACA";
  },
  "kmer-random"() {
    const f = $("#f-kmer");
    f.sequence.value = randomDna(500_000);
    f.k.value = 6;
  },
  async "batch-sample"() {
    const { by, records } = await sample();
    const f = $("#f-batch");
    f.reference.value = by.ref || "";
    f.sequences.value = records.filter((r) => r.name !== "ref").map((r) => `>${r.name}\n${r.sequence}`).join("\n");
  },
  "batch-random"() {
    const ref = randomDna(1200);
    const f = $("#f-batch");
    f.reference.value = ref;
    f.sequences.value = Array.from({ length: 16 }, (_, i) => `>variant${i + 1}\n${mutate(ref, 0.04 + i * 0.01)}`).join("\n");
  },
};
$$("[data-fill]").forEach((b) => b.addEventListener("click", async () => {
  try { await fillers[b.dataset.fill](); } catch (e) { alert(e.message); }
}));
fillers["batch-sample"]().catch(() => {});

// ---------------------------------------------------------------- live cluster panel
const workerEls = new Map();
function renderWorkers(workers) {
  const grid = $("#workers");
  const key = workers.map((w) => w.rank).join(",");
  if (grid.dataset.key !== key) {
    grid.dataset.key = key;
    grid.innerHTML = "";
    workerEls.clear();
    for (const w of workers) {
      const d = document.createElement("div");
      d.className = "worker";
      d.innerHTML = `<div class="top"><span class="swatch" style="background:${workerColor(w.rank)}"></span><span class="name">${w.rank === 0 ? "Rank 0 (local)" : "Rank " + w.rank}</span><span class="state"></span></div>
        <div class="host"></div>
        <div class="nums"><span>Jobs</span><span class="jobs"></span><span>OpenMP threads</span><span class="omp"></span><span>Compute</span><span class="cms"></span><span>MPI in/out</span><span class="io"></span></div>
        <div class="share" title="Share of all jobs"><i style="background:${workerColor(w.rank)}"></i></div>`;
      grid.appendChild(d);
      workerEls.set(w.rank, d);
    }
  }
  const total = workers.reduce((s, w) => s + w.jobs, 0) || 1;
  for (const w of workers) {
    const d = workerEls.get(w.rank);
    d.classList.toggle("busy", w.active > 0);
    $(".state", d).textContent = w.active > 0 ? `● ${w.current_task}` : "idle";
    $(".host", d).textContent = w.hostname;
    $(".omp", d).textContent = w.omp_threads;
    $(".jobs", d).textContent = fmtInt(w.jobs) + (w.failures ? ` (${w.failures} failed)` : "");
    $(".cms", d).textContent = fmtMs(w.compute_ms);
    $(".io", d).textContent = w.rank === 0 ? "—" : `${fmtBytes(w.bytes_sent)} / ${fmtBytes(w.bytes_received)}`;
    $(".share i", d).style.width = ((100 * w.jobs) / total).toFixed(1) + "%";
  }
}
function setConn(ok, text) {
  const pill = $("#conn");
  $(".dot", pill).className = "dot " + (ok ? "good" : "bad");
  $(".ptext", pill).textContent = text;
  $("#offline").classList.toggle("show", !ok);
}
async function pollStats() {
  try {
    const s = await api("/api/stats");
    const nWorkers = s.workers.filter((w) => w.rank > 0).length;
    setConn(true, `Connected · ${nWorkers ? nWorkers + " worker" + (nWorkers === 1 ? "" : "s") : "local mode"}`);
    renderWorkers(s.workers);
    const lookups = s.cache.hits + s.cache.misses;
    $("#coord-host").textContent = s.coordinator_host;
    $("#c-req").textContent = fmtInt(s.requests);
    $("#c-inflight").textContent = fmtInt(s.in_flight);
    $("#c-hit").textContent = lookups ? ((100 * s.cache.hits) / lookups).toFixed(1) + "%" : "—";
    $("#c-cache").textContent = `${fmtInt(s.cache.entries)} · ${fmtBytes(s.cache.bytes)}`;
    $("#t-workers").textContent = nWorkers || "0";
    $("#t-workers-h").textContent = nWorkers ? `${new Set(s.workers.map((w) => w.hostname)).size} host(s)` : "computing on rank 0";
    $("#t-jobs").textContent = fmtInt(s.workers.reduce((a, w) => a + w.jobs, 0));
    $("#t-bytes").textContent = fmtBytes(s.workers.reduce((a, w) => a + w.bytes_sent + w.bytes_received, 0));
    $("#t-omp").textContent = s.openmp ? "On" : "Off";
    $("#t-omp-h").textContent = s.openmp ? `${Math.max(...s.workers.map((w) => w.omp_threads))} threads per rank` : "built without OpenMP";
    $("#t-thread").textContent = s.mpi_thread_multiple ? "Multiple" : "Serialized";
    $("#t-up").textContent = fmtUptime(s.uptime_s);
    $("#t-target").textContent = "gRPC " + s.target;
    const active = s.in_flight > 0;
    $("#link-grpc").classList.toggle("active", active);
    $("#link-http").classList.toggle("active", active);
    $("#link-mpi").classList.toggle("active", s.workers.some((w) => w.rank > 0 && w.active > 0));
  } catch (e) {
    setConn(false, e.message.startsWith("Gateway") ? "Gateway offline" : "Middleware offline");
  }
}
pollStats();
setInterval(pollStats, 1000);

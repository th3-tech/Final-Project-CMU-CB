/* app.js — rendering + interactivity for the Chou-Fasman / chameleon explorer.
 * Predictions are recomputed live from each sequence via ChouFasman (chou_fasman.js);
 * data.js supplies sequences, deposited truth, per-residue chameleon flags, accuracy. */
(function () {
  "use strict";
  const DATA = PROTEIN_DATA;
  const CF = ChouFasman;
  const byId = {};
  DATA.proteins.forEach(p => (byId[p.id] = p));

  const el = s => document.querySelector(s);
  const state = { protein: null, k: 6, reveal: true, cols: 50 };
  let cur = null; // {p, pred4, pred3, cham}

  const STATE_CLASS = { H: "h", E: "e", T: "t", C: "c" };
  const STATE_NAME = { H: "helix", E: "strand", T: "turn", C: "coil" };

  /* ---------- selection & search ---------- */
  function populateSelect(filter) {
    const sel = el("#protein");
    const q = (filter || "").trim().toUpperCase();
    const matches = DATA.proteins.filter(p =>
      !q || p.id.includes(q) || p.name.toUpperCase().includes(q));
    const shown = matches.slice(0, 500);
    sel.innerHTML = shown.map(p =>
      `<option value="${p.id}">${p.id} — ${escapeHtml(trim(p.name, 46))}</option>`).join("");
    if (matches.length > shown.length) {
      sel.innerHTML += `<option disabled>…${matches.length - shown.length} more, keep typing</option>`;
    }
    if (state.protein && shown.some(p => p.id === state.protein.id)) sel.value = state.protein.id;
  }

  function chooseDefault() {
    // an illustrative default: modest length & accuracy, with visible chameleon bands
    for (const p of DATA.proteins) {
      const chamCount = countOnes(p.cham6);
      if (p.sequence.length <= 220 && p.accuracy >= 0.35 && p.accuracy <= 0.65 && chamCount >= 8)
        return p;
    }
    return DATA.proteins[0];
  }

  /* ---------- core render ---------- */
  function render() {
    const p = state.protein;
    const pred4 = CF.predict(p.sequence);
    const pred3 = pred4.replace(/T/g, "C");
    const cham = p["cham" + state.k];
    cur = { p, pred4, pred3, cham };

    renderTitle(p);
    renderStatbar(p, pred3, cham);
    renderTracks(p, pred4, pred3, cham);
    renderProteinBars(p, pred3, cham);
  }

  function renderTitle(p) {
    el("#protein-title").innerHTML =
      `${p.id} &nbsp;<small>${escapeHtml(p.name)} · ${p.sequence.length} residues</small>`;
  }

  function renderStatbar(p, pred3, cham) {
    const n = p.sequence.length;
    let correct = 0, inN = 0, inErr = 0, outN = 0, outErr = 0;
    for (let i = 0; i < n; i++) {
      const wrong = pred3[i] !== p.true_sse[i];
      if (!wrong) correct++;
      if (cham[i] === "1") { inN++; if (wrong) inErr++; }
      else { outN++; if (wrong) outErr++; }
    }
    const pct = x => (x * 100).toFixed(1) + "%";
    const inRate = inN ? pct(inErr / inN) : "—";
    const outRate = outN ? pct(outErr / outN) : "—";
    el("#statbar").innerHTML = `
      ${stat(pct(correct / n), "prediction accuracy (3-state)")}
      ${stat(`${inN} <small>(${(inN / n * 100).toFixed(0)}%)</small>`,
             `residues inside a k=${state.k} chameleon`)}
      ${stat(inRate, "error rate INSIDE chameleons", true)}
      ${stat(outRate, "error rate OUTSIDE chameleons")}`;
  }
  const stat = (num, lbl, alarm) =>
    `<div class="stat${alarm ? " alarm" : ""}"><div class="num">${num}</div><div class="lbl">${lbl}</div></div>`;

  function renderTracks(p, pred4, pred3, cham) {
    const n = p.sequence.length, cols = state.cols, blocks = [];
    for (let s = 0; s < n; s += cols) {
      let html = '<div class="block">';
      for (let i = s; i < Math.min(s + cols, n); i++) {
        const a = p.sequence[i], pr = pred4[i], tr = p.true_sse[i];
        const mm = pred3[i] !== tr;
        html += `<div class="col" data-i="${i}">`
          + `<div class="cell seq">${a}</div>`
          + `<div class="cell ${STATE_CLASS[pr]}">${pr}</div>`
          + `<div class="cell true ${STATE_CLASS[tr]}">${tr}</div>`
          + `<div class="cell mm${mm ? " on" : ""}"></div>`
          + `<div class="band${cham[i] === "1" ? " on" : ""}"></div>`
          + `</div>`;
      }
      html += "</div>";
      blocks.push(html);
    }
    const tracks = el("#tracks");
    tracks.className = "tracks" + (state.reveal ? "" : " truthhidden");
    tracks.innerHTML = blocks.join("");
  }

  function renderProteinBars(p, pred3, cham) {
    const n = p.sequence.length;
    let inN = 0, inErr = 0, outN = 0, outErr = 0;
    for (let i = 0; i < n; i++) {
      const wrong = pred3[i] !== p.true_sse[i];
      if (cham[i] === "1") { inN++; if (wrong) inErr++; } else { outN++; if (wrong) outErr++; }
    }
    const inRate = inN ? inErr / inN : 0, outRate = outN ? outErr / outN : 0;
    el("#protein-bars").innerHTML =
      barRow("inside chameleon", inRate, inN, "inside") +
      barRow("outside", outRate, outN, "outside");
  }
  function barRow(cap, rate, count, cls) {
    const w = (rate * 100).toFixed(0);
    const val = count ? (rate * 100).toFixed(0) + "%" : "n/a";
    return `<div class="bar-row"><div class="cap">${cap}<br><small>${count} res</small></div>
      <div class="bar-track"><div class="bar-fill ${cls}" style="width:${count ? w : 0}%"></div></div>
      <div class="val">${val}</div></div>`;
  }

  /* ---------- dataset-wide summary (static) ---------- */
  function renderDatasetSummary() {
    const m = DATA.meta;
    el("#ds-n").textContent = m.n_chains.toLocaleString();
    const tiles = [5, 6, 7].map(k => {
      const s = m.k_sweep[String(k)];
      const inW = (s.inside * 100).toFixed(0), outW = (s.outside * 100).toFixed(0);
      const primary = k === m.k_primary ? " primary" : "";
      return `<div class="ktile${primary}">
        <div class="khead"><b>k = ${k}${k === m.k_primary ? " (primary)" : ""}</b>
          <span class="p">${s.residues_inside.toLocaleString()} chameleon residues · p = ${fmtP(s.p_value)}</span></div>
        <div class="ksplit">
          <div class="in" style="flex:${s.inside}">${inW}%</div>
          <div class="out" style="flex:${s.outside}">${outW}%</div>
        </div>
        <div class="sub" style="margin:.5em 0 0">inside <b style="color:var(--alarm)">${inW}%</b>
          vs outside <b style="color:var(--e)">${outW}%</b> error &nbsp;(+${((s.inside - s.outside) * 100).toFixed(1)}pp)</div>
      </div>`;
    }).join("");
    el("#ds-summary").innerHTML = `
      <div class="ds-head">
        <div><div class="num">${(m.overall_accuracy * 100).toFixed(1)}%</div><div class="lbl">overall CF accuracy</div></div>
        <div><div class="num">${m.n_residues.toLocaleString()}</div><div class="lbl">residues scored</div></div>
      </div>
      <div class="ksweep">${tiles}</div>`;
  }

  /* ---------- tooltip ---------- */
  function initTooltip() {
    const tip = el("#tooltip"), tracks = el("#tracks");
    tracks.addEventListener("mousemove", e => {
      const col = e.target.closest(".col");
      if (!col) { tip.hidden = true; return; }
      const i = +col.dataset.i, p = cur.p;
      const aa = p.sequence[i], pr = pred4At(i), tr = p.true_sse[i];
      const row = CF.PROP[aa] || [1, 1, 1];
      const isCham = cur.cham[i] === "1";
      const mism = cur.pred3[i] !== tr;
      tip.innerHTML =
        `<b>${aa}${i + 1}</b> &middot; ${aaName(aa)}<hr>` +
        `<div class="trow">predicted <span class="tag ${STATE_CLASS[pr]}">${pr} ${STATE_NAME[pr]}</span></div>` +
        `<div class="trow">true <span class="tag ${STATE_CLASS[tr]}">${tr} ${STATE_NAME[tr]}</span></div>` +
        `<div class="trow" style="margin-top:3px">${mism ? '<span class="cham">✗ mismatch</span>' : '✓ match'}` +
          `${isCham ? ' <span class="cham">· chameleon</span>' : ""}</div><hr>` +
        `<div class="trow">P&alpha; (helix) <b>${row[0].toFixed(2)}</b></div>` +
        `<div class="trow">P&beta; (strand) <b>${row[1].toFixed(2)}</b></div>`;
      tip.hidden = false;
      const pad = 14, w = tip.offsetWidth, h = tip.offsetHeight;
      let x = e.clientX + pad, y = e.clientY + pad;
      if (x + w > innerWidth) x = e.clientX - w - pad;
      if (y + h > innerHeight) y = e.clientY - h - pad;
      tip.style.left = x + "px"; tip.style.top = y + "px";
    });
    tracks.addEventListener("mouseleave", () => (tip.hidden = true));
  }
  const pred4At = i => cur.pred4[i];

  /* ---------- events ---------- */
  function init() {
    state.protein = chooseDefault();
    computeCols();
    populateSelect("");
    renderDatasetSummary();
    render();
    initTooltip();

    el("#search").addEventListener("input", e => populateSelect(e.target.value));
    el("#protein").addEventListener("change", e => {
      const p = byId[e.target.value];
      if (p) { state.protein = p; render(); }
    });
    el("#kbtns").addEventListener("click", e => {
      const b = e.target.closest("button"); if (!b) return;
      state.k = +b.dataset.k;
      [...el("#kbtns").children].forEach(x => x.classList.toggle("active", x === b));
      render();
    });
    el("#reveal").addEventListener("click", e => {
      state.reveal = !state.reveal;
      e.target.setAttribute("aria-pressed", String(state.reveal));
      e.target.innerHTML = state.reveal ? "Shown — click to hide" : "Hidden — click to reveal";
      el("#tracks").className = "tracks" + (state.reveal ? "" : " truthhidden");
    });
    let t;
    addEventListener("resize", () => { clearTimeout(t); t = setTimeout(() => {
      const before = state.cols; computeCols();
      if (state.cols !== before) render();
    }, 150); });
  }

  function computeCols() {
    const w = el("#tracks").clientWidth || 900;
    state.cols = Math.max(30, Math.min(80, Math.floor(w / 15)));
  }

  /* ---------- helpers ---------- */
  function countOnes(s) { let c = 0; for (let i = 0; i < s.length; i++) if (s[i] === "1") c++; return c; }
  function trim(s, n) { return s.length > n ? s.slice(0, n - 1) + "…" : s; }
  function escapeHtml(s) { return s.replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c])); }
  function fmtP(p) { return p < 1e-4 ? p.toExponential(1) : p.toFixed(4); }
  const AANAMES = { A: "Alanine", R: "Arginine", N: "Asparagine", D: "Aspartate", C: "Cysteine",
    E: "Glutamate", Q: "Glutamine", G: "Glycine", H: "Histidine", I: "Isoleucine", L: "Leucine",
    K: "Lysine", M: "Methionine", F: "Phenylalanine", P: "Proline", S: "Serine", T: "Threonine",
    W: "Tryptophan", Y: "Tyrosine", V: "Valine" };
  function aaName(a) { return AANAMES[a] || a; }

  document.addEventListener("DOMContentLoaded", init);
})();

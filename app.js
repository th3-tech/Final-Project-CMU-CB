/* app.js — rendering + interactivity for the Chou-Fasman / chameleon explorer.
 * Predictions are recomputed live from each sequence via ChouFasman (chou_fasman.js);
 * data.js supplies sequences, deposited truth, per-residue chameleon flags, accuracy. */
(function () {
  "use strict";
  const DATA = PROTEIN_DATA, CF = ChouFasman;
  const byId = {};
  DATA.proteins.forEach(p => (byId[p.id] = p));
  const el = s => document.querySelector(s);
  const PINNED = ["7MQ7", "6PFO", "8AX7", "8F23", "6WBJ", "1WSB", "4ANQ", "3ZYA"];  // featured chains, shown first in search
  const REDUCE = matchMedia("(prefers-reduced-motion: reduce)").matches;

  const state = { protein: null, k: 6, reveal: false, confidence: false, show3d: false, cols: 50, first: true };
  const confOpacity = m => (0.35 + Math.min(Math.abs(m.pa - m.pb), 0.5) / 0.5 * 0.65).toFixed(3);
  let cur = null;
  const lastNum = {};

  // chameleon k-mer sets (shipped in data.js) so we can flag a pasted sequence
  const chamSets = {};
  [5, 6, 7].forEach(k => (chamSets[k] = new Set((DATA.meta.chameleon_kmers || {})[String(k)] || [])));
  function computeCham(seq, k) {
    const set = chamSets[k], n = seq.length, f = new Array(n).fill(false);
    for (let i = 0; i <= n - k; i++)
      if (set.has(seq.slice(i, i + k))) for (let j = i; j < i + k; j++) f[j] = true;
    let s = ""; for (let i = 0; i < n; i++) s += f[i] ? "1" : "0"; return s;
  }
  const VALID_AA = "ACDEFGHIKLMNPQRSTVWY";
  // sperm-whale myoglobin — recognizable, and overlaps the k=5 chameleon catalog
  const SAMPLE = "GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEAEMKASEDLKKHGVTVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG";
  function cleanSequence(raw) {
    const body = raw.split(/\r?\n/).filter(l => !l.trim().startsWith(">")).join("");
    const upper = body.toUpperCase().replace(/[^A-Z]/g, "");
    let seq = "", dropped = 0;
    for (const ch of upper) { if (VALID_AA.includes(ch)) seq += ch; else dropped++; }
    let note = "";
    if (dropped) note = `ignored ${dropped} non-standard character${dropped === 1 ? "" : "s"}. `;
    if (seq.length > 2000) { seq = seq.slice(0, 2000); note += "truncated to 2000 residues."; }
    return { seq, note };
  } // for count-up

  const SCLASS = { H: "h", E: "e", T: "t", C: "c" };
  const SNAME = { H: "helix", E: "strand", T: "turn", C: "coil" };

  /* ---------- combobox ---------- */
  const listEl = () => el("#protein-list");
  let activeIdx = -1, shown = [];

  function openList(filter) {
    const q = (filter || "").trim().toUpperCase();
    const matches = DATA.proteins
      .filter(p => !q || p.id.includes(q) || p.name.toUpperCase().includes(q));
    // always show these featured chains first, in this order (those that match the filter)
    const pinned = [];
    for (const id of PINNED) {
      const idx = matches.findIndex(p => p.id === id);
      if (idx >= 0) pinned.push(matches.splice(idx, 1)[0]);
    }
    shown = pinned.concat(matches).slice(0, 200);
    const L = listEl();
    if (!shown.length) {
      L.innerHTML = `<div class="none">no chains match “${escapeHtml(filter)}”</div>`;
    } else {
      L.innerHTML = shown.map((p, i) =>
        `<div class="opt${state.protein && p.id === state.protein.id ? " sel" : ""}" data-id="${p.id}" data-i="${i}" role="option">
           <span class="oid">${p.id}</span><span class="oname">${escapeHtml(p.name)}</span></div>`).join("");
    }
    L.hidden = false;
    requestAnimationFrame(() => L.classList.add("open"));
    el("#search").setAttribute("aria-expanded", "true");
    activeIdx = -1;
  }
  function closeList() {
    const L = listEl();
    L.classList.remove("open");
    el("#search").setAttribute("aria-expanded", "false");
    setTimeout(() => { if (!L.classList.contains("open")) L.hidden = true; }, 200);
  }
  function pick(id) {
    const p = byId[id];
    if (!p) return;
    state.protein = p;
    el("#search").value = "";
    el("#search").blur();
    closeList();
    render();
  }

  /* ---------- default protein ---------- */
  function chooseDefault() {
    if (byId["7MQ7"]) return byId["7MQ7"];   // pinned default
    // otherwise pick the most compelling+representative demo: decent baseline accuracy, a
    // fair number of chameleon residues, and the largest inside-vs-outside error gap.
    let best = null, bestGap = -1;
    for (const p of DATA.proteins) {
      if (p.sequence.length < 90 || p.sequence.length > 210) continue;
      if (p.accuracy < 0.45 || p.accuracy > 0.64) continue;
      if (countOnes(p.cham6) < 12) continue;
      const pred3 = CF.predict3(p.sequence), cham = p.cham6;
      let inN = 0, inE = 0, outN = 0, outE = 0;
      for (let i = 0; i < p.sequence.length; i++) {
        const w = pred3[i] !== p.true_sse[i];
        if (cham[i] === "1") { inN++; if (w) inE++; } else { outN++; if (w) outE++; }
      }
      if (inN < 12 || outN < 20) continue;
      const gap = inE / inN - outE / outN;
      if (gap > bestGap) { bestGap = gap; best = p; }
    }
    return best || DATA.proteins[0];
  }

  /* ---------- render ---------- */
  function render() {
    const p = state.protein;
    const detailed = CF.predictDetailed(p.sequence);
    const pred4 = detailed.pred, pred3 = detailed.pred;
    const cham = p.custom ? computeCham(p.sequence, state.k) : p["cham" + state.k];
    cur = { p, pred4, pred3, cham, meta: detailed.meta };
    renderTitle(p);
    renderStatbar(p, pred3, cham);
    renderTracks(p, pred4, pred3, cham);
    renderProteinBars(p, pred3, cham);
    // the truth toggle is meaningless without a truth track
    const rev = el("#reveal");
    rev.style.opacity = p.custom ? .4 : 1;
    rev.style.pointerEvents = p.custom ? "none" : "auto";
    state.first = false;
    if (state.show3d) render3DMain();
  }

  function renderTitle(p) {
    const badge = p.custom ? ' <span class="sandbox-badge">sandbox · no truth</span>' : "";
    const nm = p.custom ? "Pasted sequence" : escapeHtml(p.name);
    el("#protein-title").innerHTML = `${p.id} &nbsp;<small>${nm} · ${p.sequence.length} aa</small>${badge}`;
  }

  function tally(p, pred3, cham) {
    const n = p.sequence.length, hasTruth = p.true_sse != null;
    let correct = 0, inN = 0, inErr = 0, outN = 0, outErr = 0;
    for (let i = 0; i < n; i++) {
      const wrong = hasTruth && pred3[i] !== p.true_sse[i];
      if (hasTruth && !wrong) correct++;
      if (cham[i] === "1") { inN++; if (wrong) inErr++; }
      else { outN++; if (wrong) outErr++; }
    }
    return { n, hasTruth, acc: hasTruth ? correct / n : null, inN,
             inRate: hasTruth && inN ? inErr / inN : null,
             outN, outRate: hasTruth && outN ? outErr / outN : null };
  }

  function renderStatbar(p, pred3, cham) {
    const t = tally(p, pred3, cham);
    let cards;
    if (t.hasTruth) {
      cards = [
        { key: "acc", val: t.acc * 100, fmt: v => v.toFixed(1) + "%", sub: "", lbl: "prediction accuracy", cls: "" },
        { key: "cnt", val: t.inN, fmt: v => Math.round(v), sub: ` <small>${(t.inN / t.n * 100).toFixed(0)}% of chain</small>`,
          lbl: `residues in a k=${state.k} chameleon`, cls: "accent" },
        { key: "ein", val: t.inRate == null ? null : t.inRate * 100, fmt: v => v.toFixed(1) + "%",
          sub: "", lbl: "error INSIDE chameleons", cls: "alarm" },
        { key: "eout", val: t.outRate == null ? null : t.outRate * 100, fmt: v => v.toFixed(1) + "%",
          sub: "", lbl: "error OUTSIDE chameleons", cls: "" },
      ];
    } else {
      const cnt = s => { let c = 0; for (const x of pred3) if (x === s) c++; return c; };
      cards = [
        { key: "len", val: t.n, fmt: v => Math.round(v), sub: "", lbl: "residues pasted", cls: "" },
        { key: "ph", val: cnt("H") / t.n * 100, fmt: v => v.toFixed(0) + "%", sub: "", lbl: "predicted helix", cls: "" },
        { key: "pe", val: cnt("E") / t.n * 100, fmt: v => v.toFixed(0) + "%", sub: "", lbl: "predicted strand", cls: "" },
        { key: "cnt", val: t.inN, fmt: v => Math.round(v), sub: ` <small>${(t.inN / t.n * 100).toFixed(0)}%</small>`,
          lbl: `in a known k=${state.k} chameleon`, cls: "alarm" },
      ];
    }
    const bar = el("#statbar");
    if (!bar.children.length) {
      bar.innerHTML = cards.map((c, i) =>
        `<div class="stat" style="--i:${i}"><div class="num"><span class="num-val"></span><span class="num-sub"></span></div>
         <div class="lbl"></div></div>`).join("");
      if (!REDUCE) [...bar.children].forEach(x => x.classList.add("in"));
    }
    cards.forEach((c, i) => {
      const card = bar.children[i];
      card.classList.remove("alarm", "accent");
      if (c.cls) card.classList.add(c.cls);
      card.querySelector(".lbl").innerHTML = c.lbl;
      const numEl = card.querySelector(".num-val"), subEl = card.querySelector(".num-sub");
      subEl.innerHTML = c.val == null ? "" : c.sub;
      if (c.val == null) { numEl.textContent = "—"; lastNum[c.key] = null; return; }
      countUp(numEl, lastNum[c.key], c.val, c.fmt);
      lastNum[c.key] = c.val;
    });
  }

  function countUp(node, from, to, fmt) {
    if (REDUCE || from == null || from === to) { node.textContent = fmt(to); return; }
    const dur = 600, t0 = performance.now();
    function step(now) {
      const k = Math.min(1, (now - t0) / dur);
      const e = 1 - Math.pow(1 - k, 3);
      node.textContent = fmt(from + (to - from) * e);
      if (k < 1) requestAnimationFrame(step);
    }
    requestAnimationFrame(step);
  }

  function renderTracks(p, pred4, pred3, cham) {
    const n = p.sequence.length, cols = state.cols, blocks = [], hasTruth = p.true_sse != null;
    for (let s = 0; s < n; s += cols) {
      let html = '<div class="block">';
      for (let i = s; i < Math.min(s + cols, n); i++) {
        const a = p.sequence[i], pr = pred4[i], tr = hasTruth ? p.true_sse[i] : null;
        const mm = hasTruth && pred3[i] !== tr, isC = cham[i] === "1";
        const predStyle = state.confidence ? ` style="opacity:${confOpacity(cur.meta[i])}"` : "";
        html += `<div class="col${isC ? " chamcol" : ""}" data-i="${i}" style="--col:${i}">`
          + `<div class="cell seq">${a}</div>`
          + `<div class="cell predcell ${SCLASS[pr]}"${predStyle}>${pr}</div>`
          + (hasTruth ? `<div class="cell true ${SCLASS[tr]}">${tr}</div>` : `<div class="cell true unk"></div>`)
          + `<div class="cell mm${mm ? " on" : ""}"></div>`
          + `<div class="band${isC ? " on" + (state.first && !REDUCE ? " draw" : "") : ""}"></div>`
          + `</div>`;
      }
      blocks.push(html + "</div>");
    }
    const tr = el("#tracks");
    tr.className = "tracks" + (state.reveal ? "" : " truthhidden") + (state.first && !REDUCE ? " cascade" : "");
    tr.innerHTML = blocks.join("");
  }

  function renderProteinBars(p, pred3, cham) {
    const t = tally(p, pred3, cham);
    if (!t.hasTruth) {
      const note = t.inN > 0
        ? `A pasted sequence has no solved structure, so there's nothing to score against — no accuracy,
           no error rates. The <b>red band</b> in the tracks marks the <b>${t.inN}</b>
           residue${t.inN === 1 ? "" : "s"} whose k-mer folds <i>both</i> ways elsewhere in the PDB:
           exactly where a purely local method like Chou-Fasman is least able to be right.`
        : `A pasted sequence has no solved structure, so there's no accuracy here. No <b>k=${state.k}</b>
           chameleons were found in it — the k=6 and k=7 catalogs are small and specific to our
           2,001-chain dataset, so novel sequences rarely hit them. Switch to <b>k=5</b> (7,304 motifs)
           to see where this sequence overlaps known chameleons.`;
      el("#protein-bars").innerHTML = `<div class="sandbox-note">${note}</div>`;
      return;
    }
    el("#protein-bars").innerHTML =
      barRow("inside chameleon", t.inRate, t.inN, "inside") +
      barRow("outside", t.outRate, t.outN, "outside");
    // animate widths after paint
    requestAnimationFrame(() => el("#protein-bars").querySelectorAll(".bar-fill")
      .forEach(f => (f.style.width = f.dataset.w + "%")));
  }
  function barRow(cap, rate, count, cls) {
    const w = rate == null ? 0 : (rate * 100).toFixed(0);
    const val = rate == null ? "n/a" : (rate * 100).toFixed(0) + "%";
    return `<div class="bar-row"><div class="cap">${cap}<br><small>${count} res</small></div>
      <div class="bar-track"><div class="bar-fill ${cls}" data-w="${w}"></div></div>
      <div class="val">${val}</div></div>`;
  }

  // Feature D: shade predicted cells by decision margin (cross-fades existing cells)
  function applyConfidence() {
    document.querySelectorAll("#tracks .col").forEach(colEl => {
      const pc = colEl.querySelector(".predcell");
      if (pc) pc.style.opacity = state.confidence ? confOpacity(cur.meta[+colEl.dataset.i]) : "";
    });
  }

  /* ---------- dataset summary (static) ---------- */
  function renderDataset() {
    const m = DATA.meta;
    el("#ds-n").textContent = m.n_chains.toLocaleString();
    const tiles = [5, 6, 7].map(k => {
      const s = m.k_sweep[String(k)];
      const inW = (s.inside * 100).toFixed(0), outW = (s.outside * 100).toFixed(0);
      const prim = k === m.k_primary;
      return `<div class="ktile${prim ? " primary" : ""}">
        <div class="khead"><b>k = ${k}${prim ? " · primary" : ""}</b>
          <span class="p">${s.residues_inside.toLocaleString()} residues · p ${fmtP(s.p_value)}</span></div>
        <div class="ksplit"><div class="in" style="flex:${s.inside}">${inW}%</div>
          <div class="out" style="flex:${s.outside}">${outW}%</div></div>
        <div class="gap">inside <b>${inW}%</b> vs outside ${outW}% &nbsp;(+${((s.inside - s.outside) * 100).toFixed(1)}pp error)</div>
      </div>`;
    }).join("");
    el("#ds-summary").innerHTML = `
      <div class="ds-head">
        <div><div class="num">${(m.overall_accuracy * 100).toFixed(1)}%</div><div class="lbl">overall CF accuracy</div></div>
        <div><div class="num">${m.n_residues.toLocaleString()}</div><div class="lbl">residues scored</div></div>
      </div><div class="ksweep">${tiles}</div>`;
  }

  /* ---------- Feature C: aggregate scatter ---------- */
  function renderScatter(k) {
    const rows = (typeof AGGREGATE !== "undefined" && AGGREGATE[String(k)]) || [];
    const summ = DATA.meta.aggregate_summary[String(k)];
    const W = 520, H = 470, mL = 48, mR = 16, mT = 14, mB = 46;
    const pw = W - mL - mR, ph = H - mT - mB;
    const X = v => mL + v * pw, Y = v => mT + (1 - v) * ph;
    let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="per-chain error scatter">`;
    for (const t of [0, .25, .5, .75, 1]) {
      s += `<line class="sgrid" x1="${X(0)}" y1="${Y(t)}" x2="${X(1)}" y2="${Y(t)}"/>`
        + `<line class="sgrid" x1="${X(t)}" y1="${Y(0)}" x2="${X(t)}" y2="${Y(1)}"/>`
        + `<text class="stick" x="${mL - 8}" y="${Y(t) + 3}" text-anchor="end">${t * 100}</text>`
        + `<text class="stick" x="${X(t)}" y="${Y(0) + 18}" text-anchor="middle">${t * 100}</text>`;
    }
    s += `<line class="sdiag" x1="${X(0)}" y1="${Y(0)}" x2="${X(1)}" y2="${Y(1)}"/>`
      + `<text class="saxis" x="${X(.5)}" y="${H - 8}" text-anchor="middle">error outside chameleons (%)</text>`
      + `<text class="saxis" transform="rotate(-90 13 ${mT + ph / 2})" x="13" y="${mT + ph / 2}" text-anchor="middle">error inside (%)</text>`;
    for (const r of rows) {
      const above = r.inside > r.outside;
      s += `<circle class="sdot ${above ? "above" : "below"}" cx="${X(r.outside).toFixed(1)}" cy="${Y(r.inside).toFixed(1)}" r="3" data-id="${r.id}" data-in="${r.inside}" data-out="${r.outside}"/>`;
    }
    el("#scatter-wrap").innerHTML = s + "</svg>";
    el("#scatter-stat").innerHTML = `<b>${(summ.pct_above_line * 100).toFixed(0)}%</b> of ${rows.length.toLocaleString()}`
      + ` chains sit above the line · mean gap <b>+${(summ.mean_gap * 100).toFixed(0)} pts</b>`;
  }

  function initScatter() {
    const wrap = el("#scatter-wrap"), tip = el("#scatter-tip");
    wrap.addEventListener("mousemove", e => {
      const dot = e.target.closest(".sdot");
      if (!dot) { tip.classList.remove("show"); tip.hidden = true; return; }
      tip.innerHTML = `<span class="sid">${dot.dataset.id}</span> · inside `
        + `${(+dot.dataset.in * 100).toFixed(0)}% · outside ${(+dot.dataset.out * 100).toFixed(0)}%`;
      tip.hidden = false; tip.classList.add("show");
      tip.style.left = (e.clientX + 14) + "px"; tip.style.top = (e.clientY + 14) + "px";
    });
    wrap.addEventListener("mouseleave", () => { tip.classList.remove("show"); tip.hidden = true; });
    wrap.addEventListener("click", e => { const d = e.target.closest(".sdot"); if (d) pick(d.dataset.id); });
  }

  /* ---------- tooltip ---------- */
  function initTooltip() {
    const tip = el("#tooltip"), tracks = el("#tracks");
    const chip = el("#compare-chip");
    let lastInspI = -1;
    tracks.addEventListener("mousemove", e => {
      const col = e.target.closest(".col");
      if (!col) { tip.classList.remove("show"); tip.hidden = true; chip.classList.remove("show"); chip.hidden = true; return; }
      const i = +col.dataset.i, p = cur.p;
      if (cur.cham[i] === "1") {
        chip.hidden = false; chip.classList.add("show");
        chip.style.left = (e.clientX + 12) + "px"; chip.style.top = (e.clientY - 30) + "px";
      } else { chip.classList.remove("show"); chip.hidden = true; }
      const isC = cur.cham[i] === "1";
      if (i !== lastInspI) {        // rebuild only when the residue changes (calm, no strobe)
        lastInspI = i;
        const aa = p.sequence[i], pr = cur.pred4[i], m = cur.meta[i];
        const hasTruth = p.true_sse != null, tr = hasTruth ? p.true_sse[i] : null;
        const wH = m.pa >= m.pb;
        const wa = Math.min(100, m.pa / 2 * 100), wb = Math.min(100, m.pb / 2 * 100);
        const diff = Math.abs(m.pa - m.pb).toFixed(2);
        const stateChip = m.regionAccepted ? '<span class="ichip ok">region accepted</span>'
          : m.reverted ? '<span class="ichip rev">reverted to coil</span>'
          : '<span class="ichip dim">no nucleation</span>';
        const verdict = hasTruth
          ? `predicted <b class="s-${pr}">${SNAME[pr]}</b> · true <b class="s-${tr}">${SNAME[tr]}</b> `
            + `· <span class="${pr === tr ? "ok" : "bad"}">${pr === tr ? "✓" : "✗"}</span>`
          : `predicted <b class="s-${pr}">${SNAME[pr]}</b> · <span class="dimv">no truth</span>`;
        tip.innerHTML =
          `<div class="insp-head"><span class="ia">${aa}</span> · position ${i + 1}`
            + `${isC ? ' <span class="cham">⤢ click to compare</span>' : ""}</div>`
          + `<div class="insp-sub">${aaName(aa)}</div>`
          + `<div class="tug">`
            + `<div class="tug-row"><span class="tl">Pα</span><div class="tbar"><div class="tfill h${wH ? " win" : ""}" style="width:${wa}%"></div></div><span class="tv">${m.pa.toFixed(2)}</span></div>`
            + `<div class="tug-row"><span class="tl">Pβ</span><div class="tbar"><div class="tfill e${!wH ? " win" : ""}" style="width:${wb}%"></div></div><span class="tv">${m.pb.toFixed(2)}</span></div>`
          + `</div><div class="tug-cap">${wH ? "helix" : "strand"} favoured by ${diff}</div>`
          + `<div class="chips"><span class="ichip${m.inNucleus ? " lit" : ""}">nucleation window</span>${stateChip}</div>`
          + `<div class="verdict">${verdict}</div>`;
      }
      tip.hidden = false; tip.classList.add("show");
      const pad = 15, w = tip.offsetWidth, h = tip.offsetHeight;
      let x = e.clientX + pad, y = e.clientY + pad;
      if (x + w > innerWidth) x = e.clientX - w - pad;
      if (y + h > innerHeight) y = e.clientY - h - pad;
      tip.style.left = x + "px"; tip.style.top = y + "px";
    });
    tracks.addEventListener("mouseleave", () => {
      tip.classList.remove("show"); tip.hidden = true;
      chip.classList.remove("show"); chip.hidden = true;
      lastInspI = -1;
    });
    tracks.addEventListener("click", e => {
      const col = e.target.closest(".col"); if (!col) return;
      const i = +col.dataset.i;
      if (state.show3d && mainViewer) Structure3D.pulse(mainViewer, i + 1);
      if (cur.cham[i] === "1") openSpotlightAt(i);
    });
  }

  /* ---------- Feature A: chameleon spotlight ---------- */
  const SNAMELONG = { H: "helix", E: "strand", C: "coil" };
  let spot = null;

  function domLabel(str, pos, k) {
    const cnt = { H: 0, E: 0, C: 0 };
    for (let i = pos; i < pos + k && i < str.length; i++) cnt[str[i] === "H" ? "H" : str[i] === "E" ? "E" : "C"]++;
    let best = "H";
    for (const s of ["E", "C"]) if (cnt[s] > cnt[best]) best = s;
    return best;
  }
  function strictLabel(str, pos, k) {
    let h = 0, e = 0;
    for (let i = pos; i < pos + k && i < str.length; i++) { if (str[i] === "H") h++; else if (str[i] === "E") e++; }
    if (h > k / 2) return "H";
    if (e > k / 2) return "E";
    return "C";
  }
  const corePred = (occ, k) => domLabel(CF.predict(byId[occ.id].sequence), occ.pos, k);

  function openSpotlightAt(i) {
    const k = state.k, idx = (CHAMELEON_INDEX[String(k)]) || {}, seq = cur.p.sequence;
    let kmer = null, kpos = -1;
    for (let j = Math.max(0, i - k + 1); j <= Math.min(i, seq.length - k); j++) {
      if (idx[seq.slice(j, j + k)]) { kmer = seq.slice(j, j + k); kpos = j; break; }
    }
    if (!kmer) return;
    const entry = idx[kmer];
    const H = entry.H.map(([id, pos]) => ({ id, pos }));
    const E = entry.E.map(([id, pos]) => ({ id, pos }));
    // include the clicked protein's own occurrence if it's a clean H or E stretch
    if (!cur.p.custom) {
      const cl = strictLabel(cur.p.true_sse, kpos, k);
      if (cl === "H" && !H.some(o => o.id === cur.p.id && o.pos === kpos)) H.unshift({ id: cur.p.id, pos: kpos });
      if (cl === "E" && !E.some(o => o.id === cur.p.id && o.pos === kpos)) E.unshift({ id: cur.p.id, pos: kpos });
    }
    if (!H.length || !E.length) return;
    // default to a pair where CF's core guess matches (so "same guess" is literally true)
    const pairCount = Math.max(H.length, E.length);
    let def = 0;
    for (let p = 0; p < pairCount; p++)
      if (corePred(H[p % H.length], k) === corePred(E[p % E.length], k)) { def = p; break; }
    spot = { kmer, k, H, E, pairIndex: def };
    el("#spot-scrim").hidden = false;
    const s = el("#spotlight"); s.hidden = false; s.classList.remove("closing");
    renderSpotlight(true);
  }

  function renderSpotlight(firstOpen) {
    const { kmer, k, H, E, pairIndex } = spot;
    const hOcc = H[pairIndex % H.length], eOcc = E[pairIndex % E.length];
    el("#spot-kmer").innerHTML = [...kmer].map(c => `<span class="ktile">${c}</span>`).join("");
    const same = corePred(hOcc, k) === corePred(eOcc, k);
    el("#spot-verdict").innerHTML = same
      ? `Same sequence. Same Chou-Fasman guess (<b>${SNAMELONG[corePred(hOcc, k)]}</b>). <span class="c-alarm">Opposite reality.</span>`
      : `Same sequence, same predictor. <span class="c-alarm">Opposite reality.</span>`;
    el("#spot-panels").innerHTML = renderPanel(hOcc, "H", k, firstOpen) + renderPanel(eOcc, "E", k, firstOpen);
    renderSpot3D(hOcc, k, "H");
    renderSpot3D(eOcc, k, "E");
    const pairCount = Math.max(H.length, E.length);
    const cyc = el("#spot-cycle");
    if (pairCount > 1) {
      cyc.innerHTML = `<button id="spot-prev" aria-label="previous">‹</button>
        <span>example ${pairIndex + 1} of ${pairCount}</span><button id="spot-next" aria-label="next">›</button>`;
      el("#spot-prev").onclick = () => { spot.pairIndex = (spot.pairIndex - 1 + pairCount) % pairCount; renderSpotlight(false); };
      el("#spot-next").onclick = () => { spot.pairIndex = (spot.pairIndex + 1) % pairCount; renderSpotlight(false); };
    } else cyc.innerHTML = "";
  }

  function renderPanel(occ, foldLabel, k, firstOpen) {
    const p = byId[occ.id], seq = p.sequence, truth = p.true_sse, pred = CF.predict(seq);
    const flank = 4, start = Math.max(0, occ.pos - flank), end = Math.min(seq.length, occ.pos + k + flank);
    const stagger = firstOpen && !REDUCE;
    let cols = "";
    for (let i = start; i < end; i++) {
      const core = i >= occ.pos && i < occ.pos + k;
      cols += `<div class="scol${core ? " core" : ""}">`
        + `<div class="scell seq">${seq[i]}</div>`
        + `<div class="scell pred ${SCLASS[pred[i]]}">${pred[i]}</div>`
        + `<div class="scell tru ${SCLASS[truth[i]]}${stagger ? " pending" : ""}">${truth[i]}</div>`
        + `</div>`;
    }
    const cp = domLabel(pred, occ.pos, k), right = cp === foldLabel;
    const cap = `Chou-Fasman predicts <b>${SNAMELONG[cp]}</b> here — and it's `
      + `<span class="${right ? "ok" : "wrong"}">${right ? "right" : "wrong"}</span>`;
    const pill = foldLabel === "H" ? '<span class="spill h">folds as HELIX</span>' : '<span class="spill e">folds as STRAND</span>';
    return `<div class="spot-panel">
      <div class="spot-head"><span class="spid">${occ.id}</span><span class="spname">${escapeHtml(p.name)}</span>${pill}</div>
      <div class="spot-strip">${cols}</div><div class="spot-cap">${cap}</div>
      <div class="spot-3d" data-side="${foldLabel}"></div></div>`;
  }

  function closeSpotlight() {
    const s = el("#spotlight");
    if (s.hidden) return;
    s.classList.add("closing");
    el("#spot-scrim").hidden = true;
    setTimeout(() => { s.hidden = true; s.classList.remove("closing"); }, REDUCE ? 0 : 180);
  }

  /* ---------- Feature E: 3D structure ---------- */
  let _3dmolLoad = null, mainViewer = null;
  function ensure3Dmol() {
    if (typeof $3Dmol !== "undefined") return Promise.resolve();
    if (_3dmolLoad) return _3dmolLoad;
    _3dmolLoad = new Promise((res, rej) => {
      const s = document.createElement("script");
      s.src = "lib/3Dmol-min.js"; s.onload = res; s.onerror = rej;
      document.head.appendChild(s);
    });
    return _3dmolLoad;
  }
  function render3DMain() {
    const cont = el("#s3d-viewer"), p = state.protein;
    if (p.custom) { cont.innerHTML = '<div class="s3d-msg">not available for pasted sequences</div>'; mainViewer = null; return; }
    cont.innerHTML = '<div class="s3d-msg">loading 3Dmol…</div>';
    ensure3Dmol().then(() => Structure3D.render(cont, {
      id: p.id, sequence: p.sequence, trueSse: p.true_sse, pred: cur.pred4,
      isCham: i => cur.cham[i] === "1"
    })).then(v => { mainViewer = v; })
      .catch(() => { cont.innerHTML = '<div class="s3d-msg">3D unavailable</div>'; });
  }
  function renderSpot3D(occ, k, side) {
    const cont = document.querySelector('#spot-panels .spot-3d[data-side="' + side + '"]');
    if (!cont) return;
    const p = byId[occ.id];
    ensure3Dmol().then(() => Structure3D.render(cont, {
      id: p.id, sequence: p.sequence, trueSse: p.true_sse, pred: CF.predict(p.sequence),
      isCham: i => i >= occ.pos && i < occ.pos + k, bg: "#0e1420"
    }));
  }

  /* ---------- popover ---------- */
  function togglePop(show) {
    el("#info-pop").hidden = !show;
    el("#pop-scrim").hidden = !show;
  }

  /* ---------- custom-sequence sandbox ---------- */
  function openSeq() {
    el("#seq-modal").hidden = false;
    el("#seq-scrim").hidden = false;
    setTimeout(() => el("#seq-input").focus(), 30);
  }
  function closeSeq() {
    el("#seq-modal").hidden = true;
    el("#seq-scrim").hidden = true;
  }
  function analyzeCustom() {
    const { seq, note } = cleanSequence(el("#seq-input").value);
    const warn = el("#seq-warn");
    if (seq.length < 7) {
      warn.className = "seq-warn err";
      warn.textContent = seq.length ? "Need at least 7 residues to run k-mer analysis." : "Paste a sequence first.";
      return;
    }
    state.protein = { id: "CUSTOM", name: "Pasted sequence", sequence: seq, true_sse: null, custom: true };
    state.first = true;              // let the flourishes play for the new sequence
    closeSeq();
    computeCols();
    render();
    window.scrollTo({ top: 0, behavior: REDUCE ? "auto" : "smooth" });
  }

  /* ---------- events ---------- */
  function init() {
    state.protein = chooseDefault();
    computeCols();
    renderDataset();
    renderScatter(state.k);
    render();
    initTooltip();
    initScatter();
    // sync pill to default k
    positionPill();

    const search = el("#search");
    search.addEventListener("focus", () => openList(search.value));
    search.addEventListener("input", () => openList(search.value));
    search.addEventListener("blur", () => setTimeout(closeList, 120));
    search.addEventListener("keydown", e => {
      const opts = [...listEl().querySelectorAll(".opt")];
      if (e.key === "ArrowDown" || e.key === "ArrowUp") {
        e.preventDefault(); if (!opts.length) return;
        activeIdx = (activeIdx + (e.key === "ArrowDown" ? 1 : -1) + opts.length) % opts.length;
        opts.forEach((o, i) => o.classList.toggle("active", i === activeIdx));
        opts[activeIdx].scrollIntoView({ block: "nearest" });
      } else if (e.key === "Enter" && activeIdx >= 0) {
        pick(opts[activeIdx].dataset.id);
      } else if (e.key === "Escape") { search.blur(); closeList(); }
    });
    listEl().addEventListener("mousedown", e => {
      const o = e.target.closest(".opt"); if (o) { e.preventDefault(); pick(o.dataset.id); }
    });

    el("#kmer").addEventListener("click", e => {
      const b = e.target.closest("button"); if (!b) return;
      state.k = +b.dataset.k;
      [...el("#kmer").querySelectorAll("button")].forEach(x =>
        x.setAttribute("aria-selected", x === b ? "true" : "false"));
      positionPill();
      render();
      renderScatter(state.k);
    });

    el("#reveal").addEventListener("click", e => {
      state.reveal = !state.reveal;
      e.currentTarget.setAttribute("aria-pressed", String(state.reveal));
      el("#tracks").classList.toggle("truthhidden", !state.reveal);
    });

    el("#confidence").addEventListener("click", e => {
      state.confidence = !state.confidence;
      e.currentTarget.setAttribute("aria-pressed", String(state.confidence));
      applyConfidence();
    });

    el("#s3d-toggle").addEventListener("click", e => {
      state.show3d = !state.show3d;
      el("#s3d-body").hidden = !state.show3d;
      e.currentTarget.classList.toggle("open", state.show3d);
      e.currentTarget.innerHTML = state.show3d ? "<span>▸</span> hide 3D" : "<span>▸</span> show 3D";
      if (state.show3d) render3DMain();
    });

    el("#info").addEventListener("click", () => togglePop(true));
    el("#pop-close").addEventListener("click", () => togglePop(false));
    el("#pop-scrim").addEventListener("click", () => togglePop(false));

    el("#custom-btn").addEventListener("click", openSeq);
    el("#seq-close").addEventListener("click", closeSeq);
    el("#seq-scrim").addEventListener("click", closeSeq);
    el("#seq-go").addEventListener("click", analyzeCustom);
    el("#seq-sample").addEventListener("click", () => {
      el("#seq-input").value = SAMPLE;
      el("#seq-warn").textContent = "";
    });
    el("#seq-input").addEventListener("input", () => {
      const { note } = cleanSequence(el("#seq-input").value);
      const w = el("#seq-warn"); w.className = "seq-warn"; w.textContent = note;
    });

    el("#spot-close").addEventListener("click", closeSpotlight);
    el("#spot-scrim").addEventListener("click", closeSpotlight);
    addEventListener("keydown", e => { if (e.key === "Escape") { togglePop(false); closeSeq(); closeSpotlight(); } });

    addEventListener("scroll", () => el("#appbar").classList.toggle("scrolled", scrollY > 4), { passive: true });

    let t;
    addEventListener("resize", () => { clearTimeout(t); t = setTimeout(() => {
      const before = state.cols; computeCols();
      if (state.cols !== before) render();
    }, 150); });
  }

  function positionPill() {
    const kmer = el("#kmer");
    kmer.style.setProperty("--active", { 5: 0, 6: 1, 7: 2 }[state.k]);
  }
  function computeCols() {
    const w = el("#tracks").clientWidth || 900;
    state.cols = Math.max(30, Math.min(80, Math.floor(w / 15)));
  }

  /* ---------- helpers ---------- */
  function countOnes(s) { let c = 0; for (let i = 0; i < s.length; i++) if (s[i] === "1") c++; return c; }
  function escapeHtml(s) { return String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c])); }
  function fmtP(p) { return p < 1e-4 ? p.toExponential(1) : p.toFixed(4); }
  const AAN = { A: "Alanine", R: "Arginine", N: "Asparagine", D: "Aspartate", C: "Cysteine",
    E: "Glutamate", Q: "Glutamine", G: "Glycine", H: "Histidine", I: "Isoleucine", L: "Leucine",
    K: "Lysine", M: "Methionine", F: "Phenylalanine", P: "Proline", S: "Serine", T: "Threonine",
    W: "Tryptophan", Y: "Tyrosine", V: "Valine" };
  function aaName(a) { return AAN[a] || a; }

  document.addEventListener("DOMContentLoaded", init);
})();

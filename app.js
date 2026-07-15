/* app.js — rendering + interactivity for the Chou-Fasman / chameleon explorer.
 * Predictions are recomputed live from each sequence via ChouFasman (chou_fasman.js);
 * data.js supplies sequences, deposited truth, per-residue chameleon flags, accuracy. */
(function () {
  "use strict";
  const DATA = PROTEIN_DATA, CF = ChouFasman;
  const byId = {};
  DATA.proteins.forEach(p => (byId[p.id] = p));
  const el = s => document.querySelector(s);
  const REDUCE = matchMedia("(prefers-reduced-motion: reduce)").matches;

  const state = { protein: null, k: 6, reveal: true, cols: 50, first: true };
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
    shown = DATA.proteins
      .filter(p => !q || p.id.includes(q) || p.name.toUpperCase().includes(q))
      .slice(0, 200);
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
    // pick the most compelling+representative demo: decent baseline accuracy, a fair
    // number of chameleon residues, and the largest inside-vs-outside error gap.
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
    const pred4 = CF.predict(p.sequence);
    const pred3 = pred4.replace(/T/g, "C");
    const cham = p.custom ? computeCham(p.sequence, state.k) : p["cham" + state.k];
    cur = { p, pred4, pred3, cham };
    renderTitle(p);
    renderStatbar(p, pred3, cham);
    renderTracks(p, pred4, pred3, cham);
    renderProteinBars(p, pred3, cham);
    // the truth toggle is meaningless without a truth track
    const rev = el("#reveal");
    const hasTruth = p.true_sse != null;
    rev.style.opacity = p.custom && !hasTruth ? .4 : 1;
    rev.style.pointerEvents = p.custom && !hasTruth ? "none" : "auto";
    state.first = false;
  }

  function renderTitle(p) {
    const hasTruth = p.true_sse != null;
    const badge = p.custom && !hasTruth ? ' <span class="sandbox-badge">sandbox · no truth</span>' : "";
    const nm = p.custom && !hasTruth ? "Pasted sequence" : escapeHtml(p.name);
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
        html += `<div class="col" data-i="${i}" style="--col:${i}">`
          + `<div class="cell seq">${a}</div>`
          + `<div class="cell ${SCLASS[pr]}">${pr}</div>`
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

  /* ---------- tooltip ---------- */
  function initTooltip() {
    const tip = el("#tooltip"), tracks = el("#tracks");
    tracks.addEventListener("mousemove", e => {
      const col = e.target.closest(".col");
      if (!col) { tip.classList.remove("show"); tip.hidden = true; return; }
      const i = +col.dataset.i, p = cur.p;
      const aa = p.sequence[i], pr = cur.pred4[i], tr = p.true_sse[i];
      const row = CF.PROP[aa] || [1, 1];
      const mism = cur.pred3[i] !== tr, isC = cur.cham[i] === "1";
      tip.innerHTML =
        `<div class="thead">${aa}${i + 1}</div><div class="tsub">${aaName(aa)}</div>` +
        `<div class="trow"><span>predicted</span><span class="tag ${SCLASS[pr]}">${pr} ${SNAME[pr]}</span></div>` +
        `<div class="trow"><span>true</span><span class="tag ${SCLASS[tr]}">${tr} ${SNAME[tr]}</span></div>` +
        `<div class="trow"><span>${mism ? '<span class="cham">✗ mismatch</span>' : "✓ match"}</span>` +
          `<span>${isC ? '<span class="cham">chameleon</span>' : ""}</span></div><hr>` +
        `<div class="trow"><span>Pα helix</span><span>${row[0].toFixed(2)}</span></div>` +
        `<div class="trow"><span>Pβ strand</span><span>${row[1].toFixed(2)}</span></div>`;
      tip.hidden = false; tip.classList.add("show");
      const pad = 15, w = tip.offsetWidth, h = tip.offsetHeight;
      let x = e.clientX + pad, y = e.clientY + pad;
      if (x + w > innerWidth) x = e.clientX - w - pad;
      if (y + h > innerHeight) y = e.clientY - h - pad;
      tip.style.left = x + "px"; tip.style.top = y + "px";
    });
    tracks.addEventListener("mouseleave", () => { tip.classList.remove("show"); tip.hidden = true; });
  }

  /* ---------- popover ---------- */
  function togglePop(show) {
    el("#info-pop").hidden = !show;
    el("#pop-scrim").hidden = !show;
  }

  /* ---------- custom-sequence sandbox ---------- */
  function threeLetterToOne(letter) {
    const map = {
      ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C", GLN: "Q", GLU: "E",
      GLY: "G", HIS: "H", ILE: "I", LEU: "L", LYS: "K", MET: "M", PHE: "F",
      PRO: "P", SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V"
    };
    return map[letter] || "X";
  }
  function parseCifLoop(lines, fieldSuffixes) {
    const rows = [];
    for (let i = 0; i < lines.length; i++) {
      if (lines[i].trim() !== "loop_") continue;
      const fields = [];
      i++;
      while (i < lines.length && lines[i].trim().startsWith("_")) {
        fields.push(lines[i].trim());
        i++;
      }
      const indexes = fieldSuffixes.map(sfx => fields.findIndex(f => f.endsWith(sfx)));
      if (indexes.some(idx => idx < 0)) { continue; }
      while (i < lines.length) {
        const rowLine = lines[i].trim();
        if (!rowLine || rowLine.startsWith("#") || rowLine.startsWith("_") || rowLine.startsWith("loop_") || rowLine.startsWith("save_") || rowLine.startsWith("stop_") || rowLine.startsWith("data_")) break;
        const vals = rowLine.split(/\s+/).filter(Boolean);
        const row = {};
        indexes.forEach((idx, j) => { if (idx >= 0 && vals[idx]) row[fieldSuffixes[j]] = vals[idx]; });
        if (Object.keys(row).length) rows.push(row);
        i++;
      }
    }
    return rows;
  }
  function parsePdbStructure(text) {
    const lines = text.split(/\r?\n/);
    const seqRows = parseCifLoop(lines, ["asym_id", "entity_id", "mon_id", "seq_id"]);
    const helixRows = parseCifLoop(lines, ["beg_label_asym_id", "end_label_asym_id", "beg_label_seq_id", "end_label_seq_id"]);
    const sheetRows = parseCifLoop(lines, ["beg_label_asym_id", "end_label_asym_id", "beg_label_seq_id", "end_label_seq_id"]);
    const byChain = {};
    for (const row of seqRows) {
      const chain = row.asym_id;
      const mon = threeLetterToOne(row.mon_id);
      if (!chain || mon === "X") continue;
      if (!byChain[chain]) byChain[chain] = [];
      byChain[chain].push(mon);
    }
    const chains = Object.entries(byChain).filter(([, seq]) => seq.length >= 8);
    if (!chains.length) return null;
    chains.sort((a, b) => b[1].length - a[1].length);
    const [chain, seq] = chains[0];
    const out = Array(seq.length).fill("C");
    const applyRanges = (rows, token) => {
      for (const row of rows) {
        if (row.beg_label_asym_id !== chain && row.end_label_asym_id !== chain) continue;
        const start = Number(row.beg_label_seq_id);
        const end = Number(row.end_label_seq_id);
        if (!Number.isFinite(start) || !Number.isFinite(end)) continue;
        const lo = Math.min(start, end), hi = Math.max(start, end);
        for (let i = lo - 1; i < Math.min(hi, seq.length); i++) {
          if (i >= 0 && i < out.length) out[i] = token;
        }
      }
    };
    applyRanges(helixRows, "H");
    applyRanges(sheetRows, "E");
    return { sequence: seq.join(""), true_sse: out.join("") };
  }
  async function fetchText(url) {
    const res = await fetch(url, { cache: "no-store" });
    if (!res.ok) throw new Error(`Unable to fetch ${url}`);
    return res.text();
  }
  async function fetchPdbFasta(code) {
    const directUrl = `https://www.rcsb.org/fasta/entry/${code}`;
    const proxyUrl = `https://r.jina.ai/http://www.rcsb.org/fasta/entry/${code}`;
    try { return await fetchText(directUrl); } catch (err) { return await fetchText(proxyUrl); }
  }
  async function fetchPdbMmCif(code) {
    const directUrl = `https://files.rcsb.org/download/${code}.cif`;
    const proxyUrl = `https://r.jina.ai/http://files.rcsb.org/download/${code}.cif`;
    try { return await fetchText(directUrl); } catch (err) { return await fetchText(proxyUrl); }
  }
  function openPdb() {
    el("#pdb-modal").hidden = false;
    el("#pdb-scrim").hidden = false;
    setTimeout(() => el("#pdb-input").focus(), 30);
  }
  function openSeq() {
    el("#seq-modal").hidden = false;
    el("#seq-scrim").hidden = false;
    setTimeout(() => el("#seq-input").focus(), 30);
  }
  function closePdb() {
    el("#pdb-modal").hidden = true;
    el("#pdb-scrim").hidden = true;
  }
  function closeSeq() {
    el("#seq-modal").hidden = true;
    el("#seq-scrim").hidden = true;
  }
  async function fetchPdbFasta(code) {
    const directUrl = `https://www.rcsb.org/fasta/entry/${code}`;
    const proxyUrl = `https://r.jina.ai/http://www.rcsb.org/fasta/entry/${code}`;
    try {
      return await fetchText(directUrl, code);
    } catch (err) {
      return await fetchText(proxyUrl, code);
    }
  }
  async function fetchText(url, code) {
    const res = await fetch(url, { cache: "no-store" });
    if (!res.ok) throw new Error(`PDB ${code} not found`);
    return res.text();
  }
  async function importPdb() {
    const code = el("#pdb-input").value.trim().toUpperCase();
    const warn = el("#pdb-warn");
    warn.className = "seq-warn";
    warn.textContent = "";
    if (!/^[A-Z0-9]{4}$/.test(code)) {
      warn.className = "seq-warn err";
      warn.textContent = "Enter a valid 4-character PDB code.";
      return;
    }
    try {
      let seq = "";
      let trueSse = null;
      try {
        const cifText = await fetchPdbMmCif(code);
        const parsed = parsePdbStructure(cifText);
        if (parsed && parsed.sequence) {
          seq = parsed.sequence;
          trueSse = parsed.true_sse;
        }
      } catch (cifErr) {
        // fall back to FASTA if mmCIF parsing is unavailable
      }
      if (!seq) {
        let raw = await fetchPdbFasta(code);
        const idx = raw.indexOf(">");
        if (idx < 0) throw new Error("No FASTA data returned.");
        raw = raw.slice(idx);
        const lines = raw.split(/\r?\n/);
        const entries = [];
        let current = null;
        for (const line of lines) {
          if (line.startsWith(">")) {
            if (current !== null) entries.push(current);
            current = "";
          } else if (current !== null) {
            current += line.trim();
          }
        }
        if (current !== null) entries.push(current);
        seq = entries[0] || "";
      }
      if (!seq) throw new Error("No polymer sequence returned.");
      state.protein = { id: `PDB:${code}`, name: `PDB ${code}`, sequence: seq, true_sse: trueSse, custom: true };
      state.first = true;
      closePdb();
      computeCols();
      render();
      window.scrollTo({ top: 0, behavior: REDUCE ? "auto" : "smooth" });
    } catch (err) {
      warn.className = "seq-warn err";
      warn.textContent = err.message || "Unable to fetch that PDB entry.";
    }
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
    render();
    initTooltip();
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
    });

    el("#reveal").addEventListener("click", e => {
      state.reveal = !state.reveal;
      e.currentTarget.setAttribute("aria-pressed", String(state.reveal));
      el("#tracks").classList.toggle("truthhidden", !state.reveal);
    });

    el("#info").addEventListener("click", () => togglePop(true));
    el("#pop-close").addEventListener("click", () => togglePop(false));
    el("#pop-scrim").addEventListener("click", () => togglePop(false));

    el("#import-btn").addEventListener("click", openPdb);
    el("#pdb-close").addEventListener("click", closePdb);
    el("#pdb-scrim").addEventListener("click", closePdb);
    el("#pdb-go").addEventListener("click", importPdb);
    el("#pdb-input").addEventListener("input", () => {
      const w = el("#pdb-warn"); w.className = "seq-warn"; w.textContent = "";
    });

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

    addEventListener("keydown", e => { if (e.key === "Escape") { togglePop(false); closeSeq(); } });

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

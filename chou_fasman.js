/*
 * chou_fasman.js — browser port of Chou-Fasman/chou_fasman.py (H/E/C).
 *
 * nucleate -> extend -> ACCEPT -> resolve. The web app recomputes the prediction
 * live from each sequence, so predict() MUST match the Python reference residue-
 * for-residue (cross-checked on all chains). predictDetailed() additionally
 * exposes the per-residue algorithm state (Pa/Pb, nucleation, acceptance) that the
 * inspector (Feature B) and confidence shading (Feature D) render — it shares the
 * exact same core() so it can never diverge from predict().
 *
 * Works in the browser (window.ChouFasman) and in Node (module.exports).
 */
(function (root) {
  "use strict";

  // aa -> [Pa(alpha), Pb(beta), Pturn, f0, f1, f2, f3]; only Pa, Pb are used.
  const PROP = {
    A: [1.42, 0.83, 0.66, 0.06, 0.076, 0.035, 0.058],
    C: [0.70, 1.19, 1.19, 0.149, 0.05, 0.117, 0.128],
    D: [1.01, 0.54, 1.46, 0.147, 0.11, 0.179, 0.081],
    E: [1.51, 0.37, 0.74, 0.056, 0.06, 0.077, 0.064],
    F: [1.13, 1.38, 0.60, 0.059, 0.041, 0.065, 0.065],
    G: [0.57, 0.75, 1.56, 0.102, 0.085, 0.19, 0.152],
    H: [1.00, 0.87, 0.95, 0.14, 0.047, 0.093, 0.054],
    I: [1.08, 1.60, 0.47, 0.043, 0.034, 0.013, 0.056],
    K: [1.14, 0.74, 1.01, 0.055, 0.115, 0.072, 0.095],
    L: [1.21, 1.30, 0.59, 0.061, 0.025, 0.036, 0.07],
    M: [1.45, 1.05, 0.60, 0.068, 0.082, 0.014, 0.055],
    N: [0.67, 0.89, 1.56, 0.161, 0.083, 0.191, 0.091],
    P: [0.57, 0.55, 1.52, 0.102, 0.301, 0.034, 0.068],
    Q: [1.11, 1.10, 0.98, 0.074, 0.098, 0.037, 0.098],
    R: [0.98, 0.93, 0.95, 0.07, 0.106, 0.099, 0.085],
    S: [0.77, 0.75, 1.43, 0.12, 0.139, 0.125, 0.106],
    T: [0.83, 1.19, 0.96, 0.086, 0.108, 0.065, 0.079],
    V: [1.06, 1.70, 0.50, 0.062, 0.048, 0.028, 0.053],
    W: [1.08, 1.37, 0.96, 0.077, 0.013, 0.064, 0.167],
    Y: [0.69, 1.47, 1.14, 0.082, 0.065, 0.114, 0.125],
  };

  const PARAMS = {
    H_NUC_WINDOW: 6, H_NUC_MIN_FORMERS: 4,
    E_NUC_WINDOW: 5, E_NUC_MIN_FORMERS: 3,
    FORMER_THRESH: 1.00, EXTEND_THRESH: 1.00,
    H_ACCEPT: 1.03, E_ACCEPT: 1.10,
  };

  function col(seq, c, dflt) {
    const out = new Array(seq.length);
    for (let i = 0; i < seq.length; i++) { const r = PROP[seq[i]]; out[i] = r ? r[c] : dflt; }
    return out;
  }

  // Neumaier compensated sum — mirrors chou_fasman.py's _nsum for bit-for-bit parity.
  function nsum(v, lo, hi) {
    let total = 0.0, c = 0.0;
    for (let i = lo; i < hi; i++) {
      const x = v[i], t = total + x;
      if (Math.abs(total) >= Math.abs(x)) c += (total - t) + x;
      else c += (x - t) + total;
      total = t;
    }
    return total + c;
  }
  function mean(v, lo, hi) { return hi > lo ? nsum(v, lo, hi) / (hi - lo) : 0.0; }

  // nucleate + extend; returns candidate mask AND the nucleation-seed mask.
  function analyzeStruct(P, window, minFormers, former, extend) {
    const n = P.length, cand = new Array(n).fill(false), nuc = new Array(n).fill(false);
    let i = 0;
    while (i <= n - window) {
      let f = 0;
      for (let j = i; j < i + window; j++) if (P[j] >= former) f++;
      if (f >= minFormers) {
        for (let j = i; j < i + window; j++) nuc[j] = true;
        let s = i, e = i + window;
        while (s - 1 >= 0 && s + 3 <= n && mean(P, s - 1, s + 3) >= extend) s--;
        while (e < n && e - 3 >= 0 && mean(P, e - 3, e + 1) >= extend) e++;
        for (let j = s; j < e; j++) cand[j] = true;
        i = e;
      } else i++;
    }
    return { cand, nuc };
  }

  function runs(mark) {
    const out = [];
    let j = 0;
    const n = mark.length;
    while (j < n) {
      if (mark[j]) { const s = j; while (j < n && mark[j]) j++; out.push([s, j]); }
      else j++;
    }
    return out;
  }

  // the single source of truth: returns the prediction AND all per-residue state.
  function core(seq, params) {
    const p = Object.assign({}, PARAMS, params || {});
    const n = seq.length;
    const Pa = col(seq, 0, 1.0), Pb = col(seq, 1, 1.0);
    const ha = analyzeStruct(Pa, p.H_NUC_WINDOW, p.H_NUC_MIN_FORMERS, p.FORMER_THRESH, p.EXTEND_THRESH);
    const ea = analyzeStruct(Pb, p.E_NUC_WINDOW, p.E_NUC_MIN_FORMERS, p.FORMER_THRESH, p.EXTEND_THRESH);
    const H = ha.cand.slice(), E = ea.cand.slice();
    for (const [s, e] of runs(H))
      if (!(mean(Pa, s, e) > p.H_ACCEPT && mean(Pa, s, e) > mean(Pb, s, e)))
        for (let j = s; j < e; j++) H[j] = false;
    for (const [s, e] of runs(E))
      if (!(mean(Pb, s, e) > p.E_ACCEPT && mean(Pb, s, e) > mean(Pa, s, e)))
        for (let j = s; j < e; j++) E[j] = false;
    let pred = "";
    for (let j = 0; j < n; j++) {
      if (H[j] && E[j]) pred += Pa[j] >= Pb[j] ? "H" : "E";
      else if (H[j]) pred += "H";
      else if (E[j]) pred += "E";
      else pred += "C";
    }
    return { pred, Pa, Pb, candH: ha.cand, candE: ea.cand, nucH: ha.nuc, nucE: ea.nuc, accH: H, accE: E };
  }

  function predict(seq, params) { return seq.length ? core(seq, params).pred : ""; }
  function predict3(seq, params) { return predict(seq, params); } // already 3-state

  // per-residue algorithm state for the inspector + confidence shading (Features B/D)
  function predictDetailed(seq, params) {
    if (!seq.length) return { pred: "", meta: [] };
    const r = core(seq, params), n = seq.length, meta = new Array(n);
    for (let i = 0; i < n; i++) {
      const candidate = r.candH[i] || r.candE[i];
      const accepted = r.accH[i] || r.accE[i];
      meta[i] = {
        pa: r.Pa[i], pb: r.Pb[i],
        inNucleus: r.nucH[i] || r.nucE[i],
        regionAccepted: accepted,
        reverted: candidate && !accepted,
      };
    }
    return { pred: r.pred, meta };
  }

  const api = { predict, predict3, predictDetailed, PARAMS, PROP,
    _internals: { col, mean, nsum, analyzeStruct, runs, core } };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.ChouFasman = api;
})(typeof window !== "undefined" ? window : globalThis);

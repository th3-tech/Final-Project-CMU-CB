/*
 * chou_fasman.js — faithful browser port of Chou-Fasman/chou_fasman.py.
 *
 * The web app recomputes the prediction live from each sequence rather than
 * storing it, so this MUST match the Python reference residue-for-residue.
 * Cross-checked by cross_check.js. Keep PARAMS and the propensity table in sync
 * with the Python module and data/propensities.txt.
 *
 * Works both in the browser (attaches window.ChouFasman) and in Node
 * (module.exports) so the cross-check can import it.
 */
(function (root) {
  "use strict";

  // aa -> [alpha, beta, turn, f0, f1, f2, f3]
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
    H_FORMER_THRESH: 1.00, E_FORMER_THRESH: 1.00,
    H_EXTEND_THRESH: 1.00, E_EXTEND_THRESH: 1.00,
    H_ACCEPT: 1.05, E_ACCEPT: 1.08,  // only calibrated deviation from canonical (1.03/1.05)
    TURN_PROB_THRESH: 7.5e-5,
  };

  function col(seq, c, dflt) {
    const out = new Array(seq.length);
    for (let i = 0; i < seq.length; i++) {
      const row = PROP[seq[i]];
      out[i] = row ? row[c] : dflt;
    }
    return out;
  }

  // Neumaier compensated sum — mirrors chou_fasman.py's _nsum so the two ports are
  // bit-for-bit identical. Needed because region averages can land exactly on an
  // acceptance threshold, where a 1-ULP summation difference flips the result.
  function nsum(v, lo, hi) {
    let total = 0.0, c = 0.0;
    for (let i = lo; i < hi; i++) {
      const x = v[i];
      const t = total + x;
      if (Math.abs(total) >= Math.abs(x)) c += (total - t) + x;
      else c += (x - t) + total;
      total = t;
    }
    return total + c;
  }

  function mean(v, lo, hi) {
    return nsum(v, lo, hi) / (hi - lo);
  }

  function contiguousRegions(mask) {
    const regions = [];
    let start = null;
    for (let i = 0; i < mask.length; i++) {
      if (mask[i] && start === null) start = i;
      else if (!mask[i] && start !== null) { regions.push([start, i]); start = null; }
    }
    if (start !== null) regions.push([start, mask.length]);
    return regions;
  }

  function nucleateAndExtend(P, window, minFormers, formerThresh, extendThresh) {
    const n = P.length;
    const seed = new Array(n).fill(false);
    for (let i = 0; i <= n - window; i++) {
      let formers = 0;
      for (let j = i; j < i + window; j++) if (P[j] >= formerThresh) formers++;
      if (formers >= minFormers) for (let j = i; j < i + window; j++) seed[j] = true;
    }
    const regions = [];
    for (let [start, end] of contiguousRegions(seed)) {
      while (start - 1 >= 0 && start + 3 <= n && mean(P, start - 1, start + 3) >= extendThresh) start--;
      while (end + 1 <= n && end - 3 >= 0 && mean(P, end - 3, end + 1) >= extendThresh) end++;
      regions.push([start, end]);
    }
    return regions;
  }

  function accept(regions, Pself, Pother, acceptThresh) {
    const kept = [];
    for (const [s, e] of regions) {
      const a = mean(Pself, s, e), b = mean(Pother, s, e);
      if (a > acceptThresh && a > b) kept.push([s, e]);
    }
    return kept;
  }

  function predictTurns(seq, Pa, Pb, Pt, thresh) {
    const n = seq.length;
    const f0 = col(seq, 3, 0), f1 = col(seq, 4, 0), f2 = col(seq, 5, 0), f3 = col(seq, 6, 0);
    const turn = new Array(n).fill(false);
    for (let i = 0; i < n - 3; i++) {
      const prob = f0[i] * f1[i + 1] * f2[i + 2] * f3[i + 3];
      const at = mean(Pt, i, i + 4), aa = mean(Pa, i, i + 4), ab = mean(Pb, i, i + 4);
      if (prob > thresh && at > 1.0 && at > aa && at > ab) {
        for (let j = i; j < i + 4; j++) turn[j] = true;
      }
    }
    return turn;
  }

  function predict(seq, params) {
    const p = Object.assign({}, PARAMS, params || {});
    const n = seq.length;
    if (n === 0) return "";
    const Pa = col(seq, 0, 1.0), Pb = col(seq, 1, 1.0), Pt = col(seq, 2, 1.0);

    const helixRegions = accept(
      nucleateAndExtend(Pa, p.H_NUC_WINDOW, p.H_NUC_MIN_FORMERS, p.H_FORMER_THRESH, p.H_EXTEND_THRESH),
      Pa, Pb, p.H_ACCEPT);
    const sheetRegions = accept(
      nucleateAndExtend(Pb, p.E_NUC_WINDOW, p.E_NUC_MIN_FORMERS, p.E_FORMER_THRESH, p.E_EXTEND_THRESH),
      Pb, Pa, p.E_ACCEPT);

    const helix = new Array(n).fill(false);
    const sheet = new Array(n).fill(false);
    for (const [s, e] of helixRegions) for (let i = s; i < e; i++) helix[i] = true;
    for (const [s, e] of sheetRegions) for (let i = s; i < e; i++) sheet[i] = true;

    const overlap = new Array(n);
    for (let i = 0; i < n; i++) overlap[i] = helix[i] && sheet[i];
    for (const [s, e] of contiguousRegions(overlap)) {
      if (nsum(Pa, s, e) >= nsum(Pb, s, e)) for (let i = s; i < e; i++) sheet[i] = false;
      else for (let i = s; i < e; i++) helix[i] = false;
    }

    const turn = predictTurns(seq, Pa, Pb, Pt, p.TURN_PROB_THRESH);

    let out = "";
    for (let i = 0; i < n; i++) {
      if (helix[i]) out += "H";
      else if (sheet[i]) out += "E";
      else if (turn[i]) out += "T";
      else out += "C";
    }
    return out;
  }

  function predict3(seq, params) {
    return predict(seq, params).replace(/T/g, "C");
  }

  const api = { predict, predict3, PARAMS, PROP,
    _internals: { col, mean, nucleateAndExtend, accept } };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.ChouFasman = api;
})(typeof window !== "undefined" ? window : globalThis);

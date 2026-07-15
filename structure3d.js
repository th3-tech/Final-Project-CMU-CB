/*
 * structure3d.js — Feature E: render a protein's backbone cartoon coloured by
 * Chou-Fasman correctness, with chameleon residues highlighted.
 *
 * Structures are pre-stripped backbone PDBs in structures/<ID>.pdb (residues
 * renumbered to sequence position, so colouring is index-by-index). Secondary
 * structure isn't in the file — we set each atom's .ss from the shipped truth,
 * so 3Dmol still draws helix ribbons and sheet arrows from a light bundle.
 *
 * Structure3D.render(container, opts) -> Promise<viewer|null>
 *   opts: { id, sequence, trueSse, pred, isCham(i)->bool, spin?, bg? }
 */
(function (root) {
  "use strict";
  const OK = "#3FB37F", BAD = "#C0392B", CHAM = "#E84AC4", GREY = "#3A4354";
  const cache = {};
  const pending = {};

  // structures/<ID>.js is JSONP: __STRUCT_CB("<ID>", "<pdb text>"). We load it with an
  // injected <script> (not fetch) so it works from file:// — Chrome blocks fetch of local
  // files, which is why double-clicking index.html broke the 3D but the preview server didn't.
  window.__STRUCT_CB = function (id, pdb) {
    cache[id] = pdb;
    const r = pending[id];
    if (r) { delete pending[id]; r(pdb); }
  };

  function fetchPdb(id) {
    if (cache[id]) return Promise.resolve(cache[id]);
    return new Promise((resolve, reject) => {
      const prev = pending[id];
      pending[id] = prev ? p => { prev(p); resolve(p); } : resolve;
      if (prev) return;                       // script already injected for this id
      const s = document.createElement("script");
      s.src = "structures/" + id + ".js";
      s.onerror = () => { delete pending[id]; reject(new Error("no-structure")); };
      document.head.appendChild(s);
    });
  }

  function render(container, opts) {
    container.innerHTML = '<div class="s3d-msg">loading structure…</div>';
    return fetchPdb(opts.id).then(data => {
      container.innerHTML = "";
      const viewer = $3Dmol.createViewer(container, { backgroundColor: opts.bg || "#080b11" });
      const model = viewer.addModel(data, "pdb");
      const truth = opts.trueSse;
      model.selectedAtoms({}).forEach(a => {
        const s = truth[a.resi - 1];
        a.ss = s === "H" ? "h" : s === "E" ? "s" : "c";
      });
      const ok = [], bad = [], cham = [];
      for (let i = 0; i < opts.sequence.length; i++) {
        const resi = i + 1;
        if (opts.isCham(i)) cham.push(resi);
        else if (opts.pred[i] === truth[i]) ok.push(resi);
        else bad.push(resi);
      }
      viewer.setStyle({}, { cartoon: { color: GREY } });
      viewer.setStyle({ resi: ok }, { cartoon: { color: OK } });
      viewer.setStyle({ resi: bad }, { cartoon: { color: BAD } });
      viewer.setStyle({ resi: cham }, { cartoon: { color: CHAM, thickness: 1.0 } });
      viewer.zoomTo();
      viewer.render();
      if (opts.spin !== false) viewer.spin("y", 0.4);
      return viewer;
    }).catch(() => {
      container.innerHTML = '<div class="s3d-msg">3D structure unavailable for this chain</div>';
      return null;
    });
  }

  // pulse a single residue (2D -> 3D link): a temporary sphere shape at its CA
  function pulse(viewer, resi) {
    if (!viewer) return;
    const ca = viewer.getModel().selectedAtoms({ resi: resi, atom: "CA" })[0];
    if (!ca) return;
    const sph = viewer.addSphere({ center: { x: ca.x, y: ca.y, z: ca.z }, radius: 2.4,
      color: "#7C74FF", opacity: 0.9 });
    viewer.render();
    setTimeout(() => { viewer.removeShape(sph); viewer.render(); }, 1100);
  }

  root.Structure3D = { render, pulse, COLORS: { OK, BAD, CHAM, GREY } };
})(window);

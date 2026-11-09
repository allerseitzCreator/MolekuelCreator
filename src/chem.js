// Chemie-Kern: SMILES lesen, kekulisieren, 2D-Layout, IUPAC-Namen übersetzen
const Chem = (() => {
'use strict';

/* ======================= Graph ======================= */
class G {
  constructor() { this.atoms = []; this.bonds = []; }
  add(el, o = {}) { this.atoms.push({el, charge: o.charge || 0, arom: !!o.arom, hx: o.hx}); return this.atoms.length - 1; }
  bond(a, b, order = 1) {
    if (a === b) throw new Error('Bindung an sich selbst');
    const ex = this.bonds.find(x => (x.a === a && x.b === b) || (x.a === b && x.b === a));
    if (ex) { if (order !== 'ar' && ex.order !== 'ar') ex.order = Math.max(ex.order, order); return ex; }
    const bd = {a, b, order}; this.bonds.push(bd); return bd;
  }
  nb(i) { return this.bonds.filter(b => b.a === i || b.b === i).map(b => b.a === i ? b.b : b.a); }
  bondsOf(i) { return this.bonds.filter(b => b.a === i || b.b === i); }
  // fremden Graphen einhängen, gibt Index-Verschiebung zurück
  merge(o) { const off = this.atoms.length; o.atoms.forEach(a => { const c = {...a}; if (c.nbOrder) c.nbOrder = c.nbOrder.map(n => n >= 0 ? n + off : n); this.atoms.push(c); }); o.bonds.forEach(b => this.bonds.push({a: b.a + off, b: b.b + off, order: b.order}));
    if (o.ezSpec) this.ezSpec = [...(this.ezSpec || []), ...o.ezSpec.map(e => ({a: e.a + off, b: e.b + off, x: e.x + off, y: e.y + off, cis: e.cis}))];
    return off; }
}

/* ======================= SMILES ======================= */
function parseSmiles(s) {
  s = s.trim();
  if (!s) throw new Error('leer');
  const g = new G();
  let prev = null, pend = null, pendDir = null, i = 0;
  const stack = [], ring = {}, dirBonds = [];   // / und \ (E/Z an Doppelbindungen)
  const ORG = ['Cl', 'Br', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I', 'b', 'c', 'n', 'o', 'p', 's'];
  const nbOrder = [];          // Nachbarreihenfolge je Atom (für @/@@)
  const bondTo = (a, b, order) => {
    if (order == null) order = (g.atoms[a].arom && g.atoms[b].arom) ? 'ar' : 1;
    g.bond(a, b, order);
  };
  while (i < s.length) {
    const ch = s[i];
    if (ch === '(') { if (prev === null) throw new Error('Klammer am Anfang'); stack.push(prev); i++; continue; }
    if (ch === ')') { if (!stack.length) throw new Error('Klammer'); prev = stack.pop(); i++; continue; }
    if (ch === '.') { prev = null; i++; continue; }
    if ('-=#$:/\\'.includes(ch)) { pend = ch === '=' ? 2 : ch === '#' ? 3 : ch === ':' ? 'ar' : 1; if (ch === '/' || ch === '\\') pendDir = ch; i++; continue; }
    if (/[0-9%]/.test(ch)) {
      let n;
      if (ch === '%') { n = s.substr(i + 1, 2); i += 3; } else { n = ch; i++; }
      if (prev === null) throw new Error('Ringziffer ohne Atom');
      if (ring[n]) {
        const r = ring[n]; bondTo(r.atom, prev, pend ?? r.order);
        nbOrder[r.atom][r.slot] = prev; nbOrder[prev].push(r.atom); delete ring[n];
      } else { ring[n] = {atom: prev, order: pend, slot: nbOrder[prev].length}; nbOrder[prev].push(null); }
      pend = null; pendDir = null; continue;
    }
    let idx;
    if (ch === '[') {
      const j = s.indexOf(']', i); if (j < 0) throw new Error(']');
      const t = s.slice(i + 1, j); i = j + 1;
      const m = t.match(/^(\d*)([A-Z][a-z]?|[a-z][a-z]?)(@{0,2})(H\d*)?(\+\d*|-\d*|\++|-+)?(:\d+)?$/);
      if (!m) throw new Error('Atom [' + t + ']');
      let el = m[2]; const arom = el[0] === el[0].toLowerCase();
      el = el[0].toUpperCase() + el.slice(1);
      let charge = 0;
      if (m[5]) { const sg = m[5][0] === '+' ? 1 : -1; const rest = m[5].slice(1); charge = sg * (rest === '' ? 1 : /^\d+$/.test(rest) ? +rest : m[5].length); }
      const hx = m[4] ? (m[4].length > 1 ? +m[4].slice(1) : 1) : 0;
      idx = g.add(el, {charge, arom, hx});
      if (m[3]) g.atoms[idx].chiralTag = m[3];
      nbOrder[idx] = prev !== null ? [prev] : [];
      if (hx) for (let k = 0; k < hx; k++) nbOrder[idx].push(-1);
    } else {
      let el = s.substr(i, 2);
      if (el === 'Cl' || el === 'Br') i += 2; else { el = ch; i++; }
      if (!ORG.includes(el)) throw new Error('Unbekanntes Zeichen ' + el);
      const arom = el === el.toLowerCase();
      idx = g.add(el[0].toUpperCase() + el.slice(1), {arom});
      nbOrder[idx] = prev !== null ? [prev] : [];
    }
    if (prev !== null) { bondTo(prev, idx, pend); nbOrder[prev].push(idx); if (pendDir) dirBonds.push({u: prev, v: idx, s: pendDir}); }
    pend = null; pendDir = null; prev = idx;
  }
  if (Object.keys(ring).length) throw new Error('offener Ring');
  if (stack.length) throw new Error('offene Klammer');
  g.atoms.forEach((a, i) => { if (a.chiralTag) a.nbOrder = nbOrder[i]; });
  // cis/trans aus / und \ : „/“ = das später geschriebene Atom liegt höher
  if (dirBonds.length) {
    const sideOf = (A, d) => { const X = d.u === A ? d.v : d.u, up = d.s === '/' ? 1 : -1; return {X, side: d.u === A ? up : -up}; };
    g.ezSpec = [];
    for (const b of g.bonds) {
      if (b.order !== 2) continue;
      const da = dirBonds.find(d => (d.u === b.a || d.v === b.a) && !(d.u === b.b || d.v === b.b));
      const db = dirBonds.find(d => (d.u === b.b || d.v === b.b) && !(d.u === b.a || d.v === b.a));
      if (!da || !db) continue;
      const p = sideOf(b.a, da), q = sideOf(b.b, db);
      g.ezSpec.push({a: b.a, b: b.b, x: p.X, y: q.X, cis: p.side === q.side});
    }
  }
  return g;
}

const BASEVAL = {B: 3, C: 4, Si: 4, N: 3, P: 3, O: 2, S: 2, F: 1, Cl: 1, Br: 1, I: 1, H: 1};
// aromatische Bindungen in Einfach/Doppel umwandeln
function kekulize(g) {
  const arB = g.bonds.filter(b => b.order === 'ar');
  if (!arB.length) return g;
  const need = new Set();
  g.atoms.forEach((a, i) => {
    if (!g.bondsOf(i).some(b => b.order === 'ar')) return;
    let v = BASEVAL[a.el] ?? 4;
    if (a.el === 'C') v -= Math.abs(a.charge); else if (a.el === 'B') v -= a.charge; else v += a.charge;
    let used = (a.hx || 0);
    for (const b of g.bondsOf(i)) used += b.order === 'ar' ? 1 : b.order;
    if (v - used >= 1) need.add(i);
  });
  // perfekte Paarung der „bedürftigen“ Atome über aromatische Bindungen (Backtracking)
  const matched = new Map();
  const cands = i => arB.filter(b => (b.a === i || b.b === i)).map(b => [b, b.a === i ? b.b : b.a]).filter(([, j]) => need.has(j));
  const list = [...need].sort((x, y) => cands(x).length - cands(y).length);
  let steps = 0;
  function solve(k) {
    if (++steps > 200000) return false;
    while (k < list.length && matched.has(list[k])) k++;
    if (k >= list.length) return true;
    const i = list[k];
    for (const [b, j] of cands(i)) {
      if (matched.has(j)) continue;
      matched.set(i, b); matched.set(j, b);
      if (solve(k + 1)) return true;
      matched.delete(i); matched.delete(j);
    }
    return false;
  }
  if (!solve(0)) throw new Error('Aromatisches System lässt sich nicht kekulisieren');
  const dbl = new Set(matched.values());
  for (const b of arB) b.order = dbl.has(b) ? 2 : 1;
  return g;
}

/* ======================= Ringe ======================= */
function findRings(g) {
  const n = g.atoms.length, E = g.bonds;
  const adj = Array.from({length: n}, () => []);
  E.forEach((b, k) => { adj[b.a].push([b.b, k]); adj[b.b].push([b.a, k]); });
  // Anzahl Komponenten
  const comp = new Array(n).fill(-1); let nc = 0;
  for (let i = 0; i < n; i++) if (comp[i] < 0) { const st = [i]; comp[i] = nc; while (st.length) { const u = st.pop(); for (const [v] of adj[u]) if (comp[v] < 0) { comp[v] = nc; st.push(v); } } nc++; }
  const target = E.length - n + nc;
  if (target <= 0) return [];
  const cycles = [], seen = new Set();
  E.forEach((bd, k) => {
    // kürzester Weg a→b ohne diese Bindung
    const prev = new Map([[bd.a, null]]), q = [bd.a];
    while (q.length) {
      const u = q.shift(); if (u === bd.b) break;
      for (const [v, kk] of adj[u]) if (kk !== k && !prev.has(v)) { prev.set(v, [u, kk]); q.push(v); }
    }
    if (!prev.has(bd.b)) return;
    const path = [], edges = [k]; let c = bd.b;
    while (c !== null) { path.push(c); const p = prev.get(c); if (!p) break; edges.push(p[1]); c = p[0]; }
    const key = [...path].sort((x, y) => x - y).join(',');
    if (seen.has(key)) return; seen.add(key);
    let mask = 0n; for (const e of edges) mask |= 1n << BigInt(e);
    cycles.push({atoms: path, mask});
  });
  cycles.sort((x, y) => x.atoms.length - y.atoms.length);
  const basis = [], out = [];
  for (const c of cycles) {
    let m = c.mask;
    for (const b of basis) if ((m ^ b) < m) m ^= b;
    if (m !== 0n) { basis.push(m); basis.sort((x, y) => (y > x ? 1 : -1)); out.push(c.atoms); if (out.length >= target) break; }
  }
  return out; // Atome jeweils in Ringreihenfolge
}

/* ======================= Layout ======================= */
const TAU = 2 * Math.PI;
const nrm = x => ((x % TAU) + TAU) % TAU;
function spread(occ, k) {
  if (k <= 0) return [];
  if (!occ.length) return Array.from({length: k}, (_, i) => -Math.PI / 6 + i * TAU / k);
  const s = occ.map(nrm).sort((p, q) => p - q);
  const gaps = s.map((v, i) => ({start: v, size: (i === s.length - 1 ? s[0] + TAU : s[i + 1]) - v, m: 0}));
  for (let i = 0; i < k; i++) {
    let best = gaps[0];
    for (const g of gaps) if (g.size / (g.m + 2) > best.size / (best.m + 2) + 1e-9) best = g;
    best.m++;
  }
  const out = [];
  for (const g of gaps) for (let j = 1; j <= g.m; j++) out.push(g.start + g.size * j / (g.m + 1));
  return out;
}

function layout(g) {
  const n = g.atoms.length;
  const nb = Array.from({length: n}, (_, i) => g.nb(i));
  const bondOrder = (a, b) => { const x = g.bonds.find(q => (q.a === a && q.b === b) || (q.a === b && q.b === a)); return x ? x.order : 0; };
  const rings = findRings(g);
  // Ringsysteme: Ringe, die Atome teilen
  const sysOf = new Array(rings.length).fill(-1); const systems = [];
  rings.forEach((r, i) => {
    if (sysOf[i] >= 0) return;
    const sid = systems.length, st = [i]; sysOf[i] = sid; const members = [];
    while (st.length) { const k = st.pop(); members.push(k); rings.forEach((r2, j) => { if (sysOf[j] < 0 && r2.some(a => rings[k].includes(a))) { sysOf[j] = sid; st.push(j); } }); }
    systems.push(members.map(k => rings[k]));
  });
  const atomSys = new Array(n).fill(-1);
  systems.forEach((sys, s) => sys.forEach(r => r.forEach(a => atomSys[a] = s)));

  // Ringsystem in lokalen Koordinaten (Bindungslänge 1)
  function placeSystem(sys) {
    const P = new Map();
    const ringsLeft = [...sys].sort((x, y) => sharedCount(y) - sharedCount(x) || y.length - x.length);
    function sharedCount(r) { return sys.filter(o => o !== r && o.some(a => r.includes(a))).length; }
    const first = ringsLeft.shift();
    const R0 = 1 / (2 * Math.sin(Math.PI / first.length));
    first.forEach((a, k) => { const ang = -Math.PI / 2 + Math.PI / first.length + k * TAU / first.length; P.set(a, {x: R0 * Math.cos(ang), y: R0 * Math.sin(ang)}); });
    let guard = 0;
    while (ringsLeft.length && guard++ < 200) {
      ringsLeft.sort((x, y) => y.filter(a => P.has(a)).length - x.filter(a => P.has(a)).length);
      const r = ringsLeft.shift(), m = r.length;
      const placed = r.filter(a => P.has(a));
      const cen = () => { let x = 0, y = 0; for (const p of P.values()) { x += p.x; y += p.y; } return {x: x / P.size, y: y / P.size}; };
      const R = 1 / (2 * Math.sin(Math.PI / m)), ap = 1 / (2 * Math.tan(Math.PI / m));
      if (placed.length >= 2) {
        // zwei benachbarte, bereits platzierte Ringatome suchen
        let i0 = -1;
        for (let i = 0; i < m; i++) if (P.has(r[i]) && P.has(r[(i + 1) % m]) && !P.has(r[(i + 2) % m])) { i0 = i; break; }
        if (i0 < 0) for (let i = 0; i < m; i++) if (P.has(r[i]) && P.has(r[(i + 1) % m])) { i0 = i; break; }
        if (i0 < 0) { placeLoose(r, P); continue; }
        // Kette platzierter Atome bis zum ersten unplatzierten finden
        const A = P.get(r[i0]), B = P.get(r[(i0 + 1) % m]);
        const mx = (A.x + B.x) / 2, my = (A.y + B.y) / 2, L = Math.hypot(B.x - A.x, B.y - A.y) || 1;
        let nx = -(B.y - A.y) / L, ny = (B.x - A.x) / L;
        const c = cen();
        if ((c.x - mx) * nx + (c.y - my) * ny > 0) { nx = -nx; ny = -ny; }
        const cx = mx + nx * ap, cy = my + ny * ap;
        const a0 = Math.atan2(A.y - cy, A.x - cx), a1 = Math.atan2(B.y - cy, B.x - cx);
        let d = a1 - a0; while (d > Math.PI) d -= TAU; while (d <= -Math.PI) d += TAU;
        const st = (d < 0 ? -1 : 1) * TAU / m;
        for (let k = 0; k < m; k++) { const a = r[(i0 + k) % m]; if (!P.has(a)) P.set(a, {x: cx + R * Math.cos(a0 + k * st), y: cy + R * Math.sin(a0 + k * st)}); }
      } else if (placed.length === 1) {
        const a = placed[0], A = P.get(a), c = cen();
        let dx = A.x - c.x, dy = A.y - c.y; const l = Math.hypot(dx, dy) || 1; dx /= l; dy /= l;
        const cx = A.x + dx * R, cy = A.y + dy * R, i0 = r.indexOf(a), a0 = Math.atan2(A.y - cy, A.x - cx);
        for (let k = 1; k < m; k++) { const b = r[(i0 + k) % m]; if (!P.has(b)) P.set(b, {x: cx + R * Math.cos(a0 + k * TAU / m), y: cy + R * Math.sin(a0 + k * TAU / m)}); }
      } else ringsLeft.push(r);
    }
    return P;
  }
  function placeLoose(r, P) {
    // Notlösung für überbrückte Ringe: fehlende Atome zwischen platzierten Nachbarn verteilen
    const m = r.length;
    for (let k = 0; k < m; k++) if (!P.has(r[k])) {
      let i = k; while (!P.has(r[(i - 1 + m) % m]) && i - k > -m) i--;
      const prevA = P.get(r[(k - 1 + m) % m]) || [...P.values()][0];
      P.set(r[k], {x: prevA.x + 0.9, y: prevA.y + 0.3});
    }
  }

  const pos = new Array(n).fill(null);
  const parent = new Array(n).fill(-1);
  const done = new Array(n).fill(false);

  function placeSystemAt(s, attachAtom, ax, ay, dirAngle) {
    const P = placeSystem(systems[s]);
    const A = P.get(attachAtom);
    let cx = 0, cy = 0; for (const p of P.values()) { cx += p.x; cy += p.y; } cx /= P.size; cy /= P.size;
    const cur = Math.atan2(cy - A.y, cx - A.x);
    const rot = dirAngle === null ? 0 : dirAngle - cur;
    const cs = Math.cos(rot), sn = Math.sin(rot);
    for (const [a, p] of P) {
      const dx = p.x - A.x, dy = p.y - A.y;
      pos[a] = {x: ax + dx * cs - dy * sn, y: ay + dx * sn + dy * cs};
    }
    return [...P.keys()];
  }

  const collides = (p, self) => pos.some((q, i) => q && i !== self && Math.hypot(q.x - p.x, q.y - p.y) < 0.8);

  function expand(start) {
    const stack = [start];
    while (stack.length) {
      const u = stack.pop();
      if (done[u]) continue; done[u] = true;
      const placedN = nb[u].filter(v => pos[v]);
      const un = nb[u].filter(v => !pos[v]);
      if (!un.length) continue;
      const occ = placedN.map(v => Math.atan2(pos[v].y - pos[u].y, pos[v].x - pos[u].x));
      let angs;
      const linear = v => bondOrder(u, v) === 3 || (placedN.length === 1 && bondOrder(u, placedN[0]) === 3) ||
        (placedN.length === 1 && bondOrder(u, placedN[0]) === 2 && bondOrder(u, v) === 2);
      if (placedN.length === 1 && un.length === 1 && atomSys[u] < 0) {
        const th = occ[0];
        if (linear(un[0])) angs = [th + Math.PI];
        else {
          const c1 = th + TAU / 3, c2 = th - TAU / 3;
          const p = placedN[0], gp = parent[p];
          const pt = a => ({x: pos[u].x + Math.cos(a), y: pos[u].y + Math.sin(a)});
          if (gp >= 0 && pos[gp]) {
            const d1 = Math.hypot(pt(c1).x - pos[gp].x, pt(c1).y - pos[gp].y), d2 = Math.hypot(pt(c2).x - pos[gp].x, pt(c2).y - pos[gp].y);
            angs = [d1 >= d2 ? c1 : c2];
          } else angs = [Math.cos(c1) >= Math.cos(c2) ? c1 : c2];
        }
      } else if (!placedN.length && un.length === 1) angs = [-Math.PI / 6];
      else {
        angs = spread(occ, un.length);
        if (placedN.length === 1 && atomSys[u] < 0 && un.length > 1) {
          // Hauptkette (größter Ast) zickzackförmig fortsetzen, Substituenten auf die übrigen Plätze
          const bsz = v => { const seen = new Set([u, v]), st = [v]; let c = 0; while (st.length) { const x = st.pop(); c++; for (const y of nb[x]) if (!seen.has(y)) { seen.add(y); st.push(y); } } return c; };
          const sizes = new Map(un.map(v => [v, bsz(v)]));
          un.sort((x, y) => sizes.get(y) - sizes.get(x));
          const gp = parent[placedN[0]];
          if (gp >= 0 && pos[gp]) {
            const dist = a => Math.hypot(pos[u].x + Math.cos(a) - pos[gp].x, pos[u].y + Math.sin(a) - pos[gp].y);
            angs.sort((a1, a2) => dist(a2) - dist(a1));
          }
        }
      }
      un.forEach((v, k) => {
        let a = angs[k];
        let p = {x: pos[u].x + Math.cos(a), y: pos[u].y + Math.sin(a)};
        if (collides(p, v)) {
          for (const d of [0.35, -0.35, 0.7, -0.7, 1.05, -1.05]) {
            const q = {x: pos[u].x + Math.cos(a + d), y: pos[u].y + Math.sin(a + d)};
            if (!collides(q, v)) { a += d; p = q; break; }
          }
        }
        parent[v] = u;
        if (atomSys[v] >= 0 && atomSys[v] !== atomSys[u]) {
          const placedNow = placeSystemAt(atomSys[v], v, p.x, p.y, a);
          placedNow.forEach(x => { if (x !== v) parent[x] = v; });
          // Systematome stapeln (Anker zuletzt, damit er zuerst expandiert)
          placedNow.filter(x => x !== v).forEach(x => stack.push(x));
          stack.push(v);
        } else { pos[v] = p; stack.push(v); }
      });
    }
  }

  // Komponenten nacheinander, nebeneinander
  let offsetX = 0;
  const compDone = new Array(n).fill(false);
  for (let s0 = 0; s0 < n; s0++) {
    if (compDone[s0]) continue;
    // Komponente bestimmen
    const comp = []; const st = [s0]; compDone[s0] = true;
    while (st.length) { const u = st.pop(); comp.push(u); for (const v of nb[u]) if (!compDone[v]) { compDone[v] = true; st.push(v); } }
    let start;
    const sysIn = [...new Set(comp.map(a => atomSys[a]).filter(x => x >= 0))];
    if (sysIn.length) {
      const big = sysIn.sort((x, y) => systems[y].flat().length - systems[x].flat().length)[0];
      const anchor = systems[big][0][0];
      placeSystemAt(big, anchor, 0, 0, null);
      const members = comp.filter(a => atomSys[a] === big);
      members.forEach(a => { if (a !== anchor) parent[a] = anchor; });
      // alle Systematome expandieren
      for (const a of members) expand(a);
      start = null;
    } else {
      // Ende der längsten Kette suchen (doppelte BFS)
      const far = src => { const d = new Map([[src, 0]]), q = [src]; let last = src; while (q.length) { const u = q.shift(); last = u; for (const v of nb[u]) if (!d.has(v)) { d.set(v, d.get(u) + 1); q.push(v); } } return last; };
      start = far(far(comp[0]));
      pos[start] = {x: 0, y: 0};
      expand(start);
    }
    // nicht erreichte Atome (sollte nicht vorkommen)
    comp.forEach(a => { if (!pos[a]) pos[a] = {x: Math.random(), y: Math.random()}; });
    // verschieben, dass Komponenten nebeneinander liegen
    let minX = Infinity, maxX = -Infinity, sy = 0;
    comp.forEach(a => { minX = Math.min(minX, pos[a].x); maxX = Math.max(maxX, pos[a].x); sy += pos[a].y; });
    const dx = offsetX - minX, dy = -sy / comp.length;
    comp.forEach(a => { pos[a] = {x: pos[a].x + dx, y: pos[a].y + dy}; });
    offsetX += maxX - minX + 1.8;
  }
  // cis/trans laut SMILES: falsch liegende Seite an der Doppelbindung spiegeln
  if (g.ezSpec) {
    const ringSet = new Set(rings.flat());
    for (const e of g.ezSpec) {
      if (ringSet.has(e.a) && ringSet.has(e.b)) continue;
      const A = pos[e.a], B = pos[e.b], X = pos[e.x], Y = pos[e.y]; if (!A || !B || !X || !Y) continue;
      const dx = B.x - A.x, dy = B.y - A.y, L2 = dx * dx + dy * dy || 1;
      const side = P => Math.sign(dx * (P.y - A.y) - dy * (P.x - A.x));
      if ((side(X) === side(Y)) === e.cis) continue;
      const br = branchAtoms(g, e.b, e.a); if (br.includes(e.a)) continue;
      for (const m of br) { if (m === e.b) continue; const P = pos[m], t = ((P.x - A.x) * dx + (P.y - A.y) * dy) / L2, fx = A.x + dx * t, fy = A.y + dy * t; pos[m] = {x: 2 * fx - P.x, y: 2 * fy - P.y}; }
    }
  }
  unclash(g, pos, rings);
  peptideLayout(g, pos, rings);
  // funktionelle Gruppen (OH, =O, NH2, Cl, NO2 …) mit kurzer Bindung (¾), Gerüst mit mittlerer Bindung
  for (const fg of functionalGroups(g, rings)) {
    const A = pos[fg.attach];
    fg.members.forEach(m => { pos[m] = {x: A.x + (pos[m].x - A.x) * SHORT_BOND, y: A.y + (pos[m].y - A.y) * SHORT_BOND}; });
  }
  return pos;
}
// Peptide/Proteine: Hauptkette vom N-Terminus (links) zum C-Terminus (rechts); Reste nach oben bzw. unten
function peptideBackbone(g, rings) {
  const RS = new Set((rings || findRings(g)).flat()), n = g.atoms.length, el = i => g.atoms[i].el;
  const ord = (a, b) => (g.bonds.find(x => (x.a === a && x.b === b) || (x.a === b && x.b === a)) || {}).order;
  const heavy = i => g.nb(i).filter(j => el(j) !== 'H');
  const carbonyl = c => el(c) === 'C' && g.nb(c).some(o => el(o) === 'O' && ord(c, o) === 2);
  const nT = [], cT = new Set();
  for (let i = 0; i < n; i++) {
    if (RS.has(i)) continue;
    if (el(i) === 'N' && heavy(i).length === 1 && !g.nb(i).some(carbonyl)) nT.push(i);
    if (carbonyl(i) && g.nb(i).some(o => el(o) === 'O' && ord(i, o) === 1)) cT.add(i);
  }
  let best = null;
  for (const s0 of nT) {
    const prev = new Map([[s0, -1]]), q = [s0];
    while (q.length) { const u = q.shift(); for (const v of heavy(u)) if (!prev.has(v) && !RS.has(v)) { prev.set(v, u); q.push(v); } }
    for (const c of cT) {
      if (!prev.has(c)) continue;
      const path = []; for (let x = c; x !== -1; x = prev.get(x)) path.unshift(x);
      let amides = 0; for (let k = 0; k + 1 < path.length; k++) if (carbonyl(path[k]) && el(path[k + 1]) === 'N') amides++;
      if (amides && (!best || path.length > best.length)) best = path;
    }
  }
  if (!best) return null;
  // C-Terminus: einfach gebundenes O (OH, O⁻, Ester-O) und ggf. dessen Alkylrest verlängern die Kette
  const last = best[best.length - 1], o = g.nb(last).find(x => el(x) === 'O' && ord(last, x) === 1);
  if (o !== undefined) { best.push(o); const r = heavy(o).find(x => x !== last && !RS.has(x)); if (r !== undefined) best.push(r); }
  return best;
}
function peptideLayout(g, pos, rings) {
  const bb = peptideBackbone(g, rings); if (!bb || bb.length < 4) return;
  const inBB = new Set(bb), old = pos.map(p => p && {...p});
  // Lage der Komponente beibehalten
  const comp = branchAtoms(g, bb[0], -1);
  let minX = Infinity, sy = 0; comp.forEach(i => { minX = Math.min(minX, old[i].x); sy += old[i].y; }); sy /= comp.length;
  const np = new Map();
  bb.forEach((id, i) => np.set(id, {x: i * Math.cos(Math.PI / 6), y: i % 2 ? -0.5 : 0.5}));
  for (let i = 0; i < bb.length; i++) {
    const a = bb[i], roots = g.nb(a).filter(r => !inBB.has(r));
    const up = i % 2 === 1, base = up ? -Math.PI / 2 : Math.PI / 2;
    roots.forEach((r, k) => {
      const br = branchAtoms(g, r, a);
      if (br.some(x => inBB.has(x))) return;   // Ring über die Hauptkette (z. B. Prolin): Lage beibehalten
      const dirA = roots.length === 1 ? base : base + (k - (roots.length - 1) / 2) * Math.PI / 3;
      const A = old[a], vr = {x: old[r].x - A.x, y: old[r].y - A.y}, P = np.get(a);
      let bestP = null, bs = Infinity;
      for (const mir of [false, true]) {
        const ang0 = Math.atan2(mir ? -vr.y : vr.y, vr.x), rot = dirA - ang0, c = Math.cos(rot), sn = Math.sin(rot);
        const Q = br.map(m => { let x = old[m].x - A.x, y = old[m].y - A.y; if (mir) y = -y; return {m, x: P.x + x * c - y * sn, y: P.y + x * sn + y * c}; });
        // Rest möglichst senkrecht und ohne Zusammenstoß mit schon gesetzten Atomen
        let cx = 0; Q.forEach(q => { cx += q.x; }); cx /= Q.length;
        let sc = Math.abs(cx - P.x) * 0.5;
        for (const q of Q) for (const [id, p] of np) { const d = Math.hypot(q.x - p.x, q.y - p.y); if (d < 0.9) sc += (0.9 - d) * 10; }
        if (sc < bs) { bs = sc; bestP = Q; }
      }
      bestP.forEach(q => np.set(q.m, {x: q.x, y: q.y}));
    });
  }
  let nMin = Infinity, ny = 0, cnt = 0;
  for (const [id, p] of np) { nMin = Math.min(nMin, p.x); ny += p.y; cnt++; }
  ny /= cnt;
  for (const [id, p] of np) pos[id] = {x: p.x - nMin + minX, y: p.y - ny + sy};
}
const SHORT_BOND = 36 / 46;   // kurz (36 px) im Verhältnis zu mittel (46 px)
// Funktionelle Gruppen: zusammenhängende Heteroatome (samt H) außerhalb von Ringen, die über genau eine Bindung am
// Kohlenstoffgerüst bzw. an einem Ring hängen. Heteroatome zwischen zwei Gerüstatomen (Ester-O, Peptid-N, Ring-N) gehören zum Gerüst.
function functionalGroups(g, rings) {
  const n = g.atoms.length, ringSet = new Set((rings || findRings(g)).flat());
  if (!g.atoms.some(a => a.el === 'C')) return [];
  const skel = i => g.atoms[i].el === 'C' || ringSet.has(i);
  const seen = new Set(), out = [];
  for (let i = 0; i < n; i++) {
    if (skel(i) || seen.has(i)) continue;
    const comp = [i], att = []; seen.add(i);
    for (let k = 0; k < comp.length; k++) for (const v of g.nb(comp[k])) {
      if (skel(v)) att.push(v);
      else if (!seen.has(v)) { seen.add(v); comp.push(v); }
    }
    if (att.length === 1) out.push({attach: att[0], members: comp});
  }
  return out;
}
// Überlappungen auflösen: Seitenäste an Einfachbindungen (nicht im Ring) spiegeln oder drehen, solange es besser wird
function unclash(g, pos, rings) {
  const n = g.atoms.length; if (n < 4) return;
  const nb = Array.from({length: n}, (_, i) => g.nb(i));
  const bonded = new Set(); g.bonds.forEach(b => { bonded.add(b.a * n + b.b); bonded.add(b.b * n + b.a); });
  const ringBond = (a, b) => rings.some(r => { const i = r.indexOf(a), j = r.indexOf(b); return i >= 0 && j >= 0 && (Math.abs(i - j) === 1 || Math.abs(i - j) === r.length - 1); });
  const score = P => { let s = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { if (bonded.has(i * n + j)) continue; const d = Math.hypot(P[i].x - P[j].x, P[i].y - P[j].y); if (d < 0.8) s += (0.8 - d) * (0.8 - d) + (d < 0.3 ? 1 : 0); } return s; };
  const side = (a, b) => { const seen = new Set([a, b]), st = [b]; while (st.length) { const u = st.pop(); for (const v of nb[u]) if (!seen.has(v)) { seen.add(v); st.push(v); } } seen.delete(a); return seen; };
  const cands = [];
  for (const bd of g.bonds) {
    if (bd.order !== 1 || ringBond(bd.a, bd.b)) continue;
    for (const [a, b] of [[bd.a, bd.b], [bd.b, bd.a]]) { const S = side(a, b); if (S.has(a) || S.size > n / 2 + 1) continue; cands.push({a, b, S: [...S]}); }
  }
  let cur = score(pos);
  for (let pass = 0; pass < 40 && cur > 1e-9; pass++) {
    let best = null;
    for (const c of cands) {
      const A = pos[c.a], B = pos[c.b];
      const ux = B.x - A.x, uy = B.y - A.y, L2 = ux * ux + uy * uy || 1;
      const tfs = [p => { const t = ((p.x - A.x) * ux + (p.y - A.y) * uy) / L2, fx = A.x + ux * t, fy = A.y + uy * t; return {x: 2 * fx - p.x, y: 2 * fy - p.y}; }];
      for (const ang of [Math.PI / 3, -Math.PI / 3, 2 * Math.PI / 3, -2 * Math.PI / 3]) { const co = Math.cos(ang), si = Math.sin(ang); tfs.push(p => ({x: A.x + (p.x - A.x) * co - (p.y - A.y) * si, y: A.y + (p.x - A.x) * si + (p.y - A.y) * co})); }
      for (const tf of tfs) {
        const P = pos.slice(); for (const i of c.S) P[i] = tf(pos[i]);
        const sc = score(P);
        if (sc < cur - 1e-6 && (!best || sc < best.sc)) best = {sc, P};
      }
    }
    if (!best) break;
    for (let i = 0; i < n; i++) pos[i] = best.P[i];
    cur = best.sc;
  }
}

/* ======================= Namen ======================= */
const STEMS = {meth: 1, eth: 2, prop: 3, but: 4, pent: 5, hex: 6, hept: 7, oct: 8, non: 9, dec: 10, undec: 11, dodec: 12,
  tridec: 13, tetradec: 14, pentadec: 15, hexadec: 16, heptadec: 17, octadec: 18, nonadec: 19, icos: 20, eicos: 20};
const STEMKEYS = Object.keys(STEMS).sort((a, b) => b.length - a.length);
const MULT = {di: 2, tri: 3, tetra: 4, penta: 5, hexa: 6, bis: 2, tris: 3, tetrakis: 4};

// Trivialnamen (deutsch + englisch) → SMILES
const TRIVIAL = {
  'wasser': 'O', 'water': 'O', 'ammoniak': 'N', 'ammonia': 'N', 'kohlenstoffdioxid': 'O=C=O', 'kohlendioxid': 'O=C=O', 'carbon dioxide': 'O=C=O',
  'kohlenstoffmonoxid': '[C-]#[O+]', 'kohlenmonoxid': '[C-]#[O+]', 'carbon monoxide': '[C-]#[O+]', 'wasserstoffperoxid': 'OO', 'hydrogen peroxide': 'OO',
  'schwefelsäure': 'OS(=O)(=O)O', 'sulfuric acid': 'OS(=O)(=O)O', 'salpetersäure': 'O[N+](=O)[O-]', 'nitric acid': 'O[N+](=O)[O-]',
  'phosphorsäure': 'OP(=O)(O)O', 'phosphoric acid': 'OP(=O)(O)O', 'kohlensäure': 'OC(=O)O', 'carbonic acid': 'OC(=O)O',
  'salzsäure': 'Cl', 'chlorwasserstoff': 'Cl', 'hydrogen chloride': 'Cl', 'blausäure': 'C#N', 'cyanwasserstoff': 'C#N', 'hydrogen cyanide': 'C#N',
  'schwefelwasserstoff': 'S', 'hydrogen sulfide': 'S', 'ozon': '[O-][O+]=O', 'ozone': '[O-][O+]=O', 'harnstoff': 'NC(N)=O', 'urea': 'NC(N)=O',
  'essigsäure': 'CC(=O)O', 'acetic acid': 'CC(=O)O', 'ameisensäure': 'OC=O', 'formic acid': 'OC=O', 'propionsäure': 'CCC(=O)O', 'propionic acid': 'CCC(=O)O',
  'buttersäure': 'CCCC(=O)O', 'butyric acid': 'CCCC(=O)O', 'valeriansäure': 'CCCCC(=O)O', 'oxalsäure': 'OC(=O)C(=O)O', 'oxalic acid': 'OC(=O)C(=O)O',
  'malonsäure': 'OC(=O)CC(=O)O', 'bernsteinsäure': 'OC(=O)CCC(=O)O', 'succinic acid': 'OC(=O)CCC(=O)O', 'milchsäure': 'CC(O)C(=O)O', 'lactic acid': 'CC(O)C(=O)O',
  'brenztraubensäure': 'CC(=O)C(=O)O', 'pyruvic acid': 'CC(=O)C(=O)O', 'citronensäure': 'OC(=O)CC(O)(CC(=O)O)C(=O)O', 'zitronensäure': 'OC(=O)CC(O)(CC(=O)O)C(=O)O', 'citric acid': 'OC(=O)CC(O)(CC(=O)O)C(=O)O',
  'weinsäure': 'OC(C(O)C(=O)O)C(=O)O', 'tartaric acid': 'OC(C(O)C(=O)O)C(=O)O', 'äpfelsäure': 'OC(CC(=O)O)C(=O)O', 'malic acid': 'OC(CC(=O)O)C(=O)O',
  'fumarsäure': 'OC(=O)/C=C/C(=O)O', 'maleinsäure': 'OC(=O)/C=C\\C(=O)O', 'acrylsäure': 'C=CC(=O)O', 'acrylic acid': 'C=CC(=O)O',
  'stearinsäure': 'CCCCCCCCCCCCCCCCCC(=O)O', 'stearic acid': 'CCCCCCCCCCCCCCCCCC(=O)O', 'palmitinsäure': 'CCCCCCCCCCCCCCCC(=O)O', 'palmitic acid': 'CCCCCCCCCCCCCCCC(=O)O',
  'ölsäure': 'CCCCCCCC/C=C\\CCCCCCCC(=O)O', 'oleic acid': 'CCCCCCCC/C=C\\CCCCCCCC(=O)O',
  'aceton': 'CC(C)=O', 'acetone': 'CC(C)=O', 'formaldehyd': 'C=O', 'formaldehyde': 'C=O', 'acetaldehyd': 'CC=O', 'acetaldehyde': 'CC=O',
  'glycerin': 'OCC(O)CO', 'glycerol': 'OCC(O)CO', 'glykol': 'OCCO', 'ethylenglykol': 'OCCO', 'ethylene glycol': 'OCCO',
  'ethylen': 'C=C', 'ethylene': 'C=C', 'acetylen': 'C#C', 'acetylene': 'C#C', 'propylen': 'CC=C', 'propylene': 'CC=C', 'isopren': 'C=CC(C)=C', 'isoprene': 'C=CC(C)=C',
  'chloroform': 'ClC(Cl)Cl', 'tetrachlorkohlenstoff': 'ClC(Cl)(Cl)Cl', 'carbon tetrachloride': 'ClC(Cl)(Cl)Cl',
  'diethylether': 'CCOCC', 'diethyl ether': 'CCOCC', 'ether': 'CCOCC', 'dimethylsulfoxid': 'CS(C)=O', 'dimethyl sulfoxide': 'CS(C)=O', 'dmso': 'CS(C)=O',
  'ethylacetat': 'CCOC(C)=O', 'essigsäureethylester': 'CCOC(C)=O', 'ethyl acetate': 'CCOC(C)=O',
  'styrol': 'C=Cc1ccccc1', 'styrene': 'C=Cc1ccccc1', 'naphthalin': 'c1ccc2ccccc2c1', 'naphthalene': 'c1ccc2ccccc2c1',
  'anthracen': 'c1ccc2cc3ccccc3cc2c1', 'anthracene': 'c1ccc2cc3ccccc3cc2c1', 'phenanthren': 'c1ccc2c(c1)ccc1ccccc12', 'phenanthrene': 'c1ccc2c(c1)ccc1ccccc12',
  'biphenyl': 'c1ccc(cc1)-c1ccccc1', 'indol': 'c1ccc2[nH]ccc2c1', 'indole': 'c1ccc2[nH]ccc2c1', 'imidazol': 'c1c[nH]cn1', 'imidazole': 'c1c[nH]cn1',
  'pyrimidin': 'c1cncnc1', 'pyrimidine': 'c1cncnc1', 'purin': 'c1ncc2[nH]cnc2n1', 'purine': 'c1ncc2[nH]cnc2n1', 'chinolin': 'c1ccc2ncccc2c1', 'quinoline': 'c1ccc2ncccc2c1',
  'koffein': 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'coffein': 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'caffeine': 'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
  'aspirin': 'CC(=O)Oc1ccccc1C(=O)O', 'acetylsalicylsäure': 'CC(=O)Oc1ccccc1C(=O)O', 'acetylsalicylic acid': 'CC(=O)Oc1ccccc1C(=O)O',
  'salicylsäure': 'OC(=O)c1ccccc1O', 'salicylic acid': 'OC(=O)c1ccccc1O', 'paracetamol': 'CC(=O)Nc1ccc(O)cc1', 'acetaminophen': 'CC(=O)Nc1ccc(O)cc1',
  'ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'vanillin': 'COc1cc(C=O)ccc1O', 'nikotin': 'CN1CCCC1c1cccnc1', 'nicotine': 'CN1CCCC1c1cccnc1',
  'adrenalin': 'CNCC(O)c1ccc(O)c(O)c1', 'adrenaline': 'CNCC(O)c1ccc(O)c(O)c1', 'dopamin': 'NCCc1ccc(O)c(O)c1', 'dopamine': 'NCCc1ccc(O)c(O)c1',
  'serotonin': 'NCCc1c[nH]c2ccc(O)cc12', 'histamin': 'NCCc1c[nH]cn1', 'histamine': 'NCCc1c[nH]cn1', 'capsaicin': 'COc1cc(CNC(=O)CCCCC=CC(C)C)ccc1O',
  'menthol': 'CC(C)C1CCC(C)CC1O', 'kampfer': 'CC1(C)C2CCC1(C)C(=O)C2', 'campher': 'CC1(C)C2CCC1(C)C(=O)C2', 'camphor': 'CC1(C)C2CCC1(C)C(=O)C2',
  'glucose': 'OCC1OC(O)C(O)C(O)C1O', 'fructose': 'OCC1(O)OCC(O)C(O)C1O', 'ribose': 'OCC1OC(O)C(O)C1O', 'saccharose': 'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O', 'sucrose': 'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O',
  'cholesterin': 'CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C', 'cholesterol': 'CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C',
  'glycin': 'NCC(=O)O', 'glycine': 'NCC(=O)O', 'alanin': 'CC(N)C(=O)O', 'alanine': 'CC(N)C(=O)O', 'valin': 'CC(C)C(N)C(=O)O', 'valine': 'CC(C)C(N)C(=O)O',
  'leucin': 'CC(C)CC(N)C(=O)O', 'leucine': 'CC(C)CC(N)C(=O)O', 'isoleucin': 'CCC(C)C(N)C(=O)O', 'isoleucine': 'CCC(C)C(N)C(=O)O',
  'serin': 'OCC(N)C(=O)O', 'serine': 'OCC(N)C(=O)O', 'threonin': 'CC(O)C(N)C(=O)O', 'threonine': 'CC(O)C(N)C(=O)O', 'cystein': 'SCC(N)C(=O)O', 'cysteine': 'SCC(N)C(=O)O',
  'methionin': 'CSCCC(N)C(=O)O', 'methionine': 'CSCCC(N)C(=O)O', 'asparaginsäure': 'OC(=O)CC(N)C(=O)O', 'aspartic acid': 'OC(=O)CC(N)C(=O)O',
  'glutaminsäure': 'OC(=O)CCC(N)C(=O)O', 'glutamic acid': 'OC(=O)CCC(N)C(=O)O', 'asparagin': 'NC(=O)CC(N)C(=O)O', 'asparagine': 'NC(=O)CC(N)C(=O)O',
  'glutamin': 'NC(=O)CCC(N)C(=O)O', 'glutamine': 'NC(=O)CCC(N)C(=O)O', 'lysin': 'NCCCCC(N)C(=O)O', 'lysine': 'NCCCCC(N)C(=O)O',
  'arginin': 'NC(=N)NCCCC(N)C(=O)O', 'arginine': 'NC(=N)NCCCC(N)C(=O)O', 'phenylalanin': 'NC(Cc1ccccc1)C(=O)O', 'phenylalanine': 'NC(Cc1ccccc1)C(=O)O',
  'tyrosin': 'NC(Cc1ccc(O)cc1)C(=O)O', 'tyrosine': 'NC(Cc1ccc(O)cc1)C(=O)O', 'tryptophan': 'NC(Cc1c[nH]c2ccccc12)C(=O)O',
  'histidin': 'NC(Cc1c[nH]cn1)C(=O)O', 'histidine': 'NC(Cc1c[nH]cn1)C(=O)O', 'prolin': 'OC(=O)C1CCCN1', 'proline': 'OC(=O)C1CCCN1',
  'adenin': 'Nc1ncnc2[nH]cnc12', 'adenine': 'Nc1ncnc2[nH]cnc12', 'guanin': 'Nc1nc2[nH]cnc2c(=O)[nH]1', 'guanine': 'Nc1nc2[nH]cnc2c(=O)[nH]1',
  'cytosin': 'Nc1cc[nH]c(=O)n1', 'cytosine': 'Nc1cc[nH]c(=O)n1', 'thymin': 'Cc1c[nH]c(=O)[nH]c1=O', 'thymine': 'Cc1c[nH]c(=O)[nH]c1=O', 'uracil': 'O=c1cc[nH]c(=O)[nH]1',
  'harnsäure': 'O=c1[nH]c(=O)c2[nH]c(=O)[nH]c2[nH]1', 'uric acid': 'O=c1[nH]c(=O)c2[nH]c(=O)[nH]c2[nH]1', 'tnt': 'Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]',
  'trinitrotoluol': 'Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]', 'pikrinsäure': 'Oc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]',
  'hydrochinon': 'Oc1ccc(O)cc1', 'hydroquinone': 'Oc1ccc(O)cc1', 'resorcin': 'Oc1cccc(O)c1', 'brenzcatechin': 'Oc1ccccc1O', 'catechol': 'Oc1ccccc1O',
  'benzochinon': 'O=C1C=CC(=O)C=C1', 'benzoquinone': 'O=C1C=CC(=O)C=C1', 'phthalsäure': 'OC(=O)c1ccccc1C(=O)O', 'terephthalsäure': 'OC(=O)c1ccc(cc1)C(=O)O',
  'xylol': 'Cc1ccccc1C', 'xylene': 'Cc1ccccc1C', 'cumol': 'CC(C)c1ccccc1', 'cumene': 'CC(C)c1ccccc1', 'acetophenon': 'CC(=O)c1ccccc1', 'acetophenone': 'CC(=O)c1ccccc1',
  'benzophenon': 'O=C(c1ccccc1)c1ccccc1', 'anisol': 'COc1ccccc1', 'anisole': 'COc1ccccc1', 'benzylalkohol': 'OCc1ccccc1', 'benzyl alcohol': 'OCc1ccccc1',
  'zimtaldehyd': 'O=C/C=C/c1ccccc1', 'cinnamaldehyde': 'O=C/C=C/c1ccccc1', 'methanol': 'CO', 'ethanol': 'CCO', 'isopropanol': 'CC(C)O', 'isopropyl alcohol': 'CC(C)O',
  'tetrahydrofuran': 'C1CCOC1', 'thf': 'C1CCOC1', 'dioxan': 'C1COCCO1', 'dioxane': 'C1COCCO1', 'morpholin': 'C1COCCN1', 'morpholine': 'C1COCCN1',
  'piperidin': 'C1CCNCC1', 'piperidine': 'C1CCNCC1', 'pyrrolidin': 'C1CCNC1', 'pyrrolidine': 'C1CCNC1', 'piperazin': 'C1CNCCN1', 'piperazine': 'C1CNCCN1',
  'acetonitril': 'CC#N', 'acetonitrile': 'CC#N', 'dimethylformamid': 'CN(C)C=O', 'dimethylformamide': 'CN(C)C=O', 'dmf': 'CN(C)C=O',
  'acetamid': 'CC(N)=O', 'acetamide': 'CC(N)=O', 'formamid': 'NC=O', 'formamide': 'NC=O', 'essigsäureanhydrid': 'CC(=O)OC(C)=O', 'acetic anhydride': 'CC(=O)OC(C)=O',
  'acetylchlorid': 'CC(Cl)=O', 'acetyl chloride': 'CC(Cl)=O', 'phosgen': 'ClC(Cl)=O', 'phosgene': 'ClC(Cl)=O', 'methylamin': 'CN', 'methylamine': 'CN', 'mandelsäure': 'OC(C(=O)O)c1ccccc1', 'mandelic acid': 'OC(C(=O)O)c1ccccc1', 'phenylglycin': 'NC(C(=O)O)c1ccccc1', 'phenylglycine': 'NC(C(=O)O)c1ccccc1',
};

// natürliche Konfiguration: L-Aminosäuren, D-Zucker (offenkettig); cip: Atomindex im SMILES → R/S
const ST = (s, cip) => ({s, cip});
const STEREO_TRIVIAL = {
  'alanin|alanine': ST('CC(N)C(=O)O', {1: 'S'}), 'valin|valine': ST('CC(C)C(N)C(=O)O', {3: 'S'}),
  'leucin|leucine': ST('CC(C)CC(N)C(=O)O', {4: 'S'}), 'isoleucin|isoleucine': ST('CCC(C)C(N)C(=O)O', {2: 'S', 4: 'S'}),
  'serin|serine': ST('OCC(N)C(=O)O', {2: 'S'}), 'threonin|threonine': ST('CC(O)C(N)C(=O)O', {1: 'R', 3: 'S'}),
  'cystein|cysteine': ST('SCC(N)C(=O)O', {2: 'R'}), 'methionin|methionine': ST('CSCCC(N)C(=O)O', {4: 'S'}),
  'asparaginsäure|aspartic acid': ST('OC(=O)CC(N)C(=O)O', {4: 'S'}), 'glutaminsäure|glutamic acid': ST('OC(=O)CCC(N)C(=O)O', {5: 'S'}),
  'asparagin|asparagine': ST('NC(=O)CC(N)C(=O)O', {4: 'S'}), 'glutamin|glutamine': ST('NC(=O)CCC(N)C(=O)O', {5: 'S'}),
  'lysin|lysine': ST('NCCCCC(N)C(=O)O', {5: 'S'}), 'arginin|arginine': ST('NC(=N)NCCCC(N)C(=O)O', {7: 'S'}),
  'phenylalanin|phenylalanine': ST('NC(Cc1ccccc1)C(=O)O', {1: 'S'}), 'tyrosin|tyrosine': ST('NC(Cc1ccc(O)cc1)C(=O)O', {1: 'S'}),
  'tryptophan': ST('NC(Cc1c[nH]c2ccccc12)C(=O)O', {1: 'S'}), 'histidin|histidine': ST('NC(Cc1c[nH]cn1)C(=O)O', {1: 'S'}),
  'prolin|proline': ST('OC(=O)C1CCCN1', {3: 'S'}),
  'glucose': ST('O=CC(O)C(O)C(O)C(O)CO', {2: 'R', 4: 'S', 6: 'R', 8: 'R'}),
  'galactose': ST('O=CC(O)C(O)C(O)C(O)CO', {2: 'R', 4: 'S', 6: 'S', 8: 'R'}),
  'mannose': ST('O=CC(O)C(O)C(O)C(O)CO', {2: 'S', 4: 'S', 6: 'R', 8: 'R'}),
  'fructose': ST('OCC(=O)C(O)C(O)C(O)CO', {4: 'S', 6: 'R', 8: 'R'}),
  'ribose': ST('O=CC(O)C(O)C(O)CO', {2: 'R', 4: 'R', 6: 'R'}),
  'desoxyribose|deoxyribose|2-desoxyribose|2-deoxyribose': ST('O=CCC(O)C(O)CO', {3: 'S', 5: 'R'}),
  'xylose': ST('O=CC(O)C(O)C(O)CO', {2: 'R', 4: 'S', 6: 'R'}), 'arabinose': ST('O=CC(O)C(O)C(O)CO', {2: 'S', 4: 'R', 6: 'R'}),
  'glycerinaldehyd|glyceraldehyde|glycerinaldehyde': ST('O=CC(O)CO', {2: 'R'}),
  'erythrose': ST('O=CC(O)C(O)CO', {2: 'R', 4: 'R'}), 'threose': ST('O=CC(O)C(O)CO', {2: 'S', 4: 'R'}),
  'weinsäure|tartaric acid': ST('OC(C(O)C(=O)O)C(=O)O', {1: 'R', 2: 'R'}),
  'glucopyranose|glucose (ring)|glucose ring': ST('OCC1OC(O)C(O)C(O)C1O', {}),
};
for (const [keys, v] of Object.entries(STEREO_TRIVIAL)) for (const k of keys.split('|')) TRIVIAL[k] = v;
for (const k of Object.keys(TRIVIAL)) if (k.includes(' ')) TRIVIAL[k.replace(/\s+/g, '')] = TRIVIAL[k];
function trivialGraph(v) {
  if (typeof v === 'string') return parseSmiles(v);
  const g = parseSmiles(v.s);
  for (const [i, d] of Object.entries(v.cip)) g.atoms[+i].cip = d;
  return g;
}

// Ring-Stammsysteme: SMILES-Atomreihenfolge = IUPAC-Nummerierung
const RINGPARENTS = [
  ['benzol', 'c1ccccc1'], ['benzene', 'c1ccccc1'], ['pyridin', 'n1ccccc1'], ['pyridine', 'n1ccccc1'],
  ['furan', 'o1cccc1'], ['thiophen', 's1cccc1'], ['thiophene', 's1cccc1'], ['pyrrol', '[nH]1cccc1'], ['pyrrole', '[nH]1cccc1'],
  ['phenol', 'c1ccccc1', [['hydroxy', 1]]], ['anilin', 'c1ccccc1', [['amino', 1]]], ['aniline', 'c1ccccc1', [['amino', 1]]],
  ['benzolamin', 'c1ccccc1', [['amino', 1]]], ['benzenamine', 'c1ccccc1', [['amino', 1]]],
  ['benzoesäure', 'c1ccccc1', [['carboxy', 1]]], ['benzoicacid', 'c1ccccc1', [['carboxy', 1]]],
  ['benzolcarbonsäure', 'c1ccccc1', [['carboxy', 1]]], ['benzenecarboxylicacid', 'c1ccccc1', [['carboxy', 1]]],
  ['benzaldehyd', 'c1ccccc1', [['formyl', 1]]], ['benzaldehyde', 'c1ccccc1', [['formyl', 1]]],
  ['benzonitril', 'c1ccccc1', [['cyano', 1]]], ['benzonitrile', 'c1ccccc1', [['cyano', 1]]],
  ['toluol', 'c1ccccc1', [['methyl', 1]]], ['toluene', 'c1ccccc1', [['methyl', 1]]],
  ['naphthalin', 'c1cccc2ccccc12'], ['naphthalene', 'c1cccc2ccccc12'],
].sort((a, b) => b[0].length - a[0].length);
// Naphthalin: SMILES-Reihenfolge c1(1)c(2)c(3)c(4)c2(4a)c(5)c(6)c(7)c(8)c12(8a) → Locant → Atomindex
const NAPHTHMAP = {1: 0, 2: 1, 3: 2, 4: 3, 5: 5, 6: 6, 7: 7, 8: 8};

// ---- Substituenten ----
function chainGraph(n, cyclo) { const g = new G(); for (let i = 0; i < n; i++) { g.add('C'); if (i) g.bond(i - 1, i, 1); } if (cyclo && n > 2) g.bond(n - 1, 0, 1); return g; }
const S = (smiles, attach = 0, order = 1, charge) => ({build: () => ({g: parseSmiles(smiles), attach, order})});
const SUBST = {
  'fluoro': S('F'), 'fluor': S('F'), 'chloro': S('Cl'), 'chlor': S('Cl'), 'bromo': S('Br'), 'brom': S('Br'), 'iodo': S('I'), 'iod': S('I'), 'jod': S('I'),
  'hydroxy': S('O'), 'amino': S('N'), 'nitro': S('[N+](=O)[O-]'), 'nitroso': S('N=O'), 'oxo': S('O', 0, 2), 'cyano': S('C#N'), 'cyan': S('C#N'),
  'sulfanyl': S('S'), 'mercapto': S('S'), 'carboxy': S('C(=O)O'), 'formyl': S('C=O'), 'carbamoyl': S('C(N)=O'), 'acetyl': S('C(C)=O'),
  'imino': S('N', 0, 2), 'sulfo': S('S(=O)(=O)O'), 'trifluormethyl': S('C(F)(F)F'), 'trifluoromethyl': S('C(F)(F)F'),
  'vinyl': S('C=C'), 'ethenyl': S('C=C'), 'ethinyl': S('C#C'), 'ethynyl': S('C#C'), 'allyl': S('CC=C'),
  'phenyl': S('c1ccccc1'), 'benzyl': S('Cc1ccccc1'), 'phenoxy': S('Oc1ccccc1'), 'benzoyl': S('C(=O)c1ccccc1'),
  'isopropyl': S('C(C)C'), 'isobutyl': S('CC(C)C'), 'secbutyl': S('C(C)CC'), 'tertbutyl': S('C(C)(C)C'), 'isopentyl': S('CCC(C)C'), 'neopentyl': S('CC(C)(C)C'),
  'methyliden': S('C', 0, 2), 'methylidene': S('C', 0, 2), 'methylen': S('C', 0, 2), 'methylene': S('C', 0, 2),
};
// Alkoxy- und Alkylthio-Gruppen
for (const [st, n] of Object.entries(STEMS)) {
  if (n <= 6) SUBST[st + 'oxy'] = {build: () => { const g = chainGraph(n); const o = g.add('O'); g.bond(o, 0, 1); return {g, attach: o, order: 1}; }};
  SUBST[st + 'yloxy'] = SUBST[st + 'oxy'] || {build: () => { const g = chainGraph(n); const o = g.add('O'); g.bond(o, 0, 1); return {g, attach: o, order: 1}; }};
}
const SUBSTKEYS = Object.keys(SUBST).sort((a, b) => b.length - a.length);

// ---- Stamm (Kette/Ring) mit Mehrfachbindungen und Hauptgruppe ----
const SUFFIXES = [
  ['carbonsäure', 'cooh_c'], ['carboxylicacid', 'cooh_c'], ['carbaldehyd', 'cho_c'], ['carbaldehyde', 'cho_c'], ['carbonitril', 'cn_c'], ['carbonitrile', 'cn_c'],
  ['säure', 'acid'], ['saeure', 'acid'], ['oicacid', 'acid'], ['oate', 'oat'], ['oat', 'oat'], ['amide', 'amid'], ['amid', 'amid'],
  ['nitrile', 'nitril'], ['nitril', 'nitril'], ['amine', 'amin'], ['amin', 'amin'], ['thiol', 'thiol'], ['one', 'on'], ['on', 'on'],
  ['al', 'al'], ['ol', 'ol'],
];
function parseLocs(str) { return str.split(',').map(x => x === 'n' ? 'n' : +x); }

// liest Stamm ab Position 0 von t; sub=true → Substituent (Endung yl)
function parseParent(t, sub) {
  if (!sub) for (const [key, smi, extra] of RINGPARENTS) {
    if (t === key) return {len: key.length, spec: {ringParent: key, smiles: smi, extra: extra || []}};
  }
  let i = 0, cyclo = false;
  if (t.startsWith('cyclo')) { cyclo = true; i = 5; }
  const stem = STEMKEYS.find(k => t.startsWith(k, i));
  if (!stem) return null;
  const n = STEMS[stem]; i += stem.length;
  const spec = {cyclo, n, unsat: [], suffix: null};
  const rest = () => t.slice(i);
  // Mehrfachbindungen
  let m, sawUnsat = false;
  while ((m = rest().match(/^a?(?:-?([\dn,]+)-)?(di|tri|tetra)?(en|yn|in)(?:e(?=$|-|[bcdfghjklmnpqrstvwxz]))?/))) {
    // "in" darf nicht Teil von "amin"/"amine" sein → bereits durch Muster ausgeschlossen
    const locs = m[1] ? parseLocs(m[1]) : null;
    const count = m[2] ? MULT[m[2]] : 1;
    if (locs && locs.length !== count) return null;
    spec.unsat.push({locs, count, order: m[3] === 'en' ? 2 : 3});
    i += m[0].length; sawUnsat = true;
  }
  if (!sawUnsat) {
    if (rest().startsWith('an')) { i += 2; if (/^e(?=$|-|[bcdfghjklmnpqrstvwxz])/.test(rest()) ) i++; }
    else if (!(sub && rest().startsWith('yl'))) return null;
  }
  if (sub) {
    m = rest().match(/^(?:-?(\d+)-)?yl/);
    if (!m) return null;
    spec.attach = m[1] ? +m[1] : 1; i += m[0].length;
    return {len: i, spec};
  }
  if (rest() === '') return {len: i, spec};
  m = rest().match(/^-?(?:([\dn,]+)-)?(di|tri|tetra?|penta?|hexa?)?/);
  const MV = {di: 2, tri: 3, tetr: 4, tetra: 4, pent: 5, penta: 5, hex: 6, hexa: 6};
  const locs = m[1] ? parseLocs(m[1]) : null, count = m[2] ? MV[m[2]] : 1;
  let r2 = rest().slice(m[0].length);
  const suf = SUFFIXES.find(([k]) => r2 === k);
  if (!suf) return null;
  if (locs && locs.length !== count) return null;
  spec.suffix = {type: suf[1], locs, count};
  return {len: t.length, spec};
}

// Substituenten-Liste lesen: gibt [{locs, count, sub}] und Rest zurück
function parseSubstAt(t) {
  // komplexer Substituent in Klammern
  if (t[0] === '(' || t[0] === '[') {
    const close = t[0] === '(' ? ')' : ']';
    let depth = 0, j = 0;
    for (; j < t.length; j++) { if (t[j] === t[0]) depth++; else if (t[j] === close && --depth === 0) break; }
    if (j >= t.length) return null;
    const inner = t.slice(1, j);
    const sb = SUBST[inner] ? SUBST[inner] : complexSubst(inner);
    if (!sb) return null;
    return {len: j + 1, sub: sb};
  }
  for (const k of SUBSTKEYS) if (t.startsWith(k)) return {len: k.length, sub: SUBST[k], key: k};
  // Alkyl / Cycloalkyl / Alkenyl
  const pp = parseParent(t, true);
  if (pp) return {len: pp.len, sub: {build: () => buildSubParent(pp.spec), alkyl: true}};
  return null;
}
function complexSubst(inner) {
  const r = parseName(inner, true);
  return r ? {build: () => r(), alkyl: true} : null;
}
function buildSubParent(spec, prefixes = []) {
  const g = new G();
  const C = buildSkeleton(g, spec);
  applyPrefixes(g, C, prefixes, spec, null);
  return {g, attach: C[spec.attach - 1], order: 1};
}

function readPrefixes(t, sub) {
  const list = []; let i = 0;
  while (i < t.length) {
    const r = t.slice(i);
    let m = r.match(/^-?([\dn,]+)-/);
    let locs = null;
    let j = 0;
    if (m) { locs = parseLocs(m[1]); j = m[0].length; }
    else if (r[0] === '-') j = 1;
    const r2 = r.slice(j);
    let count = 1;
    let mm = r2.match(/^(tetrakis|tris|bis|tetra|tri|di|penta|hexa)(?=[a-z(\[])/);
    let k = 0;
    if (mm) {
      // "di" nur als Vervielfacher, wenn dahinter ein Substituent folgt
      const test = parseSubstAt(r2.slice(mm[0].length));
      if (test) { count = MULT[mm[1]]; k = mm[0].length; } else mm = null;
    }
    const ps = parseSubstAt(r2.slice(k));
    if (!ps) return {list, rest: r, pos: i};
    if (locs && locs.length !== count) {
      if (count === 1 && locs.length > 1) count = locs.length; else return {list, rest: r, pos: i};
    }
    list.push({locs, count, sub: ps.sub});
    i += j + k + ps.len;
    // Prüfen, ob der Rest bereits ein Stamm ist (dann aufhören)
    const rest = t.slice(i).replace(/^-/, '');
    if (rest) { const pp = parseParent(rest, !!sub); if (pp && pp.len === rest.length) return {list, rest, pos: i, parentFound: true}; }
  }
  return {list, rest: '', pos: i};
}

function buildSkeleton(g, spec) {
  if (spec.smiles) {
    const off = g.merge(parseSmiles(spec.smiles));
    const C = []; const n = g.atoms.length - off; for (let k = 0; k < n; k++) C.push(off + k);
    return C;
  }
  const C = []; for (let k = 0; k < spec.n; k++) { C.push(g.add('C')); if (k) g.bond(C[k - 1], C[k], 1); }
  if (spec.cyclo) { if (spec.n < 3) throw new Error('Ring zu klein'); g.bond(C[spec.n - 1], C[0], 1); }
  // Mehrfachbindungen
  for (const u of spec.unsat) {
    // ohne Lokanten: bei Endgruppe an C1 (Propenal, Propensäure, Propennitril) beginnt die Mehrfachbindung an C2
    const st = !spec.cyclo && spec.suffix && ['al', 'acid', 'nitril', 'oat', 'amid'].includes(spec.suffix.type) && (!spec.suffix.locs || spec.suffix.locs[0] === 1) ? 2 : 1;
    const locs = u.locs || Array.from({length: u.count}, (_, k) => st + 2 * k);
    for (const L of locs) {
      const a = C[L - 1], b = spec.cyclo ? C[L % spec.n] : C[L];
      if (a === undefined || b === undefined) throw new Error('Mehrfachbindung an Position ' + L + ' passt nicht zur Kettenlänge ' + spec.n + '.');
      g.bond(a, b, u.order);
    }
  }
  return C;
}

function locAtom(C, spec, L) {
  if (spec.ringParent && spec.ringParent.startsWith('naphthal')) { const k = NAPHTHMAP[L]; if (k === undefined) throw new Error('Position ' + L); return C[k]; }
  const a = C[L - 1]; if (a === undefined) throw new Error('Position ' + L + ' gibt es im Stammsystem nicht (nur 1–' + C.length + ').'); return a;
}

function applyPrefixes(g, C, prefixes, spec, nAtoms) {
  const chainLen = C.length;
  for (const p of prefixes) {
    let locs = p.locs;
    if (!locs) {
      const isRing = spec.cyclo || spec.ringParent;
      if (isRing) locs = Array.from({length: p.count}, (_, k) => 1 + k);
      else if (p.sub.alkyl && chainLen >= 3) locs = Array(p.count).fill(2);
      else {
        // ohne Locants: der Reihe nach auf Atome mit freien Valenzen verteilen
        const ord = p.sub.build().order;
        const free = C.map(c => 4 - g.bondsOf(c).reduce((s, b) => s + b.order, 0));
        locs = [];
        for (let k = 0; k < p.count; k++) {
          const idx = free.findIndex(f => f >= ord);
          if (idx < 0) throw new Error('Zu viele Substituenten für diese Kette');
          free[idx] -= ord; locs.push(idx + 1);
        }
      }
    }
    let nIdx = 0;
    for (const L of locs) {
      let target;
      if (L === 'n') { if (!nAtoms || !nAtoms.length) throw new Error('N-Locant ohne Amin'); target = nAtoms[Math.min(nIdx++, nAtoms.length - 1)]; }
      else target = locAtom(C, spec, L);
      const s = p.sub.build();
      const off = g.merge(s.g);
      g.bond(target, off + s.attach, s.order);
    }
  }
}

function applySuffix(g, C, spec, esterAlkyl) {
  const sf = spec.suffix; const nAtoms = [];
  if (!sf) return nAtoms;
  const n = C.length, ring = spec.cyclo || spec.ringParent;
  const endType = ['al', 'acid', 'nitril', 'oat', 'amid'].includes(sf.type);
  let locs = sf.locs;
  if (!locs) {
    if (ring) locs = Array.from({length: sf.count}, (_, k) => 1 + k);
    else if (endType) locs = sf.count === 1 ? [1] : [1, n];
    else if (sf.type === 'on') locs = sf.count === 1 ? [n >= 3 ? 2 : 1] : [2, n - 1];
    else locs = Array.from({length: sf.count}, (_, k) => k === 0 ? 1 : (sf.count === 2 ? n : 1 + k));
  }
  for (const L of locs) {
    const c = locAtom(C, spec, L);
    switch (sf.type) {
      case 'ol': g.bond(c, g.add('O'), 1); break;
      case 'thiol': g.bond(c, g.add('S'), 1); break;
      case 'amin': { const N = g.add('N'); g.bond(c, N, 1); nAtoms.push(N); break; }
      case 'on': g.bond(c, g.add('O'), 2); break;
      case 'al': g.bond(c, g.add('O'), 2); break;
      case 'acid': g.bond(c, g.add('O'), 2); g.bond(c, g.add('O'), 1); break;
      case 'amid': { g.bond(c, g.add('O'), 2); const N = g.add('N'); g.bond(c, N, 1); nAtoms.push(N); break; }
      case 'nitril': g.bond(c, g.add('N'), 3); break;
      case 'oat': {
        g.bond(c, g.add('O'), 2); const O = g.add('O'); g.bond(c, O, 1);
        if (esterAlkyl) { const s = esterAlkyl.build(); const off = g.merge(s.g); g.bond(O, off + s.attach, 1); }
        else g.atoms[O].charge = -1;
        break;
      }
      case 'cooh_c': { const k = g.add('C'); g.bond(c, k, 1); g.bond(k, g.add('O'), 2); g.bond(k, g.add('O'), 1); break; }
      case 'cho_c': { const k = g.add('C'); g.bond(c, k, 1); g.bond(k, g.add('O'), 2); break; }
      case 'cn_c': { const k = g.add('C'); g.bond(c, k, 1); g.bond(k, g.add('N'), 3); break; }
    }
  }
  return nAtoms;
}

// Hauptfunktion: Name → Builder (sub=true: Substituent, gibt {g, attach, order} zurück)
function parseName(t, sub) {
  let prefixes = [], rest = t;
  const whole = parseParent(t, !!sub);
  if (!whole || whole.len !== t.length) {
    const pre = readPrefixes(t, sub);
    prefixes = pre.list; rest = pre.rest.replace(/^-/, '');
    if (!prefixes.length) rest = t;
  }
  const pp = parseParent(rest, !!sub);
  if (!pp || pp.len !== rest.length) return null;
  const spec = pp.spec;
  return () => {
    const g = new G();
    const C = buildSkeleton(g, spec);
    if (!sub) g.mainChain = C;
    let esterAlkyl = null, pf = prefixes;
    if (!sub && spec.suffix && spec.suffix.type === 'oat' && pf.length && !pf[0].locs && pf[0].count === 1 && pf[0].sub.alkyl) {
      esterAlkyl = pf[0].sub; pf = pf.slice(1);
    }
    const extra = (spec.extra || []).map(([k, L]) => ({locs: [L], count: 1, sub: parseSubstAt(k).sub}));
    const nAtoms = applySuffix(g, C, spec, esterAlkyl);
    applyPrefixes(g, C, [...extra, ...pf], spec, nAtoms);
    if (sub) return {g, attach: C[(spec.attach || 1) - 1], order: 1};
    return g;
  };
}

// Amine und Ether aus Gruppennamen: "dimethylamin", "ethylmethylether"
function groupClass(t) {
  const m = t.match(/^(.*?)(amine|amin|ether)$/);
  if (!m || !m[1]) return null;
  const pre = readPrefixes(m[1]);
  if (pre.rest !== '' || !pre.list.length) return null;
  if (pre.list.some(p => p.locs)) return null;
  const total = pre.list.reduce((s, p) => s + p.count, 0);
  const isEther = m[2] === 'ether';
  if (isEther ? total !== 2 && total !== 1 : total > 3) return null;
  return () => {
    const g = new G(); const X = g.add(isEther ? 'O' : 'N');
    const items = pre.list.flatMap(p => Array(p.count).fill(p.sub));
    if (isEther && items.length === 1) items.push(items[0]);
    for (const s of items) { const b = s.build(); const off = g.merge(b.g); g.bond(X, off + b.attach, b.order); }
    return g;
  };
}

function normalizeName(raw) {
  let t = raw.trim().toLowerCase();
  t = t.replace(/[‐-―−]/g, '-').replace(/\s*,\s*/g, ',');
  t = t.replace(/\b(tert|sec|iso|neo)[-\s]+/g, '$1');
  t = t.replace(/(^|[\s(-])\([ezrs\d,]+\)-/g, '$1').replace(/(^|[\s(])(cis|trans|[ezrsdl])-(?=[a-z(\d])/g, '$1');
  t = t.replace(/\s+acid$/, 'acid').replace(/ic acid$/, 'icacid');
  return t;
}

function nameToGraph(raw) {
  let t = normalizeName(raw);
  if (TRIVIAL[t]) return {g: trivialGraph(TRIVIAL[t]), how: 'trivial'};
  const tNoSpace = t.replace(/\s+/g, '');
  if (TRIVIAL[tNoSpace]) return {g: trivialGraph(TRIVIAL[tNoSpace]), how: 'trivial'};
  // Englischer Ester: "ethyl ethanoate", "methyl benzoate"
  let m = t.match(/^(\S+yl)\s+(\S+?)(oate|ate)$/);
  if (m) {
    const alk = parseSubstAt(m[1]);
    let acidName = m[2] + (m[3] === 'oate' ? 'oicacid' : 'icacid');
    if (m[2] === 'benzo') acidName = 'benzoicacid';
    if (alk && alk.len === m[1].length) {
      const acid = tryName(acidName);
      if (acid) return {g: esterify(acid, alk.sub), how: 'name'};
    }
  }
  // Deutscher Ester: "Ethansäureethylester", "Essigsäuremethylester"
  m = tNoSpace.match(/^(.+säure)(.+?yl)ester$/);
  if (m) {
    const acid = tryName(m[1]); const alk = parseSubstAt(m[2]);
    if (acid && alk && alk.len === m[2].length) return {g: esterify(acid, alk.sub), how: 'name'};
  }
  t = tNoSpace;
  const g = tryName(t);
  if (g) return {g, how: 'name'};
  return null;
}
function tryName(t) {
  if (TRIVIAL[t]) return trivialGraph(TRIVIAL[t]);
  try {
    const b = parseName(t, false); if (b) return b();
    const gc = groupClass(t); if (gc) return gc();
  } catch (e) { throw e; }
  return null;
}
// Carbonsäure → Ester: Wasserstoff der OH-Gruppe durch Alkyl ersetzen
function esterify(g, alkSub) {
  const acidO = g.atoms.findIndex((a, i) => a.el === 'O' && g.bondsOf(i).length === 1 && g.bondsOf(i)[0].order === 1 &&
    (() => { const c = g.nb(i)[0]; return g.bondsOf(c).some(b => b.order === 2 && g.atoms[b.a === c ? b.b : b.a].el === 'O'); })());
  if (acidO < 0) throw new Error('Keine Säuregruppe gefunden');
  const s = alkSub.build(); const off = g.merge(s.g); g.bond(acidO, off + s.attach, 1);
  return g;
}

/* ======================= Stereochemie ======================= */
const ZNUM = {H: 1, B: 5, C: 6, N: 7, O: 8, F: 9, Na: 11, Si: 14, P: 15, S: 16, Cl: 17, K: 19, Br: 35, I: 53};
const VALS = {H: [1], B: [3], C: [4], Si: [4], N: [3], P: [3, 5], O: [2], S: [2, 4, 6], F: [1], Cl: [1, 3, 5, 7], Br: [1, 3, 5, 7], I: [1, 3, 5, 7]};
// implizite H-Atome – gleiche Regeln wie in der Zeichenfläche
function implicitHg(g, i) {
  const a = g.atoms[i];
  if (a.hx !== undefined && a.hx !== null) return a.hx;
  const list = VALS[a.el]; if (!list) return 0;
  const bs = g.bondsOf(i).reduce((s, b) => s + (b.order === 'ar' ? 1.5 : b.order), 0), ch = a.charge || 0;
  let v;
  if (a.el === 'C' || a.el === 'Si') v = 4 - Math.abs(ch);
  else if (a.el === 'B') v = 3 - ch;
  else if (a.el === 'H') v = 1 - Math.abs(ch);
  else { const opts = list.map(x => x + ch).filter(x => x >= 0); v = opts.find(x => x >= bs); if (v === undefined) v = opts[opts.length - 1] ?? 0; }
  return Math.max(0, Math.round(v - bs));
}
function hCounts(g) { return g.atoms.map((_, i) => implicitHg(g, i)); }

// CIP-Rangfolge (Regel 1, hierarchischer Digraph mit Duplikatatomen, sphärenweiser Vergleich)
function makeCip(g, hc) {
  const adj = g.atoms.map(() => []);
  g.bonds.forEach(b => { const o = b.order === 'ar' ? 1 : b.order; adj[b.a].push({n: b.b, o}); adj[b.b].push({n: b.a, o}); });
  const z = i => i < 0 ? 1 : (ZNUM[g.atoms[i].el] || 0);
  const ordOf = (a, b) => (adj[a].find(x => x.n === b) || {o: 1}).o;
  function childrenUnsorted(node) {
    const c = childrenBase(node);
    // für den flachen Vergleich: nach z, dann nach z-Summe der Nachfolger
    const w = n => n.dup || n.atom < 0 ? 0 : childrenBase(n).reduce((s, x) => s + x.z, 0);
    return c.sort((x, y) => (y.z - x.z) || (w(y) - w(x)));
  }
  function children(node) {
    const out = childrenBase(node);
    return out.sort((x, y) => (y.z - x.z) || -cmp(x, y, 4, true));
  }
  function childrenBase(node) {
    if (node.dup || node.atom < 0) return [];
    const a = node.atom, p = node.parent ? node.parent.atom : null, out = [];
    for (const {n, o} of adj[a]) {
      if (n === p) { for (let k = 1; k < o; k++) out.push({atom: n, z: z(n), dup: true}); continue; }
      let anc = node.parent, isAnc = false;
      while (anc) { if (anc.atom === n) { isAnc = true; break; } anc = anc.parent; }
      if (isAnc) { out.push({atom: n, z: z(n), dup: true}); continue; }
      out.push({atom: n, z: z(n), parent: node});
      for (let k = 1; k < o; k++) out.push({atom: n, z: z(n), dup: true});
    }
    for (let k = 0; k < hc[a]; k++) out.push({atom: -1, z: 1, dup: true});
    return out;
  }
  function cmp(n1, n2, maxDepth = 40, shallow = false) {
    let F1 = [[n1]], F2 = [[n2]];
    for (let depth = 0; depth < maxDepth; depth++) {
      const L = Math.max(F1.length, F2.length);
      for (let k = 0; k < L; k++) {
        const A = F1[k] || [], B = F2[k] || [], M = Math.max(A.length, B.length);
        for (let j = 0; j < M; j++) { const za = A[j] ? A[j].z : 0, zb = B[j] ? B[j].z : 0; if (za !== zb) return za - zb; }
      }
      const N1 = [], N2 = []; let size = 0, any = false;
      const ch = shallow ? childrenUnsorted : children;
      for (const set of F1) for (const nd of set) { const c = ch(nd); if (c.length) any = true; size += c.length; N1.push(c); }
      for (const set of F2) for (const nd of set) { const c = ch(nd); if (c.length) any = true; N2.push(c); }
      if (!any || size > 6000) return 0;
      F1 = N1; F2 = N2;
    }
    return 0;
  }
  // Substituenten eines Zentrums nach Priorität sortieren; null bei Gleichstand
  function rank(c) {
    const center = {atom: c, parent: null};
    const subs = adj[c].map(({n}) => ({atom: n, z: z(n), parent: center}));
    for (let k = 0; k < hc[c]; k++) subs.push({atom: -1, z: 1, dup: true});
    if (subs.length !== 4) return null;
    subs.sort((x, y) => -cmp(x, y));
    for (let k = 0; k < 3; k++) if (cmp(subs[k], subs[k + 1]) === 0) return null;
    return subs.map(x => x.atom);   // -1 = H
  }
  // zwei Substituenten (Atom bzw. -1 = H) an c vergleichen: > 0, wenn a höhere Priorität hat
  function cmpSub(c, a, b) {
    const center = {atom: c, parent: null};
    const node = x => x < 0 ? {atom: -1, z: 1, dup: true} : {atom: x, z: z(x), parent: center};
    return cmp(node(a), node(b));
  }
  return {rank, ordOf, cmpSub};
}

// asymmetrische C-Atome: sp3, vier verschiedene Substituenten
function stereoCenters(g, hc) {
  hc = hc || hCounts(g);
  const cip = makeCip(g, hc), out = [];
  g.atoms.forEach((a, i) => {
    if (a.el !== 'C') return;
    const bs = g.bondsOf(i);
    if (bs.some(b => b.order !== 1)) return;
    if (bs.length + hc[i] !== 4 || hc[i] > 1) return;
    const r = cip.rank(i);
    if (r) out.push({atom: i, ranked: r});
  });
  return out;
}

const parity = arr => { let p = 0; for (let i = 0; i < arr.length; i++) for (let j = i + 1; j < arr.length; j++) if (arr[i] > arr[j]) p ^= 1; return p; };
// SMILES-Chiralität → R/S
function cipFromSmiles(order, tag, ranked) {
  const rk = order.map(n => ranked.indexOf(n));
  if (rk.length !== 4 || rk.some(x => x < 0)) return null;
  let t = tag === '@' ? 0 : 1;
  t ^= parity(rk) ^ parity([3, 0, 1, 2]);
  return t === 0 ? 'R' : 'S';
}
// Fischer: Positionen oben/rechts/unten/links → R/S
function cipFromFischer(ranked, up, right, down, left) {
  const posOf = n => [up, right, down, left].indexOf(n);
  const i0 = posOf(ranked[0]), i1 = posOf(ranked[1]), i2 = posOf(ranked[2]), i3 = posOf(ranked[3]);
  const cw = ((i1 - i0 + 4) % 4) < ((i2 - i0 + 4) % 4);
  const lowVertical = i3 === 0 || i3 === 2;
  return (lowVertical ? cw : !cw) ? 'R' : 'S';
}
// 2D mit Keil (z>0) / Strich (z<0): Vektoren in Bildschirmkoordinaten (y nach unten)
function cipFromVectors(vecs) { // vecs: Vektoren der Prioritäten 1..3 als {x,y,z}
  const v = vecs.map(p => ({x: p.x, y: -p.y, z: p.z}));
  const c = {x: v[1].y * v[2].z - v[1].z * v[2].y, y: v[1].z * v[2].x - v[1].x * v[2].z, z: v[1].x * v[2].y - v[1].y * v[2].x};
  const t = v[0].x * c.x + v[0].y * c.y + v[0].z * c.z;
  if (Math.abs(t) < 1e-3) return null;
  return t < 0 ? 'R' : 'S';
}

/* ---------- Fischer-Projektion ---------- */
function branchAtoms(g, root, from) {
  const seen = new Set([from, root]), st = [root], out = [root];
  while (st.length) { const u = st.pop(); for (const v of g.nb(u)) if (!seen.has(v)) { seen.add(v); out.push(v); st.push(v); } }
  return out;
}
const SUBD = '₀₁₂₃₄₅₆₇₈₉';
function hillText(g, hc, atomsList) {
  const c = {};
  for (const i of atomsList) { const el = g.atoms[i].el; c[el] = (c[el] || 0) + 1; if (hc[i]) c.H = (c.H || 0) + hc[i]; }
  const keys = Object.keys(c), ord = c.C ? ['C', ...(c.H ? ['H'] : []), ...keys.filter(k => k !== 'C' && k !== 'H').sort()] : keys.sort();
  return ord.map(k => k + (c[k] > 1 ? c[k] : '')).join('');
}
// Gruppenformel (Halbstrukturformel) eines Astes, z. B. COOH, CH2OH, CH(CH3)2
function condensed(g, hc, root, from, ringSet) {
  if (root < 0) return 'H';
  if (!ringSet) ringSet = new Set(findRings(g).flat());
  const br = branchAtoms(g, root, from);
  if (ringSet.has(root)) {
    // Ringgerüst als Summenformel (z. B. C6H4), Substituenten daran als Gruppen (OH, CH3 …)
    const comp = new Set([root]), st = [root];
    while (st.length) { const u = st.pop(); for (const v of g.nb(u)) if (v !== from && ringSet.has(v) && !comp.has(v)) { comp.add(v); st.push(v); } }
    const att = [];
    for (const u of comp) for (const v of g.nb(u)) if (v !== from && !comp.has(v)) att.push(condensed(g, hc, v, u, ringSet));
    const cnt = {}; att.forEach(t => cnt[t] = (cnt[t] || 0) + 1);
    const keys = Object.keys(cnt);
    let txt = hillText(g, hc, [...comp]);
    for (const k of keys) txt += (cnt[k] > 1 || keys.length > 1) ? '(' + k + ')' + (cnt[k] > 1 ? cnt[k] : '') : k;
    return txt;
  }
  const a = g.atoms[root], h = hc[root];
  const chg = x => { const c = g.atoms[x].charge || 0; return c ? (Math.abs(c) > 1 ? String(Math.abs(c)).replace(/\d/g, d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]) : '') + (c > 0 ? '⁺' : '⁻') : ''; };
  let s = a.el + (h ? 'H' + (h > 1 ? h : '') : '') + chg(root);
  const kids = g.nb(root).filter(n => n !== from);
  const term = [], non = [];
  for (const k of kids) (g.nb(k).length === 1 ? term : non).push(k);
  const ordOf = k => { const b = g.bonds.find(x => (x.a === root && x.b === k) || (x.b === root && x.a === k)); return b ? b.order : 1; };
  term.sort((x, y) => (ordOf(y) - ordOf(x)) || (g.atoms[x].el === 'O' ? -1 : 1));
  const groups = [];
  for (const k of term) {
    const t = g.atoms[k].el + (hc[k] ? 'H' + (hc[k] > 1 ? hc[k] : '') : '') + chg(k);
    const last = groups[groups.length - 1];
    if (last && last.t === t && ordOf(k) === 1 && last.single) last.n++; else groups.push({t, n: 1, single: ordOf(k) === 1});
  }
  // Teile: endständige Gruppen, dann weiterführende Äste; alle außer dem letzten in Klammern (=O bleibt ohne)
  const parts = groups.map(gr => gr.n > 1
    ? {t: (gr.t.length > 1 || /H/.test(gr.t) ? '(' + gr.t + ')' + gr.n : gr.t + gr.n), wrap: false}
    : {t: gr.t, wrap: gr.single});
  non.forEach(k => parts.push({t: condensed(g, hc, k, root, ringSet), wrap: true}));
  parts.forEach((pt, i) => { s += (i < parts.length - 1 && pt.wrap) ? '(' + pt.t + ')' : pt.t; });
  return s;
}
const MIRROR = {OH: 'HO', NH2: 'H2N', SH: 'HS', CH3: 'H3C', COOH: 'HOOC', CHO: 'OHC', CH2OH: 'HOH2C', NO2: 'O2N', CN: 'NC', OCH3: 'H3CO', NHCH3: 'H3CHN', CH2CH3: 'H3CH2C', CH2NH2: 'H2NH2C', CH2SH: 'HSH2C', COOCH3: 'H3COOC', CONH2: 'H2NOC'};
// Gruppenformel spiegeln (links vom Molekül): Bausteine in umgekehrter Reihenfolge, H vor das Atom (CH2NH2 → H2NH2C)
const mirrorText = t => MIRROR[t] || ({'NH3⁺': '⁺H3N', 'O⁻': '⁻O', 'COO⁻': '⁻OOC'})[t] || mirrorUnits(t);
function mirrorUnits(t) {
  const re = /^([A-Z][a-z]?)(H\d*)?(\d*)([⁺⁻]\d*)?/, units = [];
  let r = t;
  while (r.length) {
    const m = r.match(re);
    if (!m || m[1] === 'H' || !m[0]) return t;          // freies H, Klammern u. Ä.: unverändert
    if (m[3] && m[2]) return t;                          // C6H5 o. Ä. (Summenformel-artig)
    units.push(m); r = r.slice(m[0].length);
  }
  if (units.length < 1) return t;
  return units.reverse().map(m => (m[4] || '') + (m[2] || '') + m[1] + (m[3] || '')).join('');
}

function fischerLayout(g, hc, cipMap) {
  hc = hc || hCounts(g);
  const centers = stereoCenters(g, hc);
  if (!centers.length) return {ok: false, reason: 'Das Molekül hat kein asymmetrisches C-Atom.'};
  const rings = findRings(g), ringAtom = new Set(rings.flat());
  const acyc = centers.filter(c => !ringAtom.has(c.atom));
  if (!acyc.length) return {ok: false, reason: 'Die Stereozentren liegen in Ringen – die Fischer-Projektion gibt es nur für offenkettige Stereozentren.'};
  const isCenter = new Map(centers.map(c => [c.atom, c]));
  // Kohlenstoffgerüst ohne Ringatome (Wald)
  const cAdj = i => g.nb(i).filter(n => g.atoms[n].el === 'C' && !ringAtom.has(n));
  const pathBetween = (a, b) => {
    const prev = new Map([[a, null]]), q = [a];
    while (q.length) { const u = q.shift(); if (u === b) break; for (const v of cAdj(u)) if (!prev.has(v)) { prev.set(v, u); q.push(v); } }
    if (!prev.has(b)) return null;
    const p = []; let c = b; while (c !== null) { p.unshift(c); c = prev.get(c); } return p;
  };
  let best = null;
  for (const a of acyc) for (const b of acyc) {
    const p = pathBetween(a.atom, b.atom); if (!p) continue;
    const n = p.filter(x => isCenter.has(x) && !ringAtom.has(x)).length;
    if (!best || n > best.n || (n === best.n && p.length > best.p.length)) best = {p, n};
  }
  let chain = best.p;
  const size = (root, from) => root < 0 ? 0 : branchAtoms(g, root, from).length;
  const subsOf = c => [...g.nb(c), ...Array(hc[c]).fill(-1)];
  // senkrechte Enden wählen
  function pickEnd(c, exclude) {
    const cand = subsOf(c).filter(n => n >= 0 && !exclude.includes(n));
    const carb = cand.filter(n => g.atoms[n].el === 'C');
    const pool = carb.length ? carb : cand;
    if (!pool.length) return null;
    const ranked = isCenter.get(c)?.ranked || [];
    pool.sort((x, y) => (carb.length ? size(y, c) - size(x, c) : ranked.indexOf(x) - ranked.indexOf(y)));
    return pool[0];
  }
  let top, bottom;
  if (chain.length === 1) {
    const c = chain[0];
    const cand = subsOf(c).filter(n => n >= 0);
    const carb = cand.filter(n => g.atoms[n].el === 'C').sort((x, y) => size(y, c) - size(x, c));
    const ranked = isCenter.get(c).ranked;
    const rest = cand.filter(n => !carb.includes(n)).sort((x, y) => ranked.indexOf(x) - ranked.indexOf(y));
    const two = [...carb, ...rest].slice(0, 2);
    if (two.length < 2) return {ok: false, reason: 'Für dieses Stereozentrum lässt sich keine senkrechte Kette bilden.'};
    [top, bottom] = two;
  } else {
    top = pickEnd(chain[0], [chain[1]]);
    bottom = pickEnd(chain[chain.length - 1], [chain[chain.length - 2]]);
  }
  // Ausrichtung: höher oxidiertes Ende nach oben, sonst kleinerer Ast oben (C1 oben)
  const ox = (r, from) => {
    if (r === null || r < 0) return -1;
    const a = g.atoms[r]; if (a.el !== 'C') return 10;
    return g.bondsOf(r).filter(b => (b.a === r ? b.b : b.a) !== from).reduce((s, b) => { const n = b.a === r ? b.b : b.a; return s + (['O', 'N', 'S', 'F', 'Cl', 'Br', 'I'].includes(g.atoms[n].el) ? b.order : 0); }, 0);
  };
  // eine Ebene tiefer mitzählen (CH2–CHO der 2-Desoxyribose ist das oxidierte Ende)
  const oxD = (r, from) => {
    const o = ox(r, from); if (r === null || r < 0 || g.atoms[r].el !== 'C') return o;
    let deep = 0; for (const n of g.nb(r)) if (n !== from && g.atoms[n].el === 'C' && !ringAtom.has(n)) deep = Math.max(deep, ox(n, r));
    return o + 0.6 * deep;
  };
  const oxT = oxD(top, chain[0]), oxB = oxD(bottom, chain[chain.length - 1]);
  const flip = oxB > oxT || (oxB === oxT && top !== null && bottom !== null && size(bottom, chain[chain.length - 1]) < size(top, chain[0]));
  if (flip) { chain = [...chain].reverse(); [top, bottom] = [bottom, top]; }
  // Keto-C (C=O mit genau einem weiteren C) gehört zur Kette: z. B. C2 der Fructose, darüber C1 als CH2OH
  const ketoO = c => { const b = g.bonds.find(x => (x.a === c || x.b === c) && x.order === 2 && g.atoms[x.a === c ? x.b : x.a].el === 'O' && g.nb(x.a === c ? x.b : x.a).length === 1); return b ? (b.a === c ? b.b : b.a) : null; };
  const isKeto = (c, from) => c !== null && c >= 0 && g.atoms[c].el === 'C' && !ringAtom.has(c) && ketoO(c) !== null &&
    g.nb(c).filter(n => n !== from && n !== ketoO(c)).length === 1 && g.nb(c).filter(n => n !== from && n !== ketoO(c) && g.atoms[n].el === 'C').length === 1;
  const oxTopOrig = Math.max(ox(top, chain[0]), oxD(top, chain[0]) > 1 ? 2 : 0);
  while (isKeto(top, chain[0])) { const c = top, next = g.nb(c).find(n => n !== chain[0] && n !== ketoO(c)); chain.unshift(c); top = next; }
  while (isKeto(bottom, chain[chain.length - 1])) { const c = bottom, next = g.nb(c).find(n => n !== chain[chain.length - 1] && n !== ketoO(c)); chain.push(c); bottom = next; }
  const newCip = {};
  const rows = chain.map((c, k) => {
    const up = k === 0 ? top : chain[k - 1], down = k === chain.length - 1 ? bottom : chain[k + 1];
    const others = subsOf(c).filter(n => n !== up && n !== down || (n === -1));
    // H (-1) kann zweimal vorkommen (CH2) – Liste korrekt bilden
    const hor = []; let hLeft = hc[c];
    for (const n of g.nb(c)) if (n !== up && n !== down) hor.push(n);
    while (hLeft-- > 0) hor.push(-1);
    let left = hor[0] ?? null, right = hor[1] ?? null;
    // nur ein waagerechter Partner (C=O in der Kette): nach rechts, als Doppelbindung
    let dbl = false;
    if (hor.length === 1) {
      right = hor[0]; left = null;
      const b = g.bonds.find(x => (x.a === c && x.b === right) || (x.b === c && x.a === right));
      dbl = !!b && b.order === 2;
    }
    const ctr = isCenter.get(c);
    let cip = null;
    if (ctr && !ringAtom.has(c) && hor.length === 2) {
      const want = cipMap ? cipMap[c] : null;
      if (want) {
        if (cipFromFischer(ctr.ranked, up, right, down, left) !== want) [left, right] = [right, left];
      } else {
        // Standard: Substituent mit höherer Priorität (nicht H) nach rechts
        if (ctr.ranked.indexOf(left) < ctr.ranked.indexOf(right)) [left, right] = [right, left];
        newCip[c] = cipFromFischer(ctr.ranked, up, right, down, left);
      }
      cip = cipFromFischer(ctr.ranked, up, right, down, left);
    }
    return {atom: c, isCenter: !!cip, cip, up, down,
      left: left === null ? null : {atom: left, text: mirrorText(condensed(g, hc, left, c, ringAtom))},
      right: right === null ? null : {atom: right, text: condensed(g, hc, right, c, ringAtom), double: dbl}};
  });
  const topT = top === null ? null : {atom: top, text: condensed(g, hc, top, chain[0], ringAtom)};
  const botT = bottom === null ? null : {atom: bottom, text: condensed(g, hc, bottom, chain[chain.length - 1], ringAtom)};
  // D/L: bei Aminosäuren das α-C-Atom (COOH oben, NH2 am ersten Zentrum), sonst das unterste Stereozentrum
  let dl = null;
  const first = rows[0], isN = x => x && x.atom >= 0 && g.atoms[x.atom].el === 'N';
  const aminoAcid = topT && topT.text === 'COOH' && first.isCenter && (isN(first.left) || isN(first.right));
  const lastC = aminoAcid ? first : [...rows].reverse().find(r => r.isCenter);
  const oxTop = oxTopOrig;
  if (lastC && lastC.left && lastC.right && oxTop >= 2) {
    const L = lastC.left.atom, Rr = lastC.right.atom, het = n => n >= 0 && ['O', 'N'].includes(g.atoms[n].el);
    if (L === -1 && het(Rr)) dl = 'D'; else if (Rr === -1 && het(L)) dl = 'L';
  }
  const skipped = centers.length - rows.filter(r => r.isCenter).length;
  return {ok: true, rows, top: topT, bottom: botT, dl, newCip, skipped};
}

/* ---------- Haworth-Projektion ---------- */
// R/S aus 3D-Vektoren (rechtshändig: x rechts, y oben, z zum Betrachter), Prioritäten 1..3
function cip3(v) {
  const c = {x: v[1].y * v[2].z - v[1].z * v[2].y, y: v[1].z * v[2].x - v[1].x * v[2].z, z: v[1].x * v[2].y - v[1].y * v[2].x};
  const t = v[0].x * c.x + v[0].y * c.y + v[0].z * c.z;
  if (Math.abs(t) < 1e-6) return null;
  return t < 0 ? 'R' : 'S';
}
const HW6 = [[0.4, -0.866], [1.15, 0], [0.55, 0.866], [-0.55, 0.866], [-1.15, 0], [-0.4, -0.866]];
const HW5 = [[0, -0.95], [1.2, -0.15], [0.5, 0.85], [-0.5, 0.85], [-1.2, -0.15]];
function haworthLayout(g, hc, cipMap, forceUp, turn) {
  hc = hc || hCounts(g);
  const satRing = r => r.every((a, k) => { const b = r[(k + 1) % r.length]; const bd = g.bonds.find(x => (x.a === a && x.b === b) || (x.a === b && x.b === a)); return bd && bd.order === 1; });
  const allRings = findRings(g);
  const rings = allRings.filter(r => (r.length === 5 || r.length === 6) && satRing(r));
  if (!rings.length) return {ok: false, reason: 'Kein gesättigter Fünf- oder Sechsring – die Haworth-Projektion ist für Ringzucker (Pyranosen, Furanosen) gedacht.'};
  const centers = stereoCenters(g, hc), isCenter = new Map(centers.map(c => [c.atom, c]));
  const ringSet = new Set(allRings.flat());
  const newCip = {};
  const sizeAll = (root, from) => root < 0 ? 0 : branchAtoms(g, root, from).length;
  // einen Ring anordnen; dir: Umlaufrichtung (0 = nach Faustregel, 1 = umgekehrt), rot: um die senkrechte Achse gedreht (Ring-O vorne)
  function layRing(ring, dir, rot, linkNb, tn) {
    const inRing = new Set(ring);
    const exoOf = a => [...g.nb(a).filter(n => !inRing.has(n)), ...Array(hc[a]).fill(-1)];
    const hetero = a => exoOf(a).filter(n => n >= 0 && g.atoms[n].el !== 'C').length;
    const size = (root, from) => linkNb && linkNb.has(root) ? 1 : sizeAll(root, from);
    // Reihenfolge: Ring-O hinten rechts, dann das anomere C (rechts)
    let start = ring.findIndex(a => g.atoms[a].el === 'O');
    let order;
    if (start >= 0) {
      const m = ring.length, n1 = ring[(start + 1) % m], n2 = ring[(start - 1 + m) % m];
      const s1 = hetero(n1) * 10 - exoOf(n1).reduce((s, n) => s + size(n, n1), 0);
      const s2 = hetero(n2) * 10 - exoOf(n2).reduce((s, n) => s + size(n, n2), 0);
      let fwd = s1 >= s2; if (dir) fwd = !fwd;
      order = [ring[start]];
      if (fwd) for (let k = 1; k < m; k++) order.push(ring[(start + k) % m]);
      else for (let k = 1; k < m; k++) order.push(ring[(start - k + m) % m]);
    } else order = dir ? [...ring].reverse() : [...ring];
    const T = ring.length === 6 ? HW6 : HW5, sg = rot ? -1 : 1;
    // tn: Ring um die senkrechte Achse weiterdrehen (Positionen der Vorlage; vordere Kante bleibt unten)
    const m0 = T.length, sh = (((tn || 0) % m0) + m0) % m0;
    const P = new Map(order.map((a, k) => [a, {x: sg * T[(k + sh) % m0][0], y: 0, z: sg * T[(k + sh) % m0][1]}]));
    const cipOut = {};
    const atomsOut = order.map(a => {
      const exo = exoOf(a);
      let up = null, down = null, cip = null;
      if (exo.length >= 2) {
        const ctr = isCenter.get(a);
        const heavy = exo.filter(n => n >= 0);
        let pair = [exo[0], exo[1]];
        if (ctr) pair.sort((x, y) => ctr.ranked.indexOf(x) - ctr.ranked.indexOf(y));
        else pair.sort((x, y) => size(y, a) - size(x, a));
        [up, down] = pair;
        if (forceUp) {
          const upHeavy = heavy.find(n => forceUp.has(n));
          if (upHeavy !== undefined) { up = upHeavy; down = exo.find(n => n !== upHeavy); }
          else if (heavy.length) { down = heavy[0]; up = exo.find(n => n !== down); }
        }
        if (ctr) {
          const vec = n => {
            if (n === up) return {x: 0, y: 1, z: 0};
            if (n === down) return {x: 0, y: -1, z: 0};
            const p = P.get(n), q = P.get(a); return {x: p.x - q.x, y: 0, z: p.z - q.z};
          };
          const desc = () => cip3(ctr.ranked.slice(0, 3).map(vec));
          const want = forceUp ? null : (cipMap ? cipMap[a] : null);
          if (want && desc() !== want) [up, down] = [down, up];
          cip = desc();
          if (!want) cipOut[a] = cip;
        }
      } else if (exo.length === 1) up = exo[0];
      const lab = n => n === null ? null : (linkNb && linkNb.has(n) ? {atom: n, text: '', link: true} : {atom: n, text: n < 0 ? 'H' : condensed(g, hc, n, a, ringSet)});
      return {atom: a, x: P.get(a).x, z: P.get(a).z, el: g.atoms[a].el, isCenter: !!cip, cip, up: lab(up), down: lab(down)};
    });
    return {atoms: atomsOut, cip: cipOut, size: ring.length, ring};
  }
  const score = r => (r.filter(a => g.atoms[a].el === 'O').length === 1 ? 10 : 0) + r.filter(a => isCenter.has(a)).length;
  // mehrere Zuckerringe, verbunden über Brücken (glykosidisches O, ggf. CH₂–O): nebeneinander zeichnen
  const multi = (() => {
    const sug = rings.filter(r => r.filter(a => g.atoms[a].el === 'O').length === 1);
    if (sug.length < 2 || forceUp) return null;
    const ringOf = new Map(); sug.forEach((r, i) => r.forEach(a => ringOf.set(a, i)));
    const links = [], seenPair = new Set();
    sug.forEach((r, i) => r.forEach(a => {
      for (const n1 of g.nb(a)) {
        if (ringSet.has(n1)) continue;
        const tryEnd = (path, last) => { for (const b of g.nb(last)) { if (b === a || path.includes(b)) continue; const j = ringOf.get(b); if (j !== undefined && j !== i) { const key = Math.min(i, j) + '-' + Math.max(i, j); if (!seenPair.has(key)) { seenPair.add(key); links.push({ri: i, a, rj: j, b, path}); } return true; } } return false; };
        if (tryEnd([n1], n1)) continue;
        for (const n2 of g.nb(n1)) { if (n2 === a || ringSet.has(n2)) continue; if (tryEnd([n1, n2], n2)) break; }
      }
    }));
    if (links.length !== sug.length - 1) return null;   // nicht zusammenhängend oder verzweigt über Mehrfachbrücken
    const deg = sug.map((_, i) => links.filter(l => l.ri === i || l.rj === i).length);
    if (deg.some(d => d > 2)) return null;              // verzweigt: nicht darstellbar
    // reduzierendes Ende (freies OH am anomeren C) kommt nach rechts
    const freeAnom = i => sug[i].some(a => { const os = g.nb(a).filter(n => g.atoms[n].el === 'O'); return os.length >= 2 && os.some(o => !ringSet.has(o) && g.nb(o).length === 1); });
    const ends = sug.map((_, i) => i).filter(i => deg[i] === 1);
    let first = ends.find(i => !freeAnom(i)); if (first === undefined) first = ends[0];
    if (ends.length === 2 && freeAnom(ends[0]) === freeAnom(ends[1])) first = sug[ends[0]].length >= sug[ends[1]].length ? ends[0] : ends[1];
    const seq = [first], used = new Set([first]), lk = [];
    while (seq.length < sug.length) {
      const cur = seq[seq.length - 1];
      const l = links.find(x => (x.ri === cur && !used.has(x.rj)) || (x.rj === cur && !used.has(x.ri)));
      if (!l) return null;
      const fwd = l.ri === cur;
      lk.push(fwd ? {a: l.a, b: l.b, path: l.path} : {a: l.b, b: l.a, path: [...l.path].reverse()});
      const nx = fwd ? l.rj : l.ri; seq.push(nx); used.add(nx);
    }
    // Lage je Ring (alle Drehstellungen und beide Umlaufrichtungen): Verbindung zum linken Nachbarn links, zum rechten rechts;
    // die Brücke C–O–C soll eine gerade Linie sein: zeigt die Bindung am linken Ring nach unten, muss sie am rechten nach oben zeigen (und umgekehrt)
    const out = [];
    let prevSide = null;
    seq.forEach((ri, k) => {
      const ring = sug[ri], inA = k > 0 ? lk[k - 1].b : null, outA = k < seq.length - 1 ? lk[k].a : null;
      const linkNb = new Set([...(inA !== null ? [lk[k - 1].path[lk[k - 1].path.length - 1]] : []), ...(outA !== null ? [lk[k].path[0]] : [])]);
      const sideOf = (L, id) => { const q = L.atoms.find(z => z.atom === id); return q.up && q.up.link ? 'up' : q.down && q.down.link ? 'down' : null; };
      let best = null;
      for (let d = 0; d < 2; d++) for (let tn = 0; tn < ring.length; tn++) {
        const L = layRing(ring, d, 0, linkNb, tn), X = id => L.atoms.find(q => q.atom === id).x;
        let sc = 10 * ((outA !== null ? X(outA) : 0) - (inA !== null ? X(inA) : 0));
        if (inA !== null && prevSide && sideOf(L, inA) === prevSide) sc -= 100;
        if (L.atoms.some(q => g.atoms[q.atom].el === 'O' && q.z < 0)) sc += 1;
        if (!d && !tn) sc += 0.2;
        if (!best || sc > best.sc) best = {sc, L};
      }
      const L = best.L;
      prevSide = outA !== null ? sideOf(L, outA) : null;
      out.push({...L, ox: 0});
    });
    const linkOut = lk.map(l => ({a: l.a, b: l.b, path: l.path.map(n => ({atom: n, el: g.atoms[n].el, text: g.atoms[n].el + (hc[n] ? 'H' + (hc[n] > 1 ? hc[n] : '') : '')}))}));
    return {rings: out, links: linkOut};
  })();
  if (multi) {
    multi.rings.forEach(r => Object.assign(newCip, r.cip));
    const inAny = new Set(multi.rings.flatMap(r => r.ring));
    const red = multi.rings[multi.rings.length - 1];
    const skipped = centers.filter(c => !inAny.has(c.atom) && !multi.links.some(l => l.path.some(q => q.atom === c.atom))).length;
    return {ok: true, atoms: red.atoms, rings: multi.rings.map(r => ({atoms: r.atoms, size: r.size})), links: multi.links, newCip, skipped, size: red.size};
  }
  const ring = [...rings].sort((x, y) => score(y) - score(x))[0];
  const L = layRing(ring, 0, 0, null, turn);
  Object.assign(newCip, L.cip);
  const inRing = new Set(ring);
  const skipped = centers.filter(c => !inRing.has(c.atom)).length;
  return {ok: true, atoms: L.atoms, newCip, skipped, size: ring.length};
}
// Ringzucker: β-D-Form über „oben“-Atome festgelegt
const RINGSUGAR = {
  'glucopyranose': ['OCC1OC(O)C(O)C(O)C1O', [5, 9, 1], 4], 'galactopyranose': ['OCC1OC(O)C(O)C(O)C1O', [5, 9, 11, 1], 4],
  'mannopyranose': ['OCC1OC(O)C(O)C(O)C1O', [5, 7, 9, 1], 4], 'ribofuranose': ['OCC1OC(O)C(O)C1O', [5, 1], 4],
  'desoxyribofuranose': ['OCC1OC(O)CC1O', [5, 1], 4], 'deoxyribofuranose': ['OCC1OC(O)CC1O', [5, 1], 4],
  '2-desoxyribofuranose': ['OCC1OC(O)CC1O', [5, 1], 4], '2-deoxyribofuranose': ['OCC1OC(O)CC1O', [5, 1], 4],
  'fructofuranose': ['OCC1OC(O)(CO)C(O)C1O', [5, 9, 1], 4],
};
function ringSugarGraph(key, anomer, dl) {
  const [smi, up, anomC] = RINGSUGAR[key];
  const g = parseSmiles(smi);
  const f = haworthLayout(g, null, null, new Set(up));
  for (const a of f.atoms) if (a.cip) g.atoms[a.atom].cip = a.cip;
  const flip = i => { const a = g.atoms[i]; if (a.cip) a.cip = a.cip === 'R' ? 'S' : 'R'; };
  if (dl === 'L') g.atoms.forEach((_, i) => flip(i));
  g.atoms[anomC].betaCip = g.atoms[anomC].cip;
  if (anomer === 'A') flip(anomC);
  return g;
}

// Deskriptoren am Graphen bestimmen: aus SMILES (@/@@) und aus dem Namen
function assignStereo(g, st) {
  const hc = hCounts(g);
  const centers = stereoCenters(g, hc);
  const warn = [];
  for (const c of centers) {
    const a = g.atoms[c.atom];
    if (a.chiralTag && a.nbOrder) { const d = cipFromSmiles(a.nbOrder, a.chiralTag, c.ranked); if (d) a.cip = d; }
  }
  if (st && st.list.length) {
    const byLoc = st.list.filter(x => x.loc), noLoc = st.list.filter(x => !x.loc);
    for (const x of byLoc) {
      const idx = g.mainChain ? g.mainChain[x.loc - 1] : undefined;
      if (idx === undefined || !centers.some(c => c.atom === idx)) { warn.push('C' + x.loc + ' ist kein Stereozentrum'); continue; }
      g.atoms[idx].cip = x.d;
    }
    if (noLoc.length) {
      if (centers.length === noLoc.length) centers.forEach((c, k) => { g.atoms[c.atom].cip = noLoc[k].d; });
      else warn.push('Deskriptor ohne Position passt nicht eindeutig');
    }
  }
  if (st && st.dl) {
    const cipMap = {}; g.atoms.forEach((a, i) => { if (a.cip) cipMap[i] = a.cip; });
    const f = fischerLayout(g, hc, cipMap);
    if (f.ok) {
      Object.entries(f.newCip).forEach(([i, d]) => { g.atoms[i].cip = d; });
      if (f.dl && f.dl !== st.dl) centers.forEach(c => { const a = g.atoms[c.atom]; if (a.cip) a.cip = a.cip === 'R' ? 'S' : 'R'; });
      if (!f.dl) warn.push('D/L lässt sich hier nicht zuordnen');
    } else warn.push('D/L nicht anwendbar');
  }
  return {centers: centers.length, warn};
}
function extractStereo(raw) {
  let s = raw.trim(); const st = {list: [], dl: null};
  let m = s.match(/^\(\s*(\d*\s*[RSrs](?:\s*,\s*\d*\s*[RSrs])*)\s*\)\s*-\s*/);
  if (m) {
    st.list = m[1].split(',').map(x => { const mm = x.trim().match(/^(\d*)\s*([RSrs])$/); return {loc: mm[1] ? +mm[1] : null, d: mm[2].toUpperCase()}; });
    s = s.slice(m[0].length);
  }
  s = s.replace(/^\(\s*[+\-±]\s*\)\s*-\s*/, '');
  m = s.match(/^(α|β|alpha|beta|a|b)\s*-\s*(?=[DLdl]\s*-)/i);
  if (m) { st.anomer = /^(α|alpha|a)$/i.test(m[1]) ? 'A' : 'B'; s = s.slice(m[0].length); }
  // „2-Desoxy-D-ribose“ → D + „2-Desoxyribose“
  m = s.match(/^(\d+\s*-\s*(?:desoxy|deoxy))\s*-\s*([DLdl])\s*-\s*/i);
  if (m) { st.dl = m[2].toUpperCase(); s = m[1].replace(/\s+/g, '') + s.slice(m[0].length); return {name: s, st}; }
  m = s.match(/^([DLdl])\s*-\s*(?=[A-Za-zÄÖÜäöü(]|\d+\s*-\s*(?:desoxy|deoxy))/i);
  if (m) { st.dl = m[1].toUpperCase(); s = s.slice(m[0].length); }
  return {name: s, st};
}

/* ======================= Eingabe ======================= */
// auch rein aromatische SMILES wie c1ccccc1 (Kleinbuchstaben mit Ringziffern)
function looksLikeSmiles(s) { return /^[A-Za-z0-9@+\-\[\]()=#$%\/\\.:*]+$/.test(s) && (/[A-Z]/.test(s) || /^[cnosp]+\d/.test(s) || /[cnosp]\d.*[cnosp]\d/.test(s)); }
// Gibt {g, pos, how} zurück oder wirft einen Fehler mit deutscher Meldung
function fromInput(raw) {
  const s0 = raw.trim();
  if (!s0) throw new Error('Bitte einen Namen oder SMILES eingeben.');
  let {name: s, st} = extractStereo(s0);
  let res = null, err = null;
  // Ringzucker: „β-D-Glucopyranose“, „α-D-Glucose“ (mit α/β → Ringform)
  let key = s.trim().toLowerCase().replace(/\s+/g, '');
  if (st.anomer && !RINGSUGAR[key]) { const alt = {glucose: 'glucopyranose', galactose: 'galactopyranose', mannose: 'mannopyranose', ribose: 'ribofuranose', fructose: 'fructofuranose', desoxyribose: 'desoxyribofuranose', deoxyribose: 'deoxyribofuranose'}[key]; if (alt) key = alt; }
  if (RINGSUGAR[key]) {
    res = {g: ringSugarGraph(key, st.anomer || 'B', st.dl || 'D'), how: 'trivial'};
    kekulize(res.g); res.stereo = {centers: stereoCenters(res.g).length, warn: []}; res.pos = layout(res.g);
    return res;
  }
  try { res = nameToGraph(s); } catch (e) { err = e; }
  if (!res && looksLikeSmiles(s0)) {
    try { res = {g: parseSmiles(s0), how: 'smiles'}; st = null; } catch (e) { err = err || e; }
  }
  if (!res) {
    const own = err && /[äöüß]|Locant|Position|Bindungen|Ring/.test(err.message);
    throw new Error(own ? err.message : 'Name nicht erkannt.');
  }
  kekulize(res.g);
  checkValence(res.g);
  res.stereo = assignStereo(res.g, st);
  res.pos = layout(res.g);
  return res;
}
function fromSmiles(s) { const g = parseSmiles(s); kekulize(g); const stereo = assignStereo(g, null); return {g, pos: layout(g), how: 'smiles', stereo}; }
const MAXV = {H: 1, B: 4, C: 4, Si: 4, N: 3, O: 2, F: 1, Cl: 7, Br: 7, I: 7, S: 6, P: 5};
function checkValence(g) {
  g.atoms.forEach((a, i) => {
    const sum = g.bondsOf(i).reduce((s, b) => s + b.order, 0);
    const max = (MAXV[a.el] ?? 8) + (a.el === 'N' || a.el === 'O' ? a.charge : 0);
    if (a.el === 'C' && sum > 4) throw new Error('Ein Kohlenstoffatom hätte ' + sum + ' Bindungen – bitte Locants prüfen.');
    if (sum > Math.max(max, 1) && a.el !== 'C') throw new Error('Ein ' + a.el + '-Atom hätte zu viele Bindungen – bitte Locants prüfen.');
  });
}

/* ======================= Struktur → IUPAC-Name (deutsch) ======================= */
// Substitutive Nomenklatur für die im Unterricht üblichen Stoffklassen. Liefert null, wenn das Molekül nicht sicher benannt werden kann.
const NM_STEM = ['', 'Meth', 'Eth', 'Prop', 'But', 'Pent', 'Hex', 'Hept', 'Oct', 'Non', 'Dec', 'Undec', 'Dodec', 'Tridec', 'Tetradec', 'Pentadec', 'Hexadec', 'Heptadec', 'Octadec', 'Nonadec', 'Icos', 'Henicos', 'Docos', 'Tricos', 'Tetracos'];
const NM_MULT = ['', '', 'di', 'tri', 'tetra', 'penta', 'hexa', 'hepta', 'octa', 'nona', 'deca', 'undeca', 'dodeca', 'trideca', 'tetradeca', 'pentadeca', 'hexadeca', 'heptadeca', 'octadeca', 'nonadeca', 'icosa'];
const NM_MULTC = ['', '', 'bis', 'tris', 'tetrakis', 'pentakis', 'hexakis', 'heptakis', 'octakis'];
const NM_HAL = {F: 'fluor', Cl: 'chlor', Br: 'brom', I: 'iod'};
// Rang der Hauptgruppen (kleiner = höher)
const NM_RANK = {acid: 1, amide: 3, nitrile: 4, al: 5, one: 6, ol: 7, thiol: 8, amine: 9};
class NameFail extends Error {}
const nmFail = () => { throw new NameFail('nicht benennbar'); };

// schwere Atome mit H-Zahl, Bindungsordnungen, Ringen
function nmGraph(g) {
  const hc = hCounts(g), map = new Map(), orig = [];
  g.atoms.forEach((a, i) => { if (a.el !== 'H') { map.set(i, orig.length); orig.push(i); } });
  const M = {n: orig.length, el: [], ch: [], h: [], adj: [], orig};
  orig.forEach((i, k) => { M.el[k] = g.atoms[i].el; M.ch[k] = g.atoms[i].charge || 0; M.h[k] = hc[i]; M.adj[k] = []; });
  for (const b of g.bonds) {
    const A = map.get(b.a), B = map.get(b.b);
    if (A === undefined && B === undefined) continue;
    if (A === undefined || B === undefined) { M.h[A === undefined ? B : A]++; continue; }
    const o = b.order === 'ar' ? 1.5 : b.order;
    M.adj[A].push({n: B, o}); M.adj[B].push({n: A, o});
  }
  M.ord = (a, b) => (M.adj[a].find(x => x.n === b) || {o: 0}).o;
  M.deg = a => M.adj[a].length;
  // Ringe (kleinste Ringe über die Bindungen)
  const rings = [];
  const seenKey = new Set();
  for (let a = 0; a < M.n; a++) for (const {n: b} of M.adj[a]) {
    if (b < a) continue;
    // kürzester Weg a→b ohne die Bindung a–b
    const prev = new Map([[a, -1]]), q = [a];
    while (q.length) { const u = q.shift(); if (u === b) break; for (const {n: v} of M.adj[u]) { if ((u === a && v === b) || prev.has(v)) continue; prev.set(v, u); q.push(v); } }
    if (!prev.has(b)) continue;
    const r = []; for (let x = b; x !== -1; x = prev.get(x)) r.push(x);
    const key = [...r].sort((p, q2) => p - q2).join(',');
    if (!seenKey.has(key)) { seenKey.add(key); rings.push(r); }
  }
  M.rings = rings;
  // Benzolringe (abwechselnd Einfach-/Doppelbindung) einheitlich als aromatisch
  for (const r of rings) {
    if (r.length !== 6 || r.some(x => M.el[x] !== 'C')) continue;
    const o = r.map((x, i) => M.ord(x, r[(i + 1) % 6]));
    if (o.every((v, i) => v === 1.5 || v + o[(i + 1) % 6] === 3)) r.forEach((x, i) => { const y = r[(i + 1) % 6]; M.adj[x].find(z => z.n === y).o = 1.5; M.adj[y].find(z => z.n === x).o = 1.5; });
  }
  M.ringOf = new Map(); rings.forEach((r, k) => r.forEach(x => { if (!M.ringOf.has(x)) M.ringOf.set(x, []); M.ringOf.get(x).push(k); }));
  M.inRing = x => M.ringOf.has(x);
  return M;
}
function nmRingInfo(M, k) {
  const r = M.rings[k], n = r.length;
  const arom = n === 6 && r.every((x, i) => { const y = r[(i + 1) % n]; return M.ord(x, y) === 1.5; }) ||
    n === 6 && (() => { const o = r.map((x, i) => M.ord(x, r[(i + 1) % n])); return o.every((v, i) => v + o[(i + 1) % n] === 3); })();
  return {r, n, arom};
}

// funktionelle Gruppe an einem C-Atom erkennen
function nmCarbonGroup(M, c) {
  if (M.el[c] !== 'C') return null;
  const nb = M.adj[c];
  const oxo = nb.find(x => M.el[x.n] === 'O' && x.o === 2 && M.deg(x.n) === 1);
  const cN = nb.filter(x => M.el[x.n] === 'C').length;
  if (oxo) {
    const oh = nb.find(x => M.el[x.n] === 'O' && x.o === 1 && M.deg(x.n) === 1 && M.h[x.n] === 1 && M.ch[x.n] === 0);
    if (oh && !M.inRing(c)) return {type: 'acid', atoms: [oxo.n, oh.n]};
    const or = nb.find(x => M.el[x.n] === 'O' && x.o === 1 && M.deg(x.n) === 2);
    if (nb.filter(x => M.el[x.n] !== 'C').length > 2) return {type: 'odd'};   // Carbonat, Carbamat …
    if (or && !M.inRing(c) && !M.inRing(or.n)) return {type: 'ester', atoms: [oxo.n, or.n], o: or.n};
    const nh = nb.find(x => M.el[x.n] === 'N' && x.o === 1);
    if (nh && !M.inRing(c)) { if (M.deg(nh.n) === 1 && M.h[nh.n] === 2) return {type: 'amide', atoms: [oxo.n, nh.n]}; return {type: 'subamide'}; }
    if (M.h[c] >= 1 && cN <= 1 && !M.inRing(c)) return {type: 'al', atoms: [oxo.n]};
    if (cN === 2) return {type: 'one', atoms: [oxo.n]};
    return {type: 'odd'};
  }
  const cn = nb.find(x => M.el[x.n] === 'N' && x.o === 3 && M.deg(x.n) === 1);
  if (cn && !M.inRing(c)) return {type: 'nitrile', atoms: [cn.n]};
  return null;
}

function nmName(g) {
  const M = nmGraph(g);
  if (!M.n) return null;
  // nur übliche Elemente; Ladungen nur in der Nitrogruppe
  for (let i = 0; i < M.n; i++) {
    if (!['C', 'H', 'O', 'N', 'S', 'F', 'Cl', 'Br', 'I'].includes(M.el[i])) return null;
    if (M.ch[i] && !(M.el[i] === 'N' && M.ch[i] === 1 && M.adj[i].filter(x => M.el[x.n] === 'O').length === 2) && !(M.el[i] === 'O' && M.ch[i] === -1 && M.adj[i].length === 1 && M.el[M.adj[i][0].n] === 'N' && M.ch[M.adj[i][0].n] === 1)) return null;
  }
  if (!M.el.includes('C')) return null;
  // Schwefel nur zweibindig (Thiol, Sulfid)
  for (let i = 0; i < M.n; i++) if (M.el[i] === 'S' && (M.deg(i) > 2 || M.adj[i].some(x => x.o !== 1))) return null;
  // nur einzelne Carbocyclen (keine Heterocyclen, keine kondensierten Ringe)
  for (const r of M.rings) if (r.some(x => M.el[x] !== 'C') || r.length < 3 || r.length > 12) return null;
  for (const [, ks] of M.ringOf) if (ks.length > 1) return null;
  // zusammenhängend?
  { const seen = new Set([0]), st = [0]; while (st.length) { const u = st.pop(); for (const {n} of M.adj[u]) if (!seen.has(n)) { seen.add(n); st.push(n); } } if (seen.size !== M.n) return null; }
  try {
    const groups = new Map();
    for (let c = 0; c < M.n; c++) { const gr = nmCarbonGroup(M, c); if (gr) { if (['subamide', 'odd'].includes(gr.type)) return null; groups.set(c, gr); } }
    const esters = [...groups].filter(([, v]) => v.type === 'ester');
    const hasAcid = [...groups.values()].some(v => v.type === 'acid');
    if (esters.length && !hasAcid) {
      if (esters.length > 1) return null;
      return {name: nmEster(M, groups, esters[0][0], esters[0][1]), locs: null};
    }
    const r = nmCore(M, groups, null);
    return {name: r.name, locs: new Map(r.seq.map((x, i) => [M.orig[x], i + 1]))};
  } catch (e) { if (e instanceof NameFail) return null; throw e; }
}

// Ester R–CO–O–R': R'-yl + Säurename auf -oat
function nmEster(M, groups, c, gr) {
  const o = gr.o;
  const r2 = M.adj[o].find(x => x.n !== c).n;
  // Säureteil: Ester-O als OH behandeln
  const M2 = {...M, h: [...M.h], adj: M.adj.map(l => [...l])};
  M2.adj[o] = M2.adj[o].filter(x => x.n !== r2); M2.adj[r2] = M2.adj[r2].filter(x => x.n !== o); M2.h[o] = 1;
  const keep = new Set(); { const st = [c]; keep.add(c); while (st.length) { const u = st.pop(); for (const {n} of M2.adj[u]) if (!keep.has(n)) { keep.add(n); st.push(n); } } }
  M2.mask = keep;
  const g2 = new Map(); for (const [k, v] of groups) if (keep.has(k)) g2.set(k, k === c ? {type: 'acid', atoms: [gr.atoms[0], o]} : v);
  const acid = nmCore(M2, g2, null).name;
  let anion;
  if (/benzoesäure$/i.test(acid)) anion = acid.replace(/oesäure$/i, 'oat');
  else if (/carbonsäure$/.test(acid)) anion = acid.replace(/carbonsäure$/, 'carboxylat');
  else if (/säure$/.test(acid)) anion = acid.replace(/säure$/, 'oat');
  else nmFail();
  const alk = nmSub(M, r2, o);
  // gebräuchliche deutsche Form „…säure…ylester“; bei Resten mit Lokanten „(…yl)…oat“
  if (!alk.complex) return nmCap(acid + alk.name + 'ester');
  return nmCap(nmWrap(alk.name) + nmLower(anion));
}
const nmLower = s => s.replace(/[A-Za-zÄÖÜäöü]/, ch => ch.toLowerCase());
// ersten Buchstaben groß – Lokanten wie „N-“, „N,N-“ am Anfang überspringen
const nmCap = s => s.replace(/^((?:[\d,]+-|N(?:,N)*-|[(\[{])*)([a-zäöü])/, (m, a, ch) => a + ch.toUpperCase());

// Hauptkette bzw. Ring wählen und benennen. fixed: Substituentenmodus {attach, from} (Kette/Ring beginnt am Anknüpfungsatom)
function nmCore(M, groups, fixed) {
  const inMask = x => !M.mask || M.mask.has(x);
  const blocked = new Set(fixed ? [fixed.from] : []);
  // Hauptgruppe bestimmen (nur im betrachteten Teil)
  const region = new Set();
  { const s0 = fixed ? fixed.attach : [...Array(M.n).keys()].find(inMask); const st = [s0]; region.add(s0); while (st.length) { const u = st.pop(); for (const {n} of M.adj[u]) if (!region.has(n) && !blocked.has(n) && inMask(n)) { region.add(n); st.push(n); } } }
  const hetGroup = x => {   // Heteroatom-Gruppe an einem C (ol, amin, thiol) – liefert Typ
    if (M.el[x] === 'O' && M.deg(x) === 1 && M.h[x] === 1 && M.adj[x][0].o === 1 && !groups.has(M.adj[x][0].n)) return 'ol';
    if (M.el[x] === 'S' && M.deg(x) === 1 && M.h[x] === 1) return 'thiol';
    if (M.el[x] === 'N' && M.ch[x] === 0 && M.adj[x].every(y => y.o === 1 && M.el[y.n] === 'C' && !groups.has(y.n)) && M.adj[x].length + M.h[x] === 3) return 'amine';
    return null;
  };
  let principal = null;
  if (!fixed) {
    for (const x of region) {
      const t = groups.get(x)?.type || hetGroup(x);
      if (t && NM_RANK[t] && (!principal || NM_RANK[t] < NM_RANK[principal])) principal = t;
    }
  }
  const isC = x => M.el[x] === 'C' && region.has(x);
  // Kandidaten: Ketten (azyklische C) und Ringe
  const cands = [];
  const acyc = x => isC(x) && !M.inRing(x);
  if (fixed) {
    if (M.inRing(fixed.attach)) cands.push({ring: M.ringOf.get(fixed.attach)[0]});
    else {
      // Ketten ab dem Anknüpfungsatom (Baum): alle Wege bis zu Blättern
      const walk = (u, prev, path) => { path.push(u); const nx = M.adj[u].filter(y => y.n !== prev && !blocked.has(y.n) && acyc(y.n) && !path.includes(y.n)); if (!nx.length) cands.push({chain: [...path]}); else { cands.push({chain: [...path]}); nx.forEach(y => walk(y.n, u, path)); } path.pop(); };
      walk(fixed.attach, fixed.from, []);
    }
  } else {
    const ac = [...region].filter(acyc);
    const leaves = ac.filter(x => M.adj[x].filter(y => acyc(y.n)).length <= 1);
    for (const a of leaves) {
      // Wege im Baum von a zu allen anderen Blättern (und a allein)
      const prev = new Map([[a, -1]]), q = [a];
      while (q.length) { const u = q.shift(); for (const {n} of M.adj[u]) if (acyc(n) && !prev.has(n)) { prev.set(n, u); q.push(n); } }
      for (const b of leaves) { if (b < a || !prev.has(b) || (b === a && M.adj[a].some(y => acyc(y.n)))) continue; const p = []; for (let x = b; x !== -1; x = prev.get(x)) p.push(x); cands.push({chain: p}); }
    }
    M.rings.forEach((r, k) => { if (r.every(x => region.has(x))) cands.push({ring: k}); });
    if (!cands.length) nmFail();
  }
  let best = null;
  for (const cd of cands) {
    for (const num of nmNumberings(M, cd, fixed)) {
      let ev = null;
      try { ev = nmEvaluate(M, groups, cd, num, principal, fixed, hetGroup, blocked); } catch (e) { if (!(e instanceof NameFail)) throw e; }
      if (!ev) continue;
      if (!best || nmCmp(ev.key, best.key) < 0) best = ev;
    }
  }
  if (!best) nmFail();
  return best;
}
function nmNumberings(M, cd, fixed) {
  if (cd.chain) {
    const p = cd.chain;
    if (fixed || p.length === 1) return [p];
    return [p, [...p].reverse()];
  }
  const r = M.rings[cd.ring], n = r.length, out = [];
  for (let s = 0; s < n; s++) for (const d of [1, -1]) {
    const seq = []; for (let k = 0; k < n; k++) seq.push(r[(s + d * k + n * 2) % n]);
    if (fixed && seq[0] !== fixed.attach) continue;
    out.push(seq);
  }
  return out;
}
// lexikographischer Vergleich von Schlüssel-Tupeln (Zahlen bzw. Zahlenlisten)
function nmCmp(a, b) {
  for (let i = 0; i < a.length; i++) {
    const x = a[i], y = b[i];
    if (Array.isArray(x)) { for (let k = 0; k < Math.max(x.length, y.length); k++) { const p = x[k] ?? Infinity, q = y[k] ?? Infinity; if (p !== q) return p < q ? -1 : 1; } }
    else if (typeof x === 'string') { if (x !== y) return x < y ? -1 : 1; }
    else if (x !== y) return x < y ? -1 : 1;
  }
  return 0;
}

// Kandidat mit Nummerierung bewerten und Namen erzeugen
function nmEvaluate(M, groups, cd, seq, principal, fixed, hetGroup, blocked) {
  const isRing = !!cd.ring || cd.ring === 0;
  const ringInfo = isRing ? nmRingInfo(M, cd.ring) : null;
  const pos = new Map(seq.map((x, i) => [x, i + 1]));
  const L = seq.length;
  const suff = [];      // Lokanten der Hauptgruppe
  const prefixes = [];  // {name, loc, complex}
  const handled = new Set(seq);
  // Mehrfachbindungen in der Hauptkette / im Ring
  const en = [], yn = [];
  if (!(ringInfo && ringInfo.arom)) {
    const pairs = isRing ? seq.map((x, i) => [x, seq[(i + 1) % L], i + 1, (i + 1) % L + 1]) : seq.slice(0, -1).map((x, i) => [x, seq[i + 1], i + 1, i + 2]);
    for (const [a, b, la, lb] of pairs) {
      const o = M.ord(a, b);
      const loc = isRing && Math.abs(la - lb) === L - 1 ? L : Math.min(la, lb);
      if (o === 2) en.push(loc); else if (o === 3) yn.push(loc); else if (o === 1.5) return null;
    }
  }
  // Substituenten und Gruppen an den Kettenatomen
  for (const x of seq) {
    const lx = pos.get(x);
    const gr = groups.get(x);
    if (gr && !isRing) {
      gr.atoms.forEach(a => handled.add(a));
      if (gr.type === principal) suff.push(lx);
      else if (gr.type === 'one') prefixes.push({name: 'oxo', loc: lx});
      else if (gr.type === 'al') { if (fixed && x === fixed.attach) prefixes.push({name: 'oxo', loc: lx}); else prefixes.push({name: 'oxo', loc: lx}); }
      else if (gr.type === 'acid') { if (fixed) return null; prefixes.push({name: 'carboxy', loc: lx}); }   // Säure innerhalb der Kette, aber nicht Hauptgruppe → nicht vorgesehen
      else if (gr.type === 'nitrile') { if (fixed && x === fixed.attach && L === 1) { prefixes.push({name: 'cyano', loc: 0}); } else return null; }
      else if (gr.type === 'amide') { handled.delete(gr.atoms[1]); prefixes.push({name: 'oxo', loc: lx}); }
      else if (gr.type === 'ester') { handled.delete(gr.atoms[1]); prefixes.push({name: 'oxo', loc: lx}); }   // –CO–OR in der Kette: oxo + Alkoxy
    } else if (gr && isRing && gr.type === 'one') {
      gr.atoms.forEach(a => handled.add(a));
      if (principal === 'one') suff.push(lx); else prefixes.push({name: 'oxo', loc: lx});
    }
    for (const {n: y, o} of M.adj[x]) {
      if (handled.has(y) || blocked.has(y) || (M.mask && !M.mask.has(y))) continue;
      if (fixed && y === fixed.from) continue;
      if (o !== 1) return null;   // =CH2, =N … nicht vorgesehen
      const ht = hetGroup(y);
      const gy = groups.get(y);
      if (ht && ht === principal) { suff.push(lx); handled.add(y); continue; }
      if (gy && gy.type === principal && isRing) { suff.push(lx); handled.add(y); gy.atoms.forEach(a => handled.add(a)); continue; }
      // Amin als Hauptgruppe mit weiteren Resten am N: N-Substituenten
      const sub = nmSub(M, y, x, groups, hetGroup);
      prefixes.push({name: sub.name, loc: lx, complex: sub.complex});
      handled.add(y);
    }
  }
  // sekundäre/tertiäre Amine: Reste am N als N-Substituenten
  if (principal === 'amine') {
    for (const x of seq) for (const {n: N} of M.adj[x]) {
      if (M.el[N] !== 'N' || hetGroup(N) !== 'amine' || seq.includes(N)) continue;
      if (!suff.includes(pos.get(x))) continue;
      for (const {n: R} of M.adj[N]) {
        if (R === x) continue;
        const sub = nmSub(M, R, N, groups, hetGroup);
        prefixes.push({name: sub.name, loc: 'N', complex: sub.complex});
      }
    }
  }
  if (fixed && suff.length) return null;
  if (!fixed && principal && !suff.length) return null;
  // Name zusammensetzen
  const numP = prefixes.filter(p => p.loc !== 'N');
  const pLocs = prefixes.map(p => p.loc === 'N' ? 0 : p.loc).sort((a, b) => a - b);
  const alpha = [...prefixes].sort((a, b) => nmAlphaKey(a.name).localeCompare(nmAlphaKey(b.name)) || (a.loc === 'N' ? -1 : 0)).map(p => p.loc === 'N' ? 0 : p.loc);
  const mb = [...en, ...yn].sort((a, b) => a - b);
  const principalCount = suff.length;
  const key = [
    -principalCount,
    isRing ? 0 : 1,                       // Ring vor Kette (bei gleicher Zahl an Hauptgruppen)
    -L,                                   // längere Kette / größerer Ring
    -mb.length, -en.length,
    -prefixes.length,
    [...suff].sort((a, b) => a - b), mb, en.sort((a, b) => a - b), pLocs, alpha
  ];
  if (fixed) key.splice(0, 2);
  const name = nmAssemble(M, {isRing, ringInfo, L, seq, suff, principal, en: en.sort((a, b) => a - b), yn: yn.sort((a, b) => a - b), prefixes, fixed});
  if (name === null) return null;
  return {key, name, seq, pos};
}
const nmAlphaKey = n => n.replace(/^[\d,N'\-]+/, '').replace(/^\(/, '').replace(/^(di|tri|tetra|penta|hexa|hepta|octa|nona|deca|bis|tris|tetrakis)(?=[a-z(])/, '').replace(/[()\d,\-]/g, '');

// Substituent ab Atom y (gebunden an x)
function nmSub(M, y, x, groups, hetGroup) {
  groups = groups || new Map();
  const el = M.el[y];
  if (NM_HAL[el]) return {name: NM_HAL[el]};
  if (el === 'O') {
    if (M.deg(y) === 1) { if (M.h[y] === 1) return {name: 'hydroxy'}; nmFail(); }
    const r = M.adj[y].find(z => z.n !== x).n;
    const gr2 = groups.get(r);
    if (gr2 && gr2.type === 'ester' && gr2.o === y) return {name: nmAcyl(M, r, y, groups) + 'oxy', complex: true};
    const s = nmSub(M, r, y, groups, hetGroup);
    return nmOxy(s);
  }
  if (el === 'S') {
    if (M.deg(y) === 1 && M.h[y] === 1) return {name: 'sulfanyl'};
    const r = M.adj[y].find(z => z.n !== x).n;
    if (M.el[r] !== 'C') nmFail();
    const s = nmSub(M, r, y, groups, hetGroup);
    return {name: (s.complex ? nmWrap(s.name) : s.name) + 'sulfanyl', complex: true};
  }
  if (el === 'N') {
    const os = M.adj[y].filter(z => M.el[z.n] === 'O');
    if (os.length === 2 && M.deg(y) === 3) return {name: 'nitro'};
    if (M.deg(y) === 1 && M.h[y] === 2) return {name: 'amino'};
    if (M.adj[y].every(z => z.o === 1 && M.el[z.n] === 'C')) {
      const rs = M.adj[y].filter(z => z.n !== x).map(z => nmSub(M, z.n, y, groups, hetGroup));
      const parts = nmPrefixString(rs.map(s => ({...s, loc: null})));
      return {name: parts + 'amino', complex: true};
    }
    nmFail();
  }
  if (el !== 'C') nmFail();
  // Kohlenstoff-Rest
  const gr = groups.get(y);
  if (gr && gr.type === 'acid') return {name: 'carboxy'};
  if (gr && gr.type === 'nitrile') return {name: 'cyano'};
  if (gr && gr.type === 'al' && M.adj[y].filter(z => M.el[z.n] === 'C').length === 1 && M.h[y] === 1) return {name: 'formyl'};
  if (gr && gr.type === 'ester' && gr.o !== x) {
    const R = M.adj[gr.o].find(z => z.n !== y).n, ox = nmOxy(nmSub(M, R, gr.o, groups, hetGroup));
    return {name: (ox.complex ? nmWrap(ox.name) : ox.name) + 'carbonyl', complex: true};
  }
  if (gr && (gr.type === 'ester' || gr.type === 'amide')) nmFail();
  const core = nmCore(M, groups, {attach: y, from: x});
  return {name: core.name, complex: /[\d,]/.test(core.name)};
}

// Acylrest R–CO– (am Carbonyl-C c, Ester-O o): Acetyl, Propanoyl, Benzoyl …
function nmAcyl(M, c, o, groups) {
  const rest = M.adj[c].filter(z => z.n !== o && !(M.el[z.n] === 'O' && z.o === 2));
  if (!rest.length) return 'formyl';
  const r = rest[0].n;
  if (M.inRing(r)) { const k = M.ringOf.get(r)[0]; if (nmRingInfo(M, k).arom && M.rings[k].every(x => M.deg(x) === 2 || x === r)) return 'benzoyl'; nmFail(); }
  let n = 1, prev = c, cur = r;
  for (;;) { if (M.el[cur] !== 'C' || M.adj[cur].some(z => z.o !== 1)) nmFail(); n++; const nx = M.adj[cur].filter(z => z.n !== prev); if (!nx.length) break; if (nx.length > 1) nmFail(); prev = cur; cur = nx[0].n; }
  return n === 2 ? 'acetyl' : nmLower(NM_STEM[n]) + 'anoyl';
}
// Alkyl → Alkoxy: Methoxy, Ethoxy, Propoxy, Butoxy, Phenoxy (auch 2-Hydroxyethoxy), sonst …yloxy
function nmOxy(s) {
  if (/(meth|eth|prop|but|phen)yl$/.test(s.name)) return {name: s.name.replace(/yl$/, 'oxy'), complex: s.complex};
  return {name: (s.complex ? nmWrap(s.name) : s.name) + 'oxy', complex: true};
}
const nmWrap = n => n.includes('[') ? '{' + n + '}' : n.includes('(') ? '[' + n + ']' : '(' + n + ')';
function nmPrefixString(prefixes) {
  // gleiche Präfixe zusammenfassen, alphabetisch ordnen
  const byName = new Map();
  for (const p of prefixes) { if (!byName.has(p.name)) byName.set(p.name, {name: p.name, complex: p.complex, locs: []}); byName.get(p.name).locs.push(p.loc); }
  const list = [...byName.values()].sort((a, b) => nmAlphaKey(a.name).localeCompare(nmAlphaKey(b.name)));
  return list.map(p => {
    const k = p.locs.length;
    const locs = p.locs.filter(l => l !== null && l !== 0).sort((a, b) => (a === 'N' ? -1 : b === 'N' ? 1 : a - b));
    const mult = k > 1 ? (p.complex ? NM_MULTC[k] : NM_MULT[k]) : '';
    const body = p.complex ? nmWrap(p.name) : p.name;
    return (locs.length ? locs.join(',') + '-' : '') + mult + body;
  }).join('-');
}

function nmAssemble(M, o) {
  const {isRing, ringInfo, L, suff, principal, en, yn, prefixes, fixed} = o;
  if (L >= NM_STEM.length && !isRing) return null;
  const nFeat = suff.length + prefixes.filter(p => p.loc !== 'N').length + en.length + yn.length;
  let pre = prefixes;
  // Lokanten weglassen, wo eindeutig
  const terminalSuffix = ['acid', 'al', 'nitrile', 'amide'].includes(principal) && !isRing;
  let omitAll = false;
  if (fixed) omitAll = !isRing && L === 1;
  else if (isRing && nFeat <= 1) omitAll = true;
  else if (!isRing && L === 1) omitAll = true;
  else if (!isRing && L <= 2 && nFeat <= 1) omitAll = true;
  else if (!isRing && L === 2 && suff.length === 0 && (prefixes.length <= 1 || (new Set(prefixes.map(p => p.name)).size === 1 && prefixes.length === (en.length ? 4 : yn.length ? 2 : 6)))) omitAll = true;
  else if (!isRing && L === 2 && suff.length === 1 && prefixes.filter(p => p.loc !== 'N').length === 0) omitAll = true;
  else if (!isRing && terminalSuffix && L <= 2 && prefixes.length + en.length + yn.length === 0) omitAll = true;
  else if (!isRing && L === 3 && nFeat === 1 && (en.length + yn.length === 1 || (principal === 'one' && suff.length === 1))) omitAll = true;
  else if (!isRing && L === 4 && nFeat === 1 && principal === 'one' && suff.length === 1) omitAll = true;
  const pstr = nmPrefixString(omitAll ? pre.map(p => ({...p, loc: p.loc === 'N' ? 'N' : 0})) : pre);
  let parent;
  if (isRing && ringInfo.arom) {
    const s1 = suff.length === 1 && (omitAll || suff[0] === 1);
    const R = {acid: ['Benzoesäure', 'carbonsäure'], al: ['Benzaldehyd', 'carbaldehyd'], ol: ['Phenol', 'ol'], amine: ['Anilin', 'amin'], nitrile: ['Benzonitril', 'carbonitril'], amide: ['Benzamid', 'carboxamid'], thiol: ['Benzolthiol', 'thiol']};
    if (fixed) parent = 'phenyl';
    else if (!principal || !suff.length) parent = 'Benzol';
    else if (suff.length === 1 && R[principal]) parent = R[principal][0];
    else if (R[principal]) parent = 'Benzol-' + suff.join(',') + '-' + NM_MULT[suff.length] + R[principal][1];
    else return null;
    if (fixed && suff.length) return null;
    // Phenol/Anilin …: Hauptgruppe an C1 (Lokant entfällt)
  } else {
    const stem = (isRing ? 'Cyclo' + nmLower(NM_STEM[L]) : NM_STEM[L]) || null;
    if (!stem || (isRing && L < 3)) return null;
    const ringSfx = {acid: 'carbonsäure', al: 'carbaldehyd', nitrile: 'carbonitril', amide: 'carboxamid'};
    // Ungesättigtheit
    let core = stem;
    const omitMB = omitAll || (!isRing && L === 2);
    const loc = arr => omitMB ? '' : arr.join(',');
    const parts = [];
    if (!en.length && !yn.length) core = stem + 'an';
    else {
      const mult = en.length > 1 || yn.length > 1;
      core = stem + (mult ? 'a' : '');
      if (en.length) parts.push((omitMB ? '' : '-' + loc(en) + '-') + NM_MULT[en.length] + 'en');
      if (yn.length) parts.push((omitMB ? '' : '-' + loc(yn) + '-') + NM_MULT[yn.length] + 'in');
      core += parts.join('');
    }
    if (fixed && !isRing && L === 1 && prefixes.length === 1 && prefixes[0].name === 'phenyl') return 'benzyl';
    if (fixed) {
      // Substituentenname: Anknüpfung an C1 → -yl
      if (!en.length && !yn.length) parent = (isRing ? core.replace(/an$/, '') : stem) + 'yl';
      else if (!isRing && L === 2) parent = stem + (en.length ? 'enyl' : 'inyl');   // Ethenyl, Ethinyl
      else parent = core + '-1-yl';
      parent = nmLower(parent);
    } else if (!principal || !suff.length) parent = core;
    else {
      const sfx = isRing && ringSfx[principal] ? ringSfx[principal] : {acid: 'säure', al: 'al', nitrile: 'nitril', amide: 'amid', one: 'on', ol: 'ol', thiol: 'thiol', amine: 'amin'}[principal];
      const k = suff.length, mult = NM_MULT[k];
      const showLoc = !omitAll && !(terminalSuffix);
      // Bindestrich-Regeln: Lokant vor der Endung ("Propan-2-ol", "Ethan-1,2-diol")
      const sLoc = showLoc ? '-' + [...suff].sort((a, b) => a - b).join(',') + '-' : '';
      if (sLoc && /(en|in)$/.test(core)) parent = core + sLoc + mult + sfx;
      else parent = core + sLoc + mult + sfx;
      if (!showLoc && k > 1 && terminalSuffix) parent = core + mult + sfx;
    }
  }
  let name = pstr ? pstr + (/^[A-Z]/.test(parent) ? nmLower(parent) : parent) : parent;
  // Substituent mit Lokanten: in Klammern (wird vom Aufrufer gesetzt)
  if (!fixed) name = nmCap(name);
  return name;
}

// Vergleich zweier Graphen (Konstitution) über verfeinerte Atomklassen
function nmHash(g, withCip) {
  const M = nmGraph(g);
  let lab = M.el.map((e, i) => e + '|' + M.ch[i] + '|' + M.h[i] + '|' + M.adj[i].map(x => x.o).sort().join('') + (withCip ? '|' + (g.atoms[M.orig[i]].cip || '') : ''));
  for (let it = 0; it < 6; it++) lab = lab.map((l, i) => l + '[' + M.adj[i].map(x => x.o + ':' + lab[x.n]).sort().join(';') + ']').map(s => { let h = 0; for (let k = 0; k < s.length; k++) h = (h * 31 + s.charCodeAt(k)) | 0; return String(h); });
  return M.n + ':' + [...lab].sort().join(',');
}
// Name erzeugen und durch Zurückübersetzen prüfen
// E/Z einer C=C-Doppelbindung aus der gezeichneten Lage (pos: Index → {x, y}); null, wenn keine Isomerie
function ezOf(g, i, j, pos, hc, cip) {
  const top = (c, other) => {
    const all = [...g.nb(c).filter(n => n !== other)]; for (let k = 0; k < hc[c]; k++) all.push(-1);
    if (all.length !== 2) return null;
    const v = cip.cmpSub(c, all[0], all[1]); if (!v) return null;
    return v > 0 ? all[0] : all[1];
  };
  const ti = top(i, j), tj = top(j, i);
  if (ti === null || tj === null || ti < 0 || tj < 0 || !pos[i] || !pos[j] || !pos[ti] || !pos[tj]) return null;
  const side = (t, a) => { const dx = pos[j].x - pos[i].x, dy = pos[j].y - pos[i].y; return Math.sign(dx * (pos[t].y - pos[a].y) - dy * (pos[t].x - pos[a].x)); };
  const s1 = side(ti, i), s2 = side(tj, j);
  if (!s1 || !s2) return null;
  return s1 === s2 ? 'Z' : 'E';
}
// opt: {pos: Index → {x, y} für E/Z, cip: Map Index → 'R'/'S'}
function iupacName(g, opt) {
  let res = null;
  try { res = nmName(g); } catch (e) { res = null; }
  if (!res) return null;
  let name = res.name;
  // Gegenprobe mit dem Namensleser: liest er eine andere Struktur, ist der Name falsch → verwerfen
  let back = null;
  try { back = fromInput(name); } catch (e) { back = null; }
  if (back && nmHash(back.g) !== nmHash(g)) return null;
  // Stereodeskriptoren (nur wenn alle an der Stammkette liegen)
  if (opt && res.locs) {
    const desc = []; let off = false;
    if (opt.cip) for (const [i, d] of opt.cip) { if (!d) continue; const l = res.locs.get(i); if (l === undefined) off = true; else desc.push({l, d}); }
    if (opt.pos) {
      const hc = hCounts(g), cip = makeCip(g, hc), ringSet = new Set(findRings(g).flat());
      for (const b of g.bonds) {
        if (b.order !== 2 || g.atoms[b.a].el !== 'C' || g.atoms[b.b].el !== 'C' || (ringSet.has(b.a) && ringSet.has(b.b))) continue;
        const e = ezOf(g, b.a, b.b, opt.pos, hc, cip); if (!e) continue;
        const la = res.locs.get(b.a), lb = res.locs.get(b.b);
        if (la === undefined || lb === undefined) { off = true; continue; }
        desc.push({l: Math.min(la, lb), d: e});
      }
    }
    if (desc.length && !off) name = '(' + desc.sort((a, b) => a.l - b.l).map(x => x.l + x.d).join(',') + ')-' + name;
  }
  return name;
}

return {iupacName, iupacNameRaw: g => { try { const r = nmName(g); return r && r.name; } catch (e) { return 'EXC ' + e.message; } }, nmHash, parseSmiles, kekulize, layout, functionalGroups, peptideBackbone, SHORT_BOND, findRings, nameToGraph, fromInput, fromSmiles, normalizeName, mirrorText, condensed,
  G, hCounts, stereoCenters, cipFromFischer, cipFromVectors, fischerLayout, assignStereo, extractStereo, haworthLayout, cip3};
})();
if (typeof module !== 'undefined') module.exports = Chem;

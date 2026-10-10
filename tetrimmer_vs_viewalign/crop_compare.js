// Compare ViewAlign's global gap-% end trim (cluster.js getTrimBoundaries) with TEtrimmer's per-row
// remove_gap_columns + crop_end_by_divergence (re-implemented from the Nature Commun. 2025 Methods text).
const fs = require('fs');
function readFa(p) {
  const out = []; let h = null, s = [];
  for (const l of fs.readFileSync(p, 'utf8').split(/\r?\n/)) {
    if (l[0] === '>') { if (h !== null) out.push({ h, seq: s.join('') }); h = l.slice(1); s = []; } else s.push(l.trim());
  }
  if (h !== null) out.push({ h, seq: s.join('') });
  return out;
}
const isG = c => c === '-' || c === '.';

function viewalignTrim(seqs, W = 15, L = 0.5, R = 0.8) {   // verbatim logic of cluster.js getTrimBoundaries
  const n = seqs.length, len = seqs[0].seq.length;
  const gaps = c => { let g = 0; for (const s of seqs) if (isG(s.seq[c])) g++; return g; };
  let left = -1, right = len, wg = 0, wc = 0;
  for (let c = 0; c < len; c++) {
    wc++; wg += gaps(c);
    if (wc > W) { wg -= gaps(c - W); wc--; }
    if (wg / (wc * n) > L) left = c; else break;
  }
  wg = 0; wc = 0; const arr = [];
  for (let c = len - 1; c >= 0; c--) {
    const g = gaps(c); arr.push(g); wc++; wg += g;
    if (wc > W) { wg -= arr.shift(); wc--; }
    if (wg / (wc * n) > R) right = c; else break;
  }
  return { left, right };
}

function colStats(seqs) {                    // per column: counts of ACGT, gap fraction, major fraction
  const len = seqs[0].seq.length, n = seqs.length, st = [];
  for (let c = 0; c < len; c++) {
    const cnt = { A: 0, C: 0, G: 0, T: 0 }; let g = 0;
    for (const s of seqs) { const ch = s.seq[c].toUpperCase(); if (isG(ch)) g++; else if (ch === 'U') cnt.T++; else if (cnt[ch] !== undefined) cnt[ch]++; }
    const tot = cnt.A + cnt.C + cnt.G + cnt.T, mx = Math.max(cnt.A, cnt.C, cnt.G, cnt.T);
    st.push({ cnt, tot, gapFrac: g / n, major: tot ? mx / tot : 0 });
  }
  return st;
}

// TEtrimmer remove_gap_columns: drop columns >80 % gap or <5 nt; then columns with 40-80 % gap and major <70 %
function tetGapCols(seqs) {
  const st = colStats(seqs), keep = [];
  for (let c = 0; c < st.length; c++) {
    const s = st[c];
    if (s.gapFrac > 0.8 || s.tot < 5) continue;
    if (s.gapFrac >= 0.4 && s.major < 0.7) continue;
    keep.push(c);
  }
  return seqs.map(x => ({ h: x.h, seq: keep.map(c => x.seq[c]).join('') }));
}

// crop_end_by_divergence: proportion matrix (cols with <5 nt -> 0; gaps -> 0); per row, a window slides in from each end
// until its mean proportion > thr; everything before the stop is deleted (turned to gaps). Two passes: (40, 0.7), (4, 1.0).
function cropDiv(seqs, win, thr) {
  const st = colStats(seqs), len = seqs[0].seq.length;
  return seqs.map(x => {
    const p = new Float64Array(len);
    for (let c = 0; c < len; c++) {
      const ch = x.seq[c].toUpperCase(); if (isG(ch) || st[c].tot < 5) continue;
      const k = ch === 'U' ? 'T' : ch; p[c] = st[c].cnt[k] !== undefined ? st[c].cnt[k] / st[c].tot : 0;
    }
    const pos = []; for (let c = 0; c < len; c++) pos.push(c);          // window measured over alignment columns of this row
    const mean = (a, b) => { let s = 0; for (let i = a; i < b; i++) s += p[i]; return s / (b - a); };
    const w = Math.min(win, len);
    let a = 0; while (a + w <= len && mean(a, a + w) < thr) a++;
    let b = len; while (b - w >= 0 && mean(b - w, b) < thr) b--;
    if (a >= b) return { h: x.h, seq: '-'.repeat(len) };
    return { h: x.h, seq: '-'.repeat(a) + x.seq.slice(a, b) + '-'.repeat(len - b) };
  });
}

function summarize(name, seqs) {
  const n = seqs.length, len = seqs[0].seq.length;
  const nt0 = seqs.reduce((s, x) => s + [...x.seq].filter(c => !isG(c)).length, 0);
  const v = viewalignTrim(seqs); const vk = Math.max(0, v.right - (v.left + 1));
  const vkept = seqs.reduce((s, x) => s + [...x.seq.slice(v.left + 1, v.right)].filter(c => !isG(c)).length, 0);
  let t = tetGapCols(seqs); const tcols = t[0].seq.length;
  t = cropDiv(t, 40, 0.7); const p1 = t.reduce((s, x) => s + [...x.seq].filter(c => !isG(c)).length, 0); console.log(`  TET pass 1 only (window 40, 0.7): residues kept ${(100 * p1 / nt0).toFixed(1)} %`); t = cropDiv(t, 4, 1.0);
  const tkept = t.reduce((s, x) => s + [...x.seq].filter(c => !isG(c)).length, 0);
  // quality of what is kept: mean major-allele fraction of retained residues (own base's column proportion), over the original alignment for VA, over the gap-cleaned one for TET
  const qual = (arr) => { const st = colStats(arr); let s = 0, m = 0; for (const x of arr) for (let c = 0; c < x.seq.length; c++) { const ch = x.seq[c].toUpperCase(); if (isG(ch) || !st[c].tot) continue; s += (st[c].cnt[ch] || 0) / st[c].tot; m++; } return s / m; };
  const vaSeqs = seqs.map(x => ({ h: x.h, seq: x.seq.slice(v.left + 1, v.right) }));
  console.log(`${name}: ${n} rows x ${len} cols, ${nt0} residues`);
  console.log(`  ViewAlign trim (global, gap %, 50/80/15): cut left ${v.left + 1}, right ${len - v.right} cols -> ${vk} cols, residues kept ${(100 * vkept / nt0).toFixed(1)} %, mean column agreement of kept residues ${qual(vaSeqs).toFixed(3)}`);
  console.log(`  TEtrimmer gap-columns: ${len} -> ${tcols} cols; + per-row divergence crop: residues kept ${(100 * tkept / nt0).toFixed(1)} %, mean column agreement ${qual(t).toFixed(3)}`);
  const lens = t.map(x => [...x.seq].filter(c => !isG(c)).length).sort((a, b) => a - b);
  console.log(`  TET per-row kept length quartiles: ${lens[0]} / ${lens[Math.floor(n / 4)]} / ${lens[Math.floor(n / 2)]} / ${lens[Math.floor(3 * n / 4)]} / ${lens[n - 1]}`);
}
for (const f of process.argv.slice(2)) summarize(f, readFa(f));

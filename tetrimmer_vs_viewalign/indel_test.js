const fs = require('fs'), R = 'C:/work/MSA-viewer-new/', Peel = require(R + 'peel.js'), { ari } = require(R + 'tests/kmer/metrics.js');
const DIR = 'C:/work/SINE_discriminator/site/alignments/';
const rd = f => fs.readFileSync(DIR + f, 'utf8').split('>').slice(1).map(b => { const l = b.split('\n'); return { header: l[0].trim(), seq: l.slice(1).join('').replace(/\s/g, '') }; });
const num = h => (h.match(/input_0*(\d+)/) || [])[1];
const sets = {
  ccr: () => { const all = rd('CURATE__ccr__subfam608_g1-g7_accr.aln.fa').filter(r => /^input_/.test(r.header));
    const grp = new Map(fs.readFileSync(DIR + 'ccr_groups.tsv', 'utf8').trim().split('\n').map(l => l.split('\t')).map(([g, n]) => [n.trim(), g]));
    return { seqs: all, truth: all.map(s => grp.get(s.header) || null) }; },
  oma: () => { const all = rd('CURATE__oma__subfam600_23seeds.aln.fa').filter(r => /^input_/.test(r.header)); const m = new Map();
    fs.readFileSync(DIR + 'oma_groups_final.tsv', 'utf8').trim().split('\n').forEach(l => { const [g, ids] = l.split('\t'); ids.trim().split(/\s+/).forEach(x => m.set(String(+x), g)); });
    return { seqs: all, truth: all.map(s => m.get(String(+num(s.header))) || null) }; } };
const variants = {
  'baseline (columns)': {},
  'indel run=1, w1': { indelWeight: 1 },
  'indel run=1, w2': { indelWeight: 2 },
  'indel run=1, w3': { indelWeight: 3 },
  'w1, minDiag 3': { indelWeight: 1, minDiag: 3 },
  'w2, minDiag 3': { indelWeight: 2, minDiag: 3 },
  'w1, refine 2 w1': { indelWeight: 1, refineIndelWeight: 1 },
  'refine w1 only': { refineIndelWeight: 1 },
};
for (const [name, load] of Object.entries(sets)) {
  const { seqs, truth } = load(), n = seqs.length, idx = truth.map((t, i) => t ? i : -1).filter(i => i >= 0), tr = idx.map(i => truth[i]);
  console.log(`\n${name}: ${n} chunks, ${idx.length} labelled`);
  for (const [vn, o] of Object.entries(variants)) {
    const r = Peel.peel(seqs, Object.assign({ metric: 'pdist' }, Peel.defaults, o));
    const lab = new Array(n).fill(null); r.groups.forEach((g, gi) => g.forEach(i => { lab[i] = 'g' + gi; })); r.unassigned.forEach(i => { lab[i] = 'u' + i; });
    // his groups split into >=2 pieces of size>=2 (calibration cases)
    const pieceOf = new Map(); r.groups.forEach((g, k) => g.forEach(i => pieceOf.set(i, k)));
    let cases = 0; for (const g of new Set(tr)) { if (g === 'ancient78') continue; const by = new Map(); idx.forEach(i => { if (truth[i] === g && pieceOf.has(i)) { const k = pieceOf.get(i); by.set(k, (by.get(k) || 0) + 1); } }); const big = [...by.values()].filter(x => x >= 2).length; if (by.size >= 2) cases += by.size - 1; }
    console.log(`  ${vn.padEnd(22)} groups ${String(r.groups.length).padStart(3)}  ARI ${ari(tr, idx.map(i => lab[i])).toFixed(3)}  unassigned ${r.unassigned.length}  split-pieces ${cases}`);
  }
}

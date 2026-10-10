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

module.exports={sets,Peel,ari};

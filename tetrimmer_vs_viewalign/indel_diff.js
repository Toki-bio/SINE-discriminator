const {sets,Peel}=require('./loader.js');
const {seqs,truth}=sets.oma(); const n=seqs.length;
const run=w=>Peel.peel(seqs,Object.assign({metric:'pdist'},Peel.defaults,{refineIndelWeight:w}));
const a=run(1), b=run(2);
const pieces=r=>{const m=new Map();r.groups.forEach((g,k)=>g.forEach(i=>m.set(i,k)));return m};
const pa=pieces(a), pb=pieces(b);
for(const g of new Set(truth.filter(Boolean))){ if(g==='ancient78')continue;
  const mem=truth.map((t,i)=>t===g?i:-1).filter(i=>i>=0);
  const cnt=p=>{const c=new Map();mem.forEach(i=>{if(p.has(i))c.set(p.get(i),(c.get(p.get(i))||0)+1)});return [...c.values()].sort((x,y)=>y-x).join('+')};
  const ca=cnt(pa), cb=cnt(pb); if(ca!==cb) console.log(g, 'his size', mem.length, '| refine w1:', ca, '| w2:', cb);
}

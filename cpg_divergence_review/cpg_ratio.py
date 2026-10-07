import sys,glob,os,statistics as st
from plate_stats import rd
TS={('C','T'),('T','C'),('A','G'),('G','A')}
def run(path):
    R=rd(path)
    # plate format: row 1 = extended consensus, row 2 = consensus; fall back to row 1 if 'extended' absent
    cons=R[1][1] if 'extended' in R[0][0].lower() else R[0][1]
    el=[i for i,c in enumerate(cons) if c.isupper() and c!='-']
    # ungapped consensus order of element columns; CpG = consensus C followed by G, or G preceded by C (in consensus without gaps)
    seq=[(i,cons[i].upper()) for i in el if cons[i]!='-']
    cpg=set()
    for a in range(len(seq)-1):
        if seq[a][1]=='C' and seq[a+1][1]=='G': cpg.add(seq[a][0]); cpg.add(seq[a+1][0])
    rows=R[2:] if 'extended' in R[0][0].lower() else R[1:]
    d_cpg=[];d_non=[];ts_cpg=[];ts_non=[]
    for n,s in rows:
        c=m=cs=ms=ct=mt=0
        for i,_ in seq:
            x=s[i].upper()
            if x in '-N' : continue
            ref=cons[i].upper()
            iscpg=i in cpg
            mism=(x!=ref)
            tr=(ref,x) in TS
            if iscpg:
                c+=1; m+=mism; ct+=(mism and tr)
            else:
                cs+=1; ms+=mism; mt+=(mism and tr)
        if c>=8 and cs>=40:
            d_cpg.append(m/c); d_non.append(ms/cs)
            ts_cpg.append(ct/m if m else float('nan')); ts_non.append(mt/ms if ms else float('nan'))
    n=len(d_cpg)
    if n<10: return None
    mc,mn=st.mean(d_cpg),st.mean(d_non)
    cpgn=len(cpg)
    return (n,cpgn,len(seq),mc,mn,(mc/mn if mn else float('nan')),
            st.mean([x for x in ts_cpg if x==x] or [float('nan')]),st.mean([x for x in ts_non if x==x] or [float('nan')]))
print('%-28s %5s %9s %10s %10s %7s %9s %9s'%('plate','copies','CpG/total','CpG dens','nonCpG','ratio','Ts@CpG','Ts@non'))
out=[]
for f in sorted(glob.glob('rsiplates/*_top100.aln.fa')):
    r=run(f)
    name=os.path.basename(f).replace('_top100.aln.fa','').replace('rsi_','rsi:')
    if r: out.append((name,r)); print('%-28s %5d %4d/%-4d %10.3f %10.3f %7.1f %9.2f %9.2f'%(name,r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7]))

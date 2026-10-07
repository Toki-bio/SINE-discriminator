import sys,itertools,collections
from plate_stats import rd
def rc(s): return s[::-1].translate(str.maketrans('ACGT','TGCA'))
def flanks(path):
    R=rd(path)[2:]; out=[]
    for n,s in R:
        raw=s.replace('-','')
        up=[i for i,c in enumerate(raw) if c.isupper()]
        if not up: continue
        out.append((n,raw[:up[0]].upper(),raw[up[-1]+1:].upper(),raw[up[0]:up[-1]+1].upper()))
    return out
K=16
def kset(s): return {s[i:i+K] for i in range(len(s)-K+1)} if len(s)>=K else set()
def shared(a,b):
    ka=kset(a[1])|kset(a[2]); 
    for x in (b[1],b[2]):
        if kset(x)&ka or kset(rc(x))&ka: return True
    return False
for fam in ['r4_32seqs','r2_3seqs','r3_58seqs']:
    S={sp:flanks('cross/%s_%s.fa'%(sp,fam)) for sp in ['rsi','rda','rre']}
    print('=== %s  copies in plates: %s'%(fam,{k:len(v) for k,v in S.items()}))
    for a,b in [('rsi','rda'),('rsi','rre'),('rda','rre')]:
        A,B=S[a],S[b]
        hits=sum(1 for x in A if any(shared(x,y) for y in B))
        print('   %s copies sharing a 16-mer in a flank with a %s copy (either strand): %d of %d (%.0f%%)'%(a,b,hits,len(A),100*hits/max(1,len(A))))
    # random control: families that should not be shared (compare r4 vs r3 plate across species)

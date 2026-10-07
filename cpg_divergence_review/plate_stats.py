import sys
def rd(f):
    r=[];n=None
    for l in open(f):
        l=l.rstrip('\n')
        if l.startswith('>'): n=l[1:]; r.append([n,[]])
        elif r: r[-1][1].append(l)
    return [(a,''.join(b)) for a,b in r]
def stats(path):
    R=rd(path); cons=R[1][1]; elem=[i for i,c in enumerate(cons) if c.isupper()]
    rows=[]
    for k,(n,s) in enumerate(R[2:],start=1):
        up=sum(1 for c in s if c.isupper()); low=sum(1 for c in s if c.islower())
        m=t=0
        for i in elem:
            c=s[i]
            if c=='-': continue
            t+=1; m+=(c.upper()==cons[i].upper())
        rows.append((k,n,up,low,(m/t if t else 0),t))
    return R,rows
if __name__=='__main__':
    R,rows=stats(sys.argv[1]); n=len(rows)
    print('consensus element columns:',sum(1 for c in R[1][1] if c.isupper()),'| copies',n)
    print('rank name                                         element_bp flank_bp identity_on_cons_cols (columns used)')
    for r in rows[:3]+rows[-int(sys.argv[2]):]:
        print('%4d %-45s %6d %6d   %.2f (%d)'%(r[0],r[1][:45],r[2],r[3],r[4],r[5]))

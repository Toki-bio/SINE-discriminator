#!/usr/bin/env python3
# Test of ways to handle SubFam's last (short) chunk, on real loci with known subfamily labels. Python 3.5 compatible (KIT).
#   python3 last_chunk_test.py <species> <population.fa> <labels.tsv> <outdir> <N> <repeats>
# population.fa: headers = locus id; labels.tsv: locus id <TAB> label.
# One replicate: sample N loci, order them as SubFam does (mafft --retree 0 --reorder), cut 50-locus chunks; r = N mod 50 loci are left over.
# Options for those r loci:  A drop (SubFam since 2026-09-05) | B keep short, plurality scaled | C last 50 loci (overlap) |
#                            D merge into the previous chunk | E balanced split (all chunks 49-50, no loss, no overlap)
# For each: how many loci are lost / double counted, purity of the end chunk, informative fraction of its consensus.
import sys, os, random, subprocess, collections

species, popfa, labtsv, out, N, REP = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
os.makedirs(out, exist_ok=True)
THREADS = os.environ.get('THREADS', '16')

def read_fa(path):
    names, seqs, cur = [], [], None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            names.append(line[1:].split()[0]); seqs.append([])
        elif names:
            seqs[-1].append(line)
    return names, [''.join(s) for s in seqs]

def write_fa(path, names, seqs):
    with open(path, 'w') as f:
        for n, s in zip(names, seqs):
            f.write('>%s\n%s\n' % (n, s))

names, seqs = read_fa(popfa)
seq_of = dict(zip(names, seqs))
label = dict(l.rstrip('\n').split('\t')[:2] for l in open(labtsv))
keep = set(sys.argv[7].split(',')) if len(sys.argv) > 7 else None
pop = [n for n in names if n in label and (keep is None or label[n] in keep)]
overall = collections.Counter(label[n] for n in pop)

def run(cmd, stdin_path=None, stdout_path=None):
    with open(stdout_path, 'w') if stdout_path else open(os.devnull, 'w') as fo:
        fi = open(stdin_path) if stdin_path else None
        subprocess.check_call(cmd, shell=True, stdin=fi, stdout=fo, stderr=subprocess.DEVNULL)

def chunk_consensus(ids, tag, plurality):
    """SubFam's per-chunk step: mafft --nuc --reorder, cons -plurality, N -> gap. Returns (informative fraction, purity, majority label)."""
    if len(ids) < 2:           # one locus: mafft needs two; its 'consensus' is the locus itself if the plurality is 1
        seq = seq_of[ids[0]]
        return ((sum(1 for ch in seq if ch not in '-nN') / float(len(seq))) if plurality <= 1 else 0.0), 1.0, label[ids[0]], (seq if plurality <= 1 else '-' * len(seq))
    fa = os.path.join(out, 'tmp_%s.fa' % tag)
    write_fa(fa, ids, [seq_of[i] for i in ids])
    al = fa + '.al'
    run('mafft --thread %s --nuc --reorder --quiet %s' % (THREADS, fa), stdout_path=al)
    con = fa + '.cons'
    run('cons -plurality %d -name x -filter < %s' % (plurality, al), stdout_path=con)
    _, cs = read_fa(con)
    c = cs[0].replace('N', '-').replace('n', '-') if cs else ''
    info = (sum(1 for ch in c if ch != '-') / float(len(c))) if c else 0.0
    cnt = collections.Counter(label[i] for i in ids)
    maj, k = cnt.most_common(1)[0]
    for p in (fa, al, con):
        try: os.remove(p)
        except OSError: pass
    return info, k / float(len(ids)), maj, c

rows = []
for rep in range(REP):
    rnd = random.Random(1000 * N + rep)
    sample = rnd.sample(pop, N)
    sfa = os.path.join(out, 'sample_%d_%d.fa' % (N, rep))
    write_fa(sfa, sample, [seq_of[i] for i in sample])
    ofa = sfa + '.ordered'
    run('mafft --thread %s --threadtb %s --threadit %s --nuc --quiet --retree 0 --reorder %s' % (THREADS, THREADS, THREADS, sfa), stdout_path=ofa)
    order, _ = read_fa(ofa)
    os.remove(sfa); os.remove(ofa)
    full, r = len(order) // 50, len(order) % 50
    tail = order[full * 50:]
    lab_tail = collections.Counter(label[i] for i in tail)
    # purity of the ordinary chunks (reference level)
    ref = []
    for c in range(min(full, 6)):
        cnt = collections.Counter(label[i] for i in order[c * 50:(c + 1) * 50]); ref.append(cnt.most_common(1)[0][1] / 50.0)
    # groups that exist in the sample but would have no locus left if the tail is dropped
    kept_labels = set(label[i] for i in order[:full * 50])
    lost_groups = [g for g in set(label[i] for i in order) if g not in kept_labels]
    base = dict(species=species, N=N, rep=rep, r=r, ref_purity=sum(ref) / len(ref) if ref else 0)
    dump = {'A_drop': {'members': [label[i] for i in tail]}}
    rows.append(dict(base, option='A_drop', lost=r, double=0, end_n=0, purity='', info='', lost_groups=len(lost_groups), tail_labels=','.join('%s:%d' % kv for kv in sorted(lab_tail.items()))))
    if r == 0:
        continue
    def add(opt, ids, plur, lost, double):
        info, pur, maj, cons_ = chunk_consensus(ids, '%s_%d_%d' % (opt, N, rep), plur)
        dump[opt] = {'members': [label[i] for i in ids], 'plurality': plur, 'info': round(info, 3), 'purity': round(pur, 3), 'consensus': cons_}
        rows.append(dict(base, option=opt, lost=lost, double=double, end_n=len(ids), purity=round(pur, 3), info=round(info, 3), lost_groups=0, tail_labels=maj))
    add('B_keep_short_plur18', tail, 18, 0, 0)                         # what SubFam did before 2026-09-05
    add('B2_keep_short_scaled', tail, max(1, int(round(0.36 * r))), 0, 0)
    add('C_last50_overlap', order[-50:], 18, 0, 50 - r)
    add('D_merge_previous', order[(full - 1) * 50:], max(1, int(round(0.36 * (50 + r)))), 0, 0)
    # E: balanced split of everything into `full + 1` chunks of 49-50 loci (no loss, no overlap)
    n_ch = full + 1; sizes = [len(order) // n_ch + (1 if k < len(order) % n_ch else 0) for k in range(n_ch)]
    pos, pur_all, info_last = 0, [], None
    for k, s in enumerate(sizes):
        ids = order[pos:pos + s]; pos += s
        cnt = collections.Counter(label[i] for i in ids); pur_all.append(cnt.most_common(1)[0][1] / float(s))
        if k == n_ch - 1:
            info_last = chunk_consensus(ids, 'E_%d_%d' % (N, rep), max(1, int(round(0.36 * s))))
            dump['E_balanced'] = {'members': [label[i] for i in ids], 'plurality': max(1, int(round(0.36 * s))), 'info': round(info_last[0], 3), 'purity': round(info_last[1], 3), 'consensus': info_last[3]}
    rows.append(dict(base, option='E_balanced', lost=0, double=0, end_n=sizes[-1], purity=round(info_last[1], 3), info=round(info_last[0], 3), lost_groups=0,
                     tail_labels='sizes %d-%d, mean purity of all chunks %.3f' % (min(sizes), max(sizes), sum(pur_all) / len(pur_all))))
    if os.environ.get('DUMP'):
        import json
        dump['_order'] = [label[i] for i in order]
        dump['_full'] = full; dump['_r'] = r; dump['_N'] = len(order); dump['_sizesE'] = sizes
        json.dump(dump, open(os.path.join(out, 'example_%s_N%d_rep%d.json' % (species, N, rep)), 'w'))
    sys.stderr.write('%s N=%d rep=%d r=%d done\n' % (species, N, rep, r))

cols = ['species', 'N', 'rep', 'r', 'option', 'lost', 'double', 'end_n', 'purity', 'info', 'lost_groups', 'ref_purity', 'tail_labels']
with open(os.path.join(out, '%s_N%d.tsv' % (species, N)), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r_ in rows:
        f.write('\t'.join(str(r_[c]) for c in cols) + '\n')
print('wrote', os.path.join(out, '%s_N%d.tsv' % (species, N)))

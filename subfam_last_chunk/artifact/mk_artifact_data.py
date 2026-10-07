import json, glob, os
D = 'C:/work/MSA-viewer/scratch/lastchunk/ex/'
pick = [
    ('mix7_N2017_rep0', 'All seven groups, 2,017 loci', 'The full ccr sample: 17 loci are left over after 40 chunks.'),
    ('mix3_N1033_rep1', 'Three groups, 1,033 loci', 'g4, g6 and g7 only: the 33 leftover loci are all g7.'),
    ('tight_N317_rep0', 'Two tight groups, 317 loci', 'g5 and new13 only: all 17 leftover loci are new13.'),
    ('small7_N117_rep0', 'A small run, 117 loci', 'The size of a per-family run: 17 of 117 loci (15 %) are dropped.'),
]
OPTS = ['A_drop', 'B_keep_short_plur18', 'B2_keep_short_scaled', 'C_last50_overlap', 'D_merge_previous', 'E_balanced']
code = {'g2': 'a', 'g3': 'b', 'g4': 'c', 'g5': 'd', 'g6': 'e', 'g7': 'f', 'new13': 'g'}
out = {'examples': []}
for key, title, sub in pick:
    d = json.load(open(D + 'example_' + key + '.json'))
    ex = {'title': title, 'sub': sub, 'N': d['_N'], 'r': d['_r'], 'full': d['_full'],
          'order': ''.join(code[x] for x in d['_order']), 'sizesE': d['_sizesE'], 'opts': {}}
    for o in OPTS:
        v = d[o]
        ex['opts'][o] = {'members': ''.join(code[x] for x in v['members']), 'info': v.get('info'), 'purity': v.get('purity'),
                         'plurality': v.get('plurality'), 'consensus': (v.get('consensus') or '')[:150]}
    out['examples'].append(ex)
json.dump(out, open('C:/work/MSA-viewer/scratch/lastchunk/artifact_data.json', 'w'), separators=(',', ':'))
print(len(json.dumps(out)), 'bytes;', [e['title'] for e in out['examples']])
for e in out['examples']:
    o = e['opts']['B2_keep_short_scaled']; print(e['N'], e['r'], 'B2 consensus start:', o['consensus'][:60])

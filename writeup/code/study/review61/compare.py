# compare own per-field sums (bins of u, width 1/96) with the authors' stored per-field sums
import sys, numpy as np
R = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/'
for r2 in map(int, sys.argv[1:]):
    m = np.load(R + 'review_162/mine_r%d.npz' % r2); a = np.load(R + 'quartic/data/hold_r%d.npz' % r2)
    import cypari2; pari = cypari2.Pari()
    rows = [l.split() for l in open(R + 'quartic/data/hold_r%d_fields.txt' % r2) if not l.startswith('#')]
    rows = sorted(set((int(D), f) for D, f in rows), key=lambda r: (abs(r[0]), r[1]))      # the authors' order (qcounts.field_list)
    assert all(int(a['D'][i]) == rows[i][0] for i in range(len(rows)))
    pos = {str(pari.polredabs(pari(f))).replace(' ', ''): i for i, (D, f) in enumerate(rows)}
    mine = [l.split() for l in open(R + 'review_162/fields_r%d.txt' % r2)]
    ia = np.array([pos[mine[j][1]] for j in m['idx']])
    assert all(int(a['D'][ia[k]]) == int(m['D'][k]) for k in range(len(ia)))
    dn = np.abs(m['n96'] - a['num'][ia]); dd = np.abs(m['d96'] - a['den'][ia])
    bad = np.nonzero((dn.max(1) > 1e-6) | (dd.max(1) > 1e-6))[0]
    print('r2=%d: %d fields compared; max |num diff| %.2e, max |den diff| %.2e; fields with any difference: %d %s'
          % (r2, len(ia), dn.max(), dd.max(), len(bad), [int(m['D'][b]) for b in bad[:5]]))

import numpy as np, sys
h = np.load(sys.argv[1]); m = np.load(sys.argv[2])
i0, i1 = int(m['i0']), int(m['i1'])
sl = slice(i0 + 1, i1 + 1)
print('cells', i0 + 1, '..', i1, 'p in (', m['X'][i0], ',', m['X'][i1], ']')
for k in ['N2', 'N4', 'NSG', 'NT026', 'NT046']:
    a, b = h[k][sl], m[k][sl]
    print(k, 'hold total', int(a.sum()), 'mine', int(b.sum()), 'cells differing', int((a != b).sum()))
for k in ['P2', 'P4', 'SG', 'T026', 'T046']:
    a, b = h[k][sl], m[k][sl]
    rel = np.abs(a - b) / np.abs(a)
    print(k, 'hold total %.6f mine %.6f' % (a.sum(), b.sum()), 'max rel diff per cell %.2e' % rel.max())
# residue classes: hold c5[:, r] = sum ln p over p = r mod 5 ?
for q in [3, 4, 5, 7, 12, 15, 20, 21, 28]:
    key = 'c%d' % q
    if key not in h.files: continue
    a = h[key][sl]
    thc = m['thc'][sl]
    b = np.stack([thc[:, [c for c in range(420) if c % q == r]].sum(1) for r in range(q)], 1)
    rel = np.abs(a - b).max() / np.abs(a).max()
    print(key, 'max abs diff / max %.2e' % rel)

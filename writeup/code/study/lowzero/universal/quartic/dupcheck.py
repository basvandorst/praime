import cypari2, os
from collections import defaultdict
pari = cypari2.Pari(); pari.allocatemem(200_000_000); H = os.path.dirname(os.path.abspath(__file__)); bad = 0; pairs = 0
for r2 in (0, 1, 2):
    g = defaultdict(list)
    for l in open(os.path.join(H, 'data', 'disc_r%d_fields.txt' % r2)):
        if not l.startswith('#'): D, f = l.split(); g[int(D)].append(f)
    for D, fs in g.items():
        for i in range(len(fs)):
            for j in range(i + 1, len(fs)):
                pairs += 1; bad += int(pari.nfisisom(pari(fs[i]), pari(fs[j])) != 0)
    print('r2=%d: %d fields, %d discriminants' % (r2, sum(len(v) for v in g.values()), len(g)))
print('pairs with equal discriminant: %d, isomorphic pairs: %d' % (pairs, bad))

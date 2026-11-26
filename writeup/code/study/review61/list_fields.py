# review 162, task 1a: list S4 quartic fields 20000 <= |D| < 60000 by r2 with one nflist call per signature (other chunking than the
# authors' 8 chunks of 5000), filter by |D|, polredabs, check S4 and no isomorphic pair; compare with the authors' field files.
import time, sys
import cypari2
pari = cypari2.Pari(); pari.allocatemem(300_000_000); pari.default('parisizemax', 600_000_000)
SRC = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/quartic/data/hold_r%d_fields.txt'
OUT = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_162/fields_r%d.txt'
LO, HI = 20000, 60000
for r2 in (0, 1, 2):
    t0 = time.time()
    fl = pari('nflist("S4", [%d, %d], %d)' % (LO, HI - 1, r2))
    mine = {}
    nret = len(fl); bad_sig = 0; bad_gal = 0
    for f in fl:
        D = int(pari.nfdisc(f))
        if not (LO <= abs(D) < HI): continue
        g = pari.polredabs(f)
        if int(pari.polsturm(g)) != 4 - 2 * r2: bad_sig += 1
        if int(pari.polgalois(g)[0]) != 24: bad_gal += 1
        key = str(g).replace(' ', '')
        if key in mine: print('DUPLICATE (isomorphic) field', key, D)
        mine[key] = D
    sgn_ok = all((D < 0) == (r2 == 1) for D in mine.values())
    # authors' list
    auth = {}
    for line in open(SRC % r2):
        if line.startswith('#'): continue
        D, fs = line.split(); g = str(pari.polredabs(pari(fs))).replace(' ', '')
        auth.setdefault(g, []).append(int(D))
    nauth_rows = sum(len(v) for v in auth.values()); auth_dup = sum(1 for v in auth.values() if len(v) > 1)
    only_mine = set(mine) - set(auth); only_auth = set(auth) - set(mine)
    disc_mismatch = sum(1 for k in set(mine) & set(auth) if auth[k][0] != mine[k])
    print('r2=%d: nflist returned %d, in range %d (distinct after polredabs %d), sign ok %s, wrong signature %d, not S4 %d; '
          'authors: %d rows, %d distinct fields, %d fields listed twice; only mine %d, only authors %d, disc mismatch %d; %.1f s'
          % (r2, nret, len(mine), len(set(mine)), sgn_ok, bad_sig, bad_gal, nauth_rows, len(auth), auth_dup, len(only_mine), len(only_auth),
             disc_mismatch, time.time() - t0), flush=True)
    with open(OUT % r2, 'w') as o:
        for k, D in sorted(mine.items(), key=lambda kv: (abs(kv[1]), kv[0])): o.write('%d %s\n' % (D, k))

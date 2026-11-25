import sys, os, json, time
REPO = '/home/bas/praime/writeup/code/study/lowzero/universal/blind4'
sys.path.insert(0, REPO); sys.path.insert(0, os.path.dirname(REPO))
import cypari2
from bz4_lib import params
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
G = json.load(open(os.path.join(REPO, 'genus2.json'))); TS = json.load(open(os.path.join(REPO, 'testset4.json')))
drawn = set((g[0], g[1], g[2]) for g in TS['g2'])
todo = [r for r in G['primitive'] if (r[0], r[1], r[2]) not in drawn]
for idx in (13, 12, 14):
    N, f, h, m2, feq = todo[idx]; print(idx, N, f, h, m2, feq)
    L = pari('lfungenus2([%s,%s])' % (f, h)); print('  params', params(pari, L), 'lfunparams', pari.lfunparams(L))
for T in (40, 80, 120, 153):
    t0 = time.time()
    try:
        N, f, h, m2, feq = todo[13]; L = pari('lfungenus2([%s,%s])' % (f, h))
        Li = pari.lfuninit(L, [1, T]); print('lfuninit [1,%d] ok %.0f s, stack used %.0f MB' % (T, time.time() - t0, (pari.stacksize() - pari.getavma() if hasattr(pari, 'getavma') else 0) / 1e6)); del Li
    except cypari2.PariError as e:
        print('lfuninit [1,%d] FAIL %.0f s: %s' % (T, time.time() - t0, str(e)[:80])); break

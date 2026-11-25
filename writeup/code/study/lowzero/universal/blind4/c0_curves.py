"""semistable elliptic-curve isogeny classes with N <= 250 from the box a1, a3 in {0,1}, a2 in {-1,0,1}, |a4| <= 60, |a6| <= 150
(the minimal model of the first box curve represents the class; class key = (N, a_p for p <= 97) on MINIMAL models); excluded: every class used in round 1
(11a, 14a, 15a, 17a, 19a, 37a) or round 151 (Sym^2 and Sym^3 tests and proxies, the probed 30a). Writes curves.json."""
import os, json, time
import cypari2
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
PR = [int(p) for p in pari.primes([2, 97])]
t0 = time.time(); classes = {}
for a1 in (0, 1):
    for a2 in (-1, 0, 1):
        for a3 in (0, 1):
            for a4 in range(-60, 61):
                for a6 in range(-150, 151):
                    try:
                        E = pari.ellinit([a1, a2, a3, a4, a6])
                        if len(E) == 0: continue
                        N = int(pari.ellglobalred(E)[0])
                    except Exception:
                        continue
                    if N > 250 or not pari.issquarefree(N): continue
                    E = pari.ellminimalmodel(E)                      # keys from MINIMAL models (ellap needs minimality at p)
                    key = (N,) + tuple(int(pari.ellap(E, p)) for p in PR)
                    if key not in classes:
                        classes[key] = dict(N=N, ainv=[int(x) for x in E[:5]], w=int(pari.ellrootno(E)))
print('classes', len(classes), 'secs', round(time.time() - t0, 1), flush=True)
R151 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'blind', 'testset.json')))
used = [[0, -1, 1, -10, -20], [1, 0, 1, 4, -6], [1, 1, 1, -10, -10], [1, -1, 1, -1, -14], [0, 1, 1, -9, -15], [0, 0, 1, -1, 0], [1, 0, 1, -19, 26]]
for cls in ('T3', 'T4a'):
    used += [x[1] for x in R151[cls]['test'] + R151[cls]['proxy']]
ukeys = set()
for a in used:
    E = pari.ellminimalmodel(pari.ellinit(a)); N = int(pari.ellglobalred(E)[0]); ukeys.add((N,) + tuple(int(pari.ellap(E, p)) for p in PR))
free = sorted([v for k, v in classes.items() if k not in ukeys], key=lambda v: (v['N'], v['ainv']))
print('free classes', len(free), ' w=-1:', sum(1 for v in free if v['w'] == -1), ' w=+1:', sum(1 for v in free if v['w'] == 1))
for lim in (100, 150, 200, 250):
    print('  N <= %d: w=-1 %d, w=+1 %d' % (lim, sum(1 for v in free if v['w'] == -1 and v['N'] <= lim), sum(1 for v in free if v['w'] == 1 and v['N'] <= lim)))
json.dump(dict(free=free, n_used_excluded=len(ukeys)), open('curves.json', 'w'), indent=1)

# Reviewer's own curve data: PARI ellinit from the a-invariants, conductor, root number, analytic rank,
# a_p for p <= P (all primes, good and bad), isogeny check for equal conductors.
import os, json, sys, time
import numpy as np
import cypari2

P = int(sys.argv[1]) if len(sys.argv) > 1 else 500000
pari = cypari2.Pari()
pari.allocatemem(300 * 10**6)
d = json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))
cur = d['curves']
primes = np.array([int(x) for x in pari("primes(primepi(%d))" % P)], dtype=np.int64)
ids, Ns, ws, ranks, semis = [], [], [], [], []
AP = np.zeros((len(cur), len(primes)), dtype=np.int16)
t0 = time.time()
for k, c in enumerate(cur):
    a = c['ainv']
    E = pari.ellinit(a)
    gr = pari.ellglobalred(E)
    N = int(gr[0])
    w = int(pari.ellrootno(E))
    rk = int(pari.ellanalyticrank(E)[0])
    # semistable: N squarefree
    semis.append(bool(pari.issquarefree(N)))
    pari('global(EE); EE = ellinit(%s)' % str(a))
    v = pari('my(v=vectorsmall(%d),i=0);forprime(p=2,%d,i++;v[i]=ellap(EE,p));v' % (len(primes), P))
    AP[k] = np.array([int(x) for x in v], dtype=np.int16)
    ids.append(c['id']); Ns.append(N); ws.append(w); ranks.append(rk)
    assert N == c['N'], (c['id'], N)
    assert w == c['w'], (c['id'], w)
print('curves done', time.time() - t0, 's; all semistable:', all(semis))
ids = np.array(ids); Ns = np.array(Ns); ws = np.array(ws); ranks = np.array(ranks)
print('rank counts', {r: int((ranks == r).sum()) for r in set(ranks.tolist())})
print('w vs rank parity ok:', bool(np.all(ws == (-1) ** ranks)))
# equal-conductor pairs: isogenous iff a_p agree for all p (check first 2000 primes)
for i in range(len(ids)):
    for j in range(i + 1, len(ids)):
        if Ns[i] == Ns[j]:
            same = int((AP[i, :2000] == AP[j, :2000]).sum())
            corr = float(np.mean((AP[i] / np.sqrt(primes)) * (AP[j] / np.sqrt(primes))))
            print('equal N', Ns[i], ids[i], ids[j], 'same a_p in first 2000 primes:', same, ' mean lam1*lam2 (p<=P):', round(corr, 4))
# any two curves (any N) with identical a_p on the first 200 primes (excluding bad)
nid = 0
for i in range(len(ids)):
    for j in range(i + 1, len(ids)):
        if np.array_equal(AP[i, :200], AP[j, :200]):
            nid += 1; print('identical first 200 a_p:', ids[i], ids[j])
print('pairs with identical a_p on first 200 primes:', nid)
np.savez_compressed(os.path.dirname(os.path.abspath(__file__))+'/rv_curves.npz',
                    ids=ids, N=Ns, w=ws, rank=ranks, primes=primes, ap=AP)
print('saved', time.time() - t0)

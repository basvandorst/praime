"""H1: own recomputation of c = prod_p G_p(1/p) for genus-2 L-functions.
Exact local factors: p = 2 and bad p from lfuneuler (PARI object passes lfuncheckfeq); odd good p <= P0 from
hyperellcharpoly. Good p in (P0, Pmax]: a_p by own point count on (2y + h)^2 = F = 4f + h^2 (numpy), first-order factor
1 - (e1^2 - 1)/p, e1 = a_p/sqrt(p). Writes cown/<key>.npz: primes, log G_p, e1 (all p <= Pmax), e2 (p <= P0 exact, else nan).
usage: c_own.py Pmax key1 key2 ...  (keys from g2data.load())"""
import sys, os, math, time, json
import numpy as np, cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g2data import load
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6); pari.default('realprecision', 38)
Pmax = int(float(sys.argv[1])); keys = sys.argv[2:]; P0 = 1000
OUT = 'cown'; os.makedirs(OUT, exist_ok=True)
rows = {r['key']: r for r in load()}
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
PR = primes_upto(Pmax)
cur = []
for key in keys:
    r = rows[key]; f, h = r['spec'][1], r['spec'][2]
    F = [int(x) for x in pari.Vecrev(pari('4*(%s)+(%s)^2' % (f, h)))]
    L = pari('lfungenus2([%s,%s])' % (f, h)); N = int(pari.lfunparams(L)[0])
    cur.append(dict(key=key, f=f, h=h, F=F, L=L, N=N, logG=np.zeros(len(PR)), e1=np.zeros(len(PR)), e2=np.full(len(PR), np.nan), sgn=np.ones(len(PR)),
                    bad=[int(p) for p in pari.factor(N)[0]]))
def Gfrom(P):
    """P = analytic coefficients [1, c1, .., cd] of 1/L_p(X); returns F_p(1/p)/(1 - 1/p) (G_p(1/p) at good p), e1, e2."""
    d = len(P) - 1; J = 12
    lam = [1.0] + [0.0] * J            # coefficients of 1/P
    for j in range(1, J + 1):
        lam[j] = -sum(P[i] * lam[j - i] for i in range(1, min(j, d) + 1))
    return None, lam
t0 = time.time()
for i, p in enumerate(PR):
    p = int(p)
    for C in cur:
        if p == 2 or C['N'] % p == 0 or p <= P0:
            if p == 2 or C['N'] % p == 0:
                e = pari.lfuneuler(C['L'], p); Pp = 1 / e
                co = [float(c) for c in pari.Vecrev(Pp)] if str(pari.type(Pp)) == 't_POL' else [float(Pp)]
            else:
                cp = pari.hyperellcharpoly(pari.Mod(1, p) * pari.Pol(pari.Vec(C['F'][::-1])))
                co = [float(c) for c in pari.Vec(cp)]          # x^4 + c3 x^3 + ... : L-poly coefficients 1, c3, c2, c1, c0
            co = [c / co[0] for c in co]
            Pa = [c * p ** (-j / 2) for j, c in enumerate(co)]   # analytic normalisation
            _, lam = Gfrom(Pa)
            Fv = sum(lam[j] * Pa[j] * p ** (-j) for j in range(len(Pa)))   # real coefficients: conj = identity
            C['logG'][i] = math.log(abs(Fv / (1 - 1 / p))); C['sgn'][i] = 1 if Fv > 0 else -1; C['e1'][i] = -Pa[1] if len(Pa) > 1 else 0.0
            C['e2'][i] = Pa[2] if len(Pa) > 2 else 0.0
        else:
            Fc = C['F'] + [0] * (7 - len(C['F']))
            if cur[0].get('qr_p') != p:
                x = np.arange(p, dtype=np.int64); qr = np.full(p, -1, dtype=np.int64); qr[(x * x) % p] = 1; qr[0] = 0
                hx = np.arange(1, (p - 1) // 2 + 1, dtype=np.int64); y = (hx * hx) % p
                cur[0]['qr'] = qr; cur[0]['qr_p'] = p; cur[0]['hx'] = hx; cur[0]['y'] = y
            qr, hx, y = cur[0]['qr'], cur[0]['hx'], cur[0]['y']
            # F(x) = E(x^2) + x O(x^2); x and p - x together; x = 0 separately
            E = np.full(len(y), Fc[6] % p, dtype=np.int64)
            for a in (Fc[4], Fc[2], Fc[0]): E = (E * y + a) % p
            O = np.full(len(y), Fc[5] % p, dtype=np.int64)
            for a in (Fc[3], Fc[1]): O = (O * y + a) % p
            O = (O * hx) % p
            S = int(qr[(E + O) % p].sum() + qr[(E - O) % p].sum() + qr[Fc[0] % p])
            Fc = C['F']
            lead = Fc[-1] % p if len(Fc) == 7 else 0
            ap = -S - (int(qr[lead]) if len(Fc) == 7 else 0)
            if len(Fc) == 7 and lead == 0: raise RuntimeError('degree drops mod %d' % p)
            e1 = ap / math.sqrt(p); C['e1'][i] = e1; C['logG'][i] = math.log(1 - (e1 * e1 - 1) / p)
for C in cur:
    # own check of the point counts against lfunan (p <= 3000, good p > P0 only where own counts were used: none <= P0;
    # so check the count routine separately on 1000 < p <= 3000)
    an = pari.lfunan(C['L'], min(Pmax, 20000))
    sel = (PR > P0) & (PR <= min(Pmax, 20000)) & np.array([C['N'] % int(p) != 0 for p in PR])
    mism = int(np.sum(np.abs(C['e1'][sel] * np.sqrt(PR[sel]) - np.array([float(an[int(p) - 1]) for p in PR[sel]])) > 1e-6))
    selx = (PR <= P0)
    mism2 = int(np.sum(np.abs(C['e1'][selx] * np.sqrt(PR[selx]) - np.array([float(an[int(p) - 1]) for p in PR[selx]])) > 1e-6))
    np.savez_compressed(os.path.join(OUT, C['key'] + '.npz'), p=PR, logG=C['logG'], e1=C['e1'], e2=C['e2'], sgn=C['sgn'], N=C['N'], bad=C['bad'])
    print('%-14s N %4d bad %s  c(P0) %.5f c(Pmax) %.5f  a_p mismatches vs lfunan: own-count %d, exact-part %d  old c %.5f' % (
        C['key'], C['N'], C['bad'], np.prod(C['sgn']) * math.exp(C['logG'][PR <= P0].sum()), np.prod(C['sgn']) * math.exp(C['logG'].sum()), mism, mism2, rows[C['key']]['c']), flush=True)
print('secs %.1f' % (time.time() - t0))

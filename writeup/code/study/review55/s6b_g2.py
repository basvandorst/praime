"""Genus-2 constants: exact G_p(1/p) for p <= PEX (hyperellcharpoly; lfuneuler at bad p and p = 2), first order beyond
(a_p from own point count of y^2 = 4f + h^2 over F_p, numpy), second-order mean correction from the exact range,
tail nu/(sqrt P ln P) with nu = 2 (Ad = std5 + ad10, both orthogonal)."""
import json, math, os, sys, time
import numpy as np, cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s6_const import F_of
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'; D = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(B + '/predictions_blind.json'))
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
PEX = 1000; P1 = int(float(sys.argv[1])); keys = sys.argv[2:]

def count_ap(Fco, p):
    """a_p = p + 1 - #C(F_p) for y^2 = F(x), deg F = 5 or 6, p odd good."""
    x = np.arange(p, dtype=np.int64); v = np.zeros(p, dtype=np.int64)
    for c in Fco:                       # highest first
        v = (v * x + (c % p)) % p
    sq = np.zeros(p, bool); sq[(x * x) % p] = True
    chi = np.where(v == 0, 0, np.where(sq[v], 1, -1))
    aff = p + int(chi.sum())
    deg = len(Fco) - 1
    while deg >= 0 and Fco[len(Fco) - 1 - deg] % p == 0: deg -= 1
    lead = Fco[len(Fco) - 1 - deg] % p
    if deg == 6: inf = 1 + (1 if sq[lead] else -1)
    elif deg == 5: inf = 1
    else: inf = None
    return p + 1 - (aff + inf)

res = {}
for key in keys:
    e = P['L'][key]; t0 = time.time()
    f = pari(e['spec'][1]); h = pari(e['spec'][2]); Fp = 4 * f + h ** 2
    Fco = [int(c) for c in pari.Vec(Fp)]
    L = pari('lfungenus2([%s,%s])' % (e['spec'][1], e['spec'][2])); N = int(e['Q'])
    ps = np.array([int(x) for x in pari('primes([2,%d])' % P1)], dtype=np.int64)
    logG = np.zeros(len(ps)); chk = []; sgn = 1
    for i, p in enumerate(ps):
        p = int(p)
        if p <= PEX:
            if N % p == 0 or p == 2:
                eu = pari.lfuneuler(L, p); co = np.array([complex(x) for x in pari.Vec(pari.denominator(eu))][::-1]); co = co / co[0]
                coA = co * (p ** -0.5) ** np.arange(len(co))
                al = 1 / np.roots(coA[::-1]) if len(coA) > 1 else np.array([])
            else:
                cp = pari.hyperellcharpoly(Fp * pari.Mod(1, p)); v = [float(x) for x in pari.Vec(cp)]
                al = np.roots(v) / math.sqrt(p)
                if i % 40 == 5: chk.append((p, int(round(-v[1])), count_ap(Fco, p)))
            a = np.zeros((1, 4), complex); a[0, :len(al)] = al
            G = F_of(a, np.array([1.0 / p])) / (1 - 1.0 / p)
            logG[i] = math.log(abs(G.real[0])); sgn *= (1 if G.real[0] > 0 else -1)
        else:
            if N % p == 0:
                eu = pari.lfuneuler(L, p); co = np.array([complex(x) for x in pari.Vec(pari.denominator(eu))][::-1]); co = co / co[0]
                coA = co * (p ** -0.5) ** np.arange(len(co)); al = 1 / np.roots(coA[::-1]) if len(coA) > 1 else np.array([])
                a = np.zeros((1, 4), complex); a[0, :len(al)] = al
                gg = (F_of(a, np.array([1.0 / p])) / (1 - 1.0 / p)).real[0]; logG[i] = math.log(abs(gg)); sgn *= (1 if gg > 0 else -1); continue
            ap = count_ap(Fco, p); lam2 = ap * ap / p
            logG[i] = math.log(1 - (lam2 - 1) / p)
    # second-order correction: mean of p^2 (logG_exact - log(1 - (lam^2-1)/p)) on good p in (300, PEX]
    sel = []
    for i, p in enumerate(ps):
        p = int(p)
        if 300 < p <= PEX and N % p:
            sel.append(i)
    # recompute first-order there for the difference
    diffs = []
    for i in sel:
        p = int(ps[i]); ap = count_ap(Fco, p); diffs.append(p * p * (logG[i] - math.log(1 - (ap * ap / p - 1) / p)))
    c2 = float(np.mean(diffs)); pr = ps[ps > PEX].astype(float); tail2 = float(np.sum(c2 / pr ** 2)) + c2 / (P1 * math.log(P1))
    cum = np.cumsum(logG); nu = 2
    out = {}
    for Pc in (1e4, 3e4, 1e5, 2e5):
        if Pc > P1 * 1.01: continue
        i = np.searchsorted(ps, Pc, side='right') - 1
        sel2 = (ps >= Pc / 10) & (ps <= Pc)
        lc = cum[sel2] + nu / (np.sqrt(ps[sel2]) * np.log(ps[sel2]))
        ces = math.exp(np.sum(lc / ps[sel2]) / np.sum(1 / ps[sel2]))
        out[int(Pc)] = (sgn * math.exp(cum[i]), sgn * math.exp(cum[i] + nu / (math.sqrt(Pc) * math.log(Pc))), sgn * ces)
    Pm = max(out)
    best = out[Pm][2] * math.exp(tail2)
    print('%-12s reg c = %8.5f | ' % (key, e['c']) + '  '.join('P=%.0e raw %.5f corr %.5f ces %.5f' % (k, *v) for k, v in sorted(out.items()))
          + ' | 2nd-order corr %+.5f -> best %.5f, best/reg %.4f | ap check %s %.0fs' % (tail2, best, best / e['c'], chk[:3], time.time() - t0), flush=True)
    res[key] = dict(out={str(k): v for k, v in out.items()}, tail2=tail2, best=best, reg=e['c'])
json.dump(res, open(D + '/s6b_g2_%s.json' % ('_'.join(k.split('_')[1] for k in keys)), 'w'))

"""DRAFT (not run as a test). Blind test of the genus-2 deficit: candidate selection, resumable.
Boxes: B' = {y^2 + h y = f: f in {-2..2}^7, deg f in {5, 6}, at least one coefficient of absolute value 2} and the old box
B = {f in {-1,0,1}^7} (rounds 151/152, review 58) above the old conductor bound; h one of the 16 0/1-polynomials of degree <= 3.
(Dry run, B' alone with 1000 <= N_L <= 3000: 4 accepted L-functions in the first 22.7 % of B', so the conductor range is
widened to 8000; the expected se is fixed at the freeze from the drawn constants.)
Rules, in order (all before any zero is computed): squarefree 4f + h^2 of degree >= 5; genus2red N_red < 10000; one per
fingerprint (N_red, hyperellcharpoly at the six smallest good odd primes); lfungenus2: 1000 <= N_L <= 8000, lfuncheckfeq <= -20,
one per (N_L, a_1..a_60); every L-function used before excluded (rounds 1, 151, 152, review-58 replication, proxies: by the
same key); primitive: mean |lambda_p|^2 over good p <= 3000 in [0.75, 1.25]; generic Sato-Tate: |z(M4)| <= 4 with
M4 = mean (a_p/sqrt p)^4 over good 100 < p <= 3e4, z = (M4 - 3)/sqrt(75/n); root number eps = lfunrootres; constant c by the
rule of H1 (own code: exact local factors p <= 1000, first order to 1e5 from own point counts, log-Cesaro over [1e4, 1e5],
nu = 2). Writes candidates_g2b.jsonl (one line per accepted L-function, all strata).
usage: (ulimit -v 1500000; timeout 600 nice .venv/bin/python -I g2b_select.py BUDGET_SECONDS)"""
import sys, os, json, time, itertools, math
import numpy as np, cypari2
HERE = os.path.dirname(os.path.abspath(__file__))
B4 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind4/'; B1 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind/'
RV = '/home/bas/praime/writeup/code/study/review58/rep/'
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
budget = float(sys.argv[1]); t0 = time.time()
OUT = os.path.join(HERE, 'candidates_g2b.jsonl'); STATE = os.path.join(HERE, 'select_state.json')
hs = ['+'.join(['x^%d' % j for j in range(4) if bits >> j & 1]) or '0' for bits in range(16)]
def old_keys():
    old = []
    for k, e in json.load(open(B4 + 'predictions_blind4.json'))['L'].items():
        if e['cls'] == 'g2': old.append(e['spec'][1:3])
    for k, e in json.load(open(B1 + 'predictions_blind.json'))['L'].items():
        if e['cls'] == 'T4b': old.append(e['spec'][1:3])
    T1 = json.load(open(B1 + 'testset.json'))['T4b']; old += [[r[1], r[2]] for r in T1['test'] + T1['proxy'] + T1['all_valid']]
    old += [[r[1], r[2]] for r in json.load(open(B4 + 'genus2.json'))['valid']]
    old += [['x^2+x', 'x^3+1'], ['x^5-2*x^3+2*x^2-x', '1'], ['x^5+x^4', 'x^3+x+1']]
    keys = set()
    for f, h in old:
        L = pari('lfungenus2([%s,%s])' % (f, h)); an = pari.lfunan(L, 60); keys.add((int(pari.lfunparams(L)[0]),) + tuple(int(x) for x in an))
    return keys
st = json.load(open(STATE)) if os.path.exists(STATE) else dict(i=0, fps=[], keys=None)
if st['keys'] is None: st['keys'] = [list(k) for k in old_keys()]
seen = set(tuple(k) for k in st['keys']); fps = set(tuple(x) for x in st['fps'])
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
PR = primes_upto(100000)
def c_rule(L, Fc, N):
    """H1 rule (own code): exact p <= 1000 (lfuneuler at 2 and bad p, hyperellcharpoly otherwise), own counts to 1e5."""
    lg = np.zeros(len(PR)); sg = 1.0; e1s = np.zeros(len(PR))
    Fpol = pari.Pol(pari.Vec(Fc[::-1])); Fz = Fc + [0] * (7 - len(Fc))
    for i, p in enumerate(PR):
        p = int(p)
        if p <= 1000 or N % p == 0:
            if p == 2 or N % p == 0:
                Pp = 1 / pari.lfuneuler(L, p); co = [float(c) for c in pari.Vecrev(Pp)] if str(pari.type(Pp)) == 't_POL' else [float(Pp)]
            else:
                co = [float(c) for c in pari.Vec(pari.hyperellcharpoly(pari.Mod(1, p) * Fpol))]
            co = [c / co[0] * p ** (-j / 2) for j, c in enumerate(co)]
            lam = [1.0] + [0.0] * 12
            for j in range(1, 13): lam[j] = -sum(co[k] * lam[j - k] for k in range(1, min(j, len(co) - 1) + 1))
            Fv = sum(lam[j] * co[j] * p ** (-j) for j in range(len(co))); lg[i] = math.log(abs(Fv / (1 - 1 / p))); sg *= (1 if Fv > 0 else -1)
            e1s[i] = -co[1] if len(co) > 1 else 0.0
        else:
            hx = np.arange(1, (p - 1) // 2 + 1, dtype=np.int64); y = (hx * hx) % p
            qr = np.full(p, -1, dtype=np.int64); qr[(np.arange(p, dtype=np.int64) ** 2) % p] = 1; qr[0] = 0
            E = np.full(len(y), Fz[6] % p, dtype=np.int64)
            for a in (Fz[4], Fz[2], Fz[0]): E = (E * y + a) % p
            O = np.full(len(y), Fz[5] % p, dtype=np.int64)
            for a in (Fz[3], Fz[1]): O = (O * y + a) % p
            O = (O * hx) % p
            S = int(qr[(E + O) % p].sum() + qr[(E - O) % p].sum() + qr[Fz[0] % p])
            ap = -S - (int(qr[Fc[-1] % p]) if len(Fc) == 7 else 0); e1 = ap / math.sqrt(p); e1s[i] = e1; lg[i] = math.log(1 - (e1 * e1 - 1) / p)
    cum = np.cumsum(lg); sel = (PR >= 10000); w = 1.0 / PR[sel]
    c = sg * math.exp(float(np.sum((cum[sel] + 2 / (np.sqrt(PR[sel]) * np.log(PR[sel]))) * w) / w.sum()))
    good = np.array([N % int(q) != 0 for q in PR]) & (PR > 100) & (PR <= 30000)
    M4 = float(np.mean(e1s[good] ** 4)); z4 = (M4 - 3) / math.sqrt(75 / good.sum())
    return c, M4, z4
coefs_all = [c for c in itertools.product((-2, -1, 0, 1, 2), repeat=7) if (c[5] or c[6]) and max(abs(x) for x in c) == 2] + \
            [c for c in itertools.product((-1, 0, 1), repeat=7) if (c[5] or c[6])]
while st['i'] < len(coefs_all) * 16:
    if time.time() - t0 > budget - 60: break
    ci, hi = divmod(st['i'], 16); st['i'] += 1; coefs = coefs_all[ci]; h = hs[hi]
    f = '+'.join('(%d)*x^%d' % (c, j) for j, c in enumerate(coefs) if c != 0)
    try:
        F = pari('4*(%s)+(%s)^2' % (f, h))
        if int(pari.poldegree(F)) < 5 or int(pari.poldegree(pari.gcd(F, pari.deriv(F)))) > 0: continue
        Nr = int(pari('genus2red([%s,%s])' % (f, h))[0])
    except cypari2.PariError:
        continue
    if Nr >= 10000: continue
    fpl = []
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        if Nr % p == 0 or len(fpl) == 6: continue
        try: fpl.append(str(pari('hyperellcharpoly(Mod(1,%d)*(4*(%s)+(%s)^2))' % (p, f, h))))
        except cypari2.PariError: fpl.append('sing%d' % p)          # model singular at p although p does not divide N_red
    fp = (Nr,) + tuple(fpl)
    if fp in fps: continue
    fps.add(fp); st['fps'].append(list(fp))
    try:
        L = pari('lfungenus2([%s,%s])' % (f, h)); NL = int(pari.lfunparams(L)[0])
        if not (1000 <= NL <= 8000): continue
        an = pari.lfunan(L, 3000); key = (NL,) + tuple(int(an[i]) for i in range(60))
        if key in seen: continue
        seen.add(key); st['keys'].append(list(key))
        feq = float(pari.lfuncheckfeq(L))
        if feq > -20: continue
        ps = [int(p) for p in pari.primes([2, 3000]) if NL % int(p)]
        m2 = float(np.mean([(float(an[p - 1]) / math.sqrt(p)) ** 2 for p in ps]))
        if not (0.75 <= m2 <= 1.25): continue
        eps = float(pari.real(pari.lfunrootres(L)[2]))
        Fc = [int(x) for x in pari.Vecrev(F)]
        c, M4, z4 = c_rule(L, Fc, NL)
    except cypari2.PariError as e:
        continue
    row = dict(N=NL, f=f, h=h, feq=feq, M2=m2, M4=M4, zM4=z4, eps=eps, c=c)
    open(OUT, 'a').write(json.dumps(row) + '\n'); print(row, flush=True)
json.dump(st, open(STATE, 'w'))
print('state i = %d of %d (%.1f %%), %.0f s' % (st['i'], len(coefs_all) * 16, 100 * st['i'] / (len(coefs_all) * 16), time.time() - t0))

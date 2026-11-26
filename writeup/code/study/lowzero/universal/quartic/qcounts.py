"""qcounts.py: field lists and per-field bin sums for the S4 quartic murmuration (used by quartic_test.py).
Per field (S4 quartic, PARI nflist, r2 = 0, 1, 2): a_K(p) = #roots of the defining quartic mod p - 1 for p not dividing disc(f);
for p | index (p not dividing D): number of degree-1 primes above p minus 1 (idealprimedec); p | D dropped.
Fine bins in u = y^(1/3), y = p/|D|, width 1/96 on (0, UMAX]: num = sum ln p sqrt|D| a_K(p), den = sum ln p.
Also the splitting type of l = 2, 3, 5, 7 (sorted list of (e, f))."""
import os, time, math, numpy as np
UMAX = 2.0; NF = 192; FEDGES = np.linspace(0.0, UMAX, NF + 1)
CHUNK = 5000
SMALL_L = (2, 3, 5, 7, 11, 13)

def pari_init():
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(200_000_000); pari.default('parisizemax', 400_000_000)
    pari('acnt(f, P) = apply(p -> #polrootsmod(f, p) - 1, P)')
    pari('acnt1(nf, p) = #select(q -> q.f == 1, idealprimedec(nf, p)) - 1')
    pari('ltype(nf, l) = vecsort(apply(q -> [q.e, q.f], idealprimedec(nf, l)))')
    return pari

def type_name(ef):
    """[(e,f),...] -> symbol like '1^2 11' (as in pred4.TYPES)."""
    ef = sorted(ef, key=lambda x: (-x[0], x[1]))
    ram = [(e, f) for e, f in ef if e > 1]; un = sorted([f for e, f in ef if e == 1], reverse=True)
    if not ram: return ''.join(str(f) for f in un) if un != [1, 1, 1, 1] else '1111'
    s = ' '.join('%d^%d' % (f, e) for e, f in ram)
    if un: s += ' ' + ''.join(str(f) for f in sorted(un))
    return s

def field_list(pari, fn, lo, hi, r2, log=print):
    done = {}
    if os.path.exists(fn):
        for line in open(fn):
            if line.startswith('#chunk'): done[int(line.split()[1])] = True
    with open(fn, 'a') as out:
        for c0 in range(lo, hi, CHUNK):
            if c0 in done: continue
            t0 = time.time(); c1 = min(hi, c0 + CHUNK) - 1
            fl = pari('nflist("S4", [%d, %d], %d)' % (c0, c1, r2))
            nin = 0
            for f in fl:   # nflist("S4", [a, b]) also returns fields with |D| < a: keep a <= |D| <= b only
                D = int(pari.nfdisc(f))
                if c0 <= abs(D) <= c1: out.write('%d %s\n' % (D, str(f).replace(' ', ''))); nin += 1
            out.write('#chunk %d %d fields %.1f s\n' % (c0, nin, time.time() - t0)); out.flush()
            log('  nflist S4 r2=%d [%d, %d]: %d fields (%d returned), %.1f s' % (r2, c0, c1, nin, len(fl), time.time() - t0))
    rows = [l.split() for l in open(fn) if not l.startswith('#')]
    rows = sorted(set((int(D), f) for D, f in rows), key=lambda r: (abs(r[0]), r[1]))
    for D, f in rows: assert (D < 0) == (r2 == 1) and lo <= abs(D) < hi
    return rows

def field_sums(pari, fstr, D):
    f = pari(fstr); A = abs(D); ymax = UMAX**3; P = pari.primes([2, int(ymax * A)])
    a = np.array([int(v) for v in pari('acnt')(f, P)], float); p = np.array([int(v) for v in P], np.int64)
    keep = (A % p != 0); ind2 = abs(int(pari.poldisc(f))) // A    # disc(f) = index^2 D
    nf = pari.nfinit(f)
    idx = [i for i in np.nonzero(keep)[0] if ind2 % int(p[i]) == 0]
    for i in idx:
        a[i] = int(pari('acnt1')(nf, int(p[i])))
    p = p[keep]; a = a[keep]; u = np.cbrt(p / A)
    j = np.searchsorted(FEDGES, u, side='left') - 1
    ok = (j >= 0) & (j < NF); lp = np.log(p[ok].astype(float))
    num = np.bincount(j[ok], weights=lp * math.sqrt(A) * a[ok], minlength=NF)
    den = np.bincount(j[ok], weights=lp, minlength=NF)
    types = [type_name([(int(e), int(ff)) for e, ff in pari('ltype')(nf, l)]) for l in SMALL_L]
    return num, den, len(p), types

def counts(stem, lo, hi, r2, budget, log=print):
    pari = pari_init(); t0 = time.time()
    rows = field_list(pari, stem + '_fields.txt', lo, hi, r2, log)
    st = stem + '.npz'
    if os.path.exists(st):
        z = np.load(st, allow_pickle=False); num, den, D, npr, done, typ = z['num'], z['den'], z['D'], z['npr'], z['done'], z['types']
        assert len(D) == len(rows) and all(int(D[i]) == rows[i][0] for i in range(len(rows)))
        typ = typ.astype('<U12')
    else:
        num = np.zeros((len(rows), NF)); den = np.zeros((len(rows), NF)); D = np.array([r[0] for r in rows], np.int64)
        npr = np.zeros(len(rows), np.int64); done = np.zeros(len(rows), bool); typ = np.full((len(rows), len(SMALL_L)), '', dtype='<U12')
    todo = np.nonzero(~done)[0]; nrun = 0; tr = time.time(); nroot = 0
    for i in todo:
        if time.time() - t0 > budget: break
        num[i], den[i], npr[i], ty = field_sums(pari, rows[i][1], rows[i][0]); typ[i] = ty; done[i] = True; nrun += 1; nroot += npr[i]
        if nrun % 200 == 0: np.savez(st, num=num, den=den, D=D, npr=npr, done=done, types=typ)
    np.savez(st, num=num, den=den, D=D, npr=npr, done=done, types=typ)
    dt = time.time() - tr
    log('%s: %d fields, %d done (%d this run, %.1f s, %.2f us per prime)%s' % (os.path.basename(stem), len(rows), int(done.sum()), nrun, dt,
        1e6 * dt / max(1, nroot), '  DONE' if done.all() else '  (run again to continue)'))
    return bool(done.all())

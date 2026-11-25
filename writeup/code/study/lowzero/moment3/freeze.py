"""freeze.py: the registered sets and the frozen predictions (run BEFORE any |L| value of a reserved character exists).

  freeze.py draw                 -> sets.json  (dev set, reserved set R: 56 fundamental discriminants, seeded draw;
                                                the 1540 pairs of R; a seeded sample of 2000 of the 27720 triples of R)
  freeze.py singles|pairs|triples [--budget 540]   -> parts/*.jsonl (resumable chunks, each run < 10 minutes)
  freeze.py assemble             -> predictions_triple.json + its sha256 (printed and written to predictions_triple.sha256)
  freeze.py dev                  -> predictions_dev.json (the same quantities for the development set, for the dry run)

Draw: all fundamental discriminants D with 100 <= |D| <= 1300 (both signs), minus the 8 development characters;
numpy default_rng(20261125).permutation of that list; walk it and keep D unless |D| is already kept (distinct conductors),
until 56 are kept.  Triples: default_rng(20261126).choice(27720, 2000, replace=False) on the lexicographic list of the
triples of the kept set sorted by D.  Predictions: tmlib (pmax = 1e5; 1e5 -> 1e6 changes P by < 1e-8 relative).
"""
import sys, os, json, time, math, itertools, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import tmlib as T

DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
NRES, NTRI = 56, 2000
T1, T2 = 100.0, 1000.0
PM = 10**5
PARTS = os.path.join(HERE, 'parts')


def draw():
    cand = [D for D in range(-1300, 1301) if 100 <= abs(D) <= 1300 and T.is_fundamental(D) and D not in DEV]
    perm = np.random.default_rng(20261125).permutation(len(cand))
    kept, absd = [], set()
    for i in perm:
        D = cand[i]
        if abs(D) in absd:
            continue
        kept.append(D); absd.add(abs(D))
        if len(kept) == NRES:
            break
    R = sorted(kept)
    pairs = list(itertools.combinations(R, 2))
    alltri = list(itertools.combinations(R, 3))
    pick = np.sort(np.random.default_rng(20261126).choice(len(alltri), NTRI, replace=False))
    tri = [alltri[i] for i in pick]
    sets = {'dev': DEV, 'reserved': R, 'n_candidates': len(cand), 'pairs': pairs, 'triples': tri,
            'rule': __doc__.split('Draw:')[1].strip()}
    json.dump(sets, open(os.path.join(HERE, 'sets.json'), 'w'), indent=0)
    print('candidates %d; reserved %d: %s' % (len(cand), len(R), R))
    print('pairs %d; triples %d of %d' % (len(pairs), len(tri), len(alltri)))
    cnt = {D: sum(D in t for t in tri) for D in R}
    print('triples per character: min %d max %d' % (min(cnt.values()), max(cnt.values())))


def single_pred(D):
    tr = T.Triple((D,), pmax=PM)
    c4 = T.P4_coeffs(D)
    f4 = lambda x: sum(c4[k] * np.log(x) ** k for k in range(5))
    zl = T.Triple((1, D), pmax=PM)
    return {'D': D, 'P1': T.integrate(tr.P, T1, T2), 'P4': T.integrate(f4, T1, T2), 'c4': c4,
            'zetaL': {'full': T.integrate(zl.P, T1, T2), 'diag': T.integrate(zl.diag, T1, T2), 'lead': T.integrate(zl.lead, T1, T2),
                      'hm': list(map(float, zl.hm))},
            'hm1': list(map(float, tr.hm))}


def pair_pred(Dp):
    tr = T.Triple(Dp, pmax=PM)
    return {'Ds': list(Dp), 'full': T.integrate(tr.P, T1, T2), 'diag': T.integrate(tr.diag, T1, T2), 'lead': T.integrate(tr.lead, T1, T2),
            'hm': list(map(float, tr.hm)), 'g': list(map(float, tr.g))}


def p2config(Ds):
    v = sorted(int(T.pari.kronecker(D, 2)) for D in Ds)
    return ''.join({-1: 'm', 0: 'r', 1: 'p'}[x] for x in v)


def triple_pred(Dt, singles, pairs):
    tr = T.Triple(Dt, pmax=PM)
    ml = tr.ml; z = [ml.var(j) for j in range(3)]
    hn = T.h_jet_no3(ml, Dt, z, PM)
    Ls = lambda t: tr.logs(t)

    def P_no3(t):
        L = Ls(t); out = 0.0
        for i, m in enumerate(ml.monos):
            term = hn[i]
            for j in range(3):
                if m[j] == 0:
                    term = term * (L[j] + 2 * T.EULER_GAMMA)
            out = out + term
        return out
    sP = [singles[D] for D in Dt]
    pP = [pairs[(Dt[i], Dt[j])] for i in range(3) for j in range(i + 1, 3)]
    P_ind = lambda t: np.prod([s.P(t) for s in sP], axis=0)
    P_kirk = lambda t: np.prod([p.P(t) for p in pP], axis=0) / np.prod([s.P(t) for s in sP], axis=0)
    out = {'Ds': list(Dt)}
    for name, f in [('full', tr.P), ('diag', tr.diag), ('lead', tr.lead), ('no3', P_no3), ('kirk', P_kirk), ('ind', P_ind)]:
        out[name] = T.integrate(f, T1, T2)
    out['hm'] = list(map(float, tr.hm)); out['hm_no3'] = list(map(float, hn)); out['g'] = list(map(float, tr.g))
    out['A3_over_Aij'] = float(tr.hm[0] / hn[0]); out['p2'] = p2config(Dt)
    tm = np.array([100.0, 316.2277660168, 1000.0])
    out['P_at'] = {'t': tm.tolist(), 'full': list(map(float, tr.P(tm))), 'diag': list(map(float, tr.diag(tm)))}
    return out


def run_stage(stage, budget):
    os.makedirs(PARTS, exist_ok=True)
    sets = json.load(open(os.path.join(HERE, 'sets.json')))
    fn = os.path.join(PARTS, stage + '.jsonl')
    done = set()
    if os.path.exists(fn):
        for line in open(fn):
            done.add(tuple(json.loads(line)['key']))
    if stage == 'singles':
        items = [(D,) for D in sets['reserved']]
    elif stage == 'pairs':
        items = [tuple(p) for p in sets['pairs']]
    else:
        items = [tuple(t) for t in sets['triples']]
    start = time.time(); n = 0
    cache_s, cache_p = {}, {}
    with open(fn, 'a') as fo:
        for it in items:
            if it in done:
                continue
            if time.time() - start > budget:
                break
            if stage == 'singles':
                r = single_pred(it[0])
            elif stage == 'pairs':
                r = pair_pred(it)
            else:
                for D in it:
                    if D not in cache_s:
                        cache_s[D] = T.Triple((D,), pmax=PM)
                for i in range(3):
                    for j in range(i + 1, 3):
                        k = (it[i], it[j])
                        if k not in cache_p:
                            cache_p[k] = T.Triple(k, pmax=PM)
                r = triple_pred(it, cache_s, cache_p)
            r['key'] = list(it)
            fo.write(json.dumps(r) + '\n'); fo.flush(); n += 1
    left = len(items) - len(done) - n
    print('%s: %d computed this run, %d left (%.0f s)' % (stage, n, left, time.time() - start))


def assemble():
    sets = json.load(open(os.path.join(HERE, 'sets.json')))
    out = {'sets': sets, 'T1': T1, 'T2': T2, 'pmax': PM, 'method': 'tmlib.py (recipe = sum over the 8 swap terms, multilinear jets; '
           'diag = AFE diagonal with sharp cutoff at sqrt(C); kirk = prod pairs / prod singles; ind = prod singles; '
           'no3 = A_3 replaced by prod A_ij; lead = h(0) prod ln(q_j t / 2 pi)); integrals by 96-point Gauss-Legendre in ln t'}
    for stage in ['singles', 'pairs', 'triples']:
        rows = [json.loads(l) for l in open(os.path.join(PARTS, stage + '.jsonl'))]
        want = {'singles': len(sets['reserved']), 'pairs': len(sets['pairs']), 'triples': len(sets['triples'])}[stage]
        keys = {tuple(r['key']) for r in rows}
        assert len(keys) == want, (stage, len(keys), want)
        out[stage] = sorted(rows, key=lambda r: r['key'])
    fn = os.path.join(HERE, 'predictions_triple.json')
    json.dump(out, open(fn, 'w'), indent=0)
    h = hashlib.sha256(open(fn, 'rb').read()).hexdigest()
    open(os.path.join(HERE, 'predictions_triple.sha256'), 'w').write(h + '  predictions_triple.json\n')
    print('predictions_triple.json sha256', h)


def dev():
    """the same quantities for the development set (dry run of moment_test.py)."""
    D8 = DEV
    singles = [single_pred(D) for D in D8]
    pairs = [pair_pred(p) for p in itertools.combinations(D8, 2)]
    cs = {D: T.Triple((D,), pmax=PM) for D in D8}
    cp = {p: T.Triple(p, pmax=PM) for p in itertools.combinations(D8, 2)}
    triples = [triple_pred(t, cs, cp) for t in itertools.combinations(D8, 3)]
    for r in singles: r['key'] = [r['D']]
    for r in pairs + triples: r['key'] = r['Ds']
    sets = {'dev': D8, 'reserved': D8, 'pairs': [list(p) for p in itertools.combinations(D8, 2)],
            'triples': [list(t) for t in itertools.combinations(D8, 3)]}
    json.dump({'sets': sets, 'T1': T1, 'T2': T2, 'pmax': PM, 'singles': singles, 'pairs': pairs, 'triples': triples},
              open(os.path.join(HERE, 'predictions_dev.json'), 'w'), indent=0)
    print('predictions_dev.json written')


if __name__ == '__main__':
    st = sys.argv[1]
    budget = float(sys.argv[sys.argv.index('--budget') + 1]) if '--budget' in sys.argv else 540.0
    if st == 'draw':
        draw()
    elif st in ('singles', 'pairs', 'triples'):
        run_stage(st, budget)
    elif st == 'assemble':
        assemble()
    elif st == 'dev':
        dev()

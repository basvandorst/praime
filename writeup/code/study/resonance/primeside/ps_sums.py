"""ps_sums.py: the prime sums of the registered test (resumable, bounded).  Primes only, no zeros.

usage: ps_sums.py RUN [budget_seconds]      RUN = dry (x <= 2e7, the seen range; NOT blind) or fresh (the registered range)
Each chunk of (group, [L, U)) is computed once and appended to sums_RUN.jsonl; a rerun skips finished chunks and stops
when the time budget (default 540 s) is spent.  Run repeatedly until it prints ALL DONE.

Statistics per series (ps_lib.seg_stats): sharp sums sum lam(p) ln p e(-a d (p/Q)^(1/d)) for a = 1 (resonance) and the 8
control factors a in ps_lib.CTRL; smooth sums with weight p^(-1/2+1/(2d)) W(p^(1/d)) for the main Hann window and 8
sub-windows (same frequencies); noise sums.  Group A also: the Pythagorean zone sums (ps_lib.zone_sums, m = 0, 4).
Memory: chunks of 2e7 (A, B, C: arrays of a few 1e6 int64), PARI chunks of ~1e6 primes (D), mfcoefs to 2e6 (Delta)."""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ps_lib as L

HERE = os.path.dirname(os.path.abspath(__file__))

def geom(x0, x1, n):
    return [int(round(x0 * (x1 / x0) ** (j / n))) // 2 * 2 for j in range(n + 1)]

# ------------------------------------------------------------------ registered configuration (also frozen in predictions_primeside.json)
RUNS = {
    'dry': {   # seen range; readings on (6e5, 2e7]; DD on (1e4, 3e5]
        'ABC': dict(read=(600000, 20000000), grid=sorted(set([0, 10000, 100000, 1000000, 10000000] + geom(600000, 20000000, 32)))),
        'D': dict(read=(600000, 20000000), grid=sorted(set([0, 10000, 100000, 1000000, 10000000] + geom(600000, 20000000, 32)))),
        'DD': dict(read=(10000, 300000), grid=sorted(set([0] + geom(10000, 300000, 16)))),
    },
    'fresh': {
        'ABC': dict(read=(30000000, 1000000000), grid=geom(30000000, 1000000000, 32)),
        'D': dict(read=(30000000, 300000000), grid=geom(30000000, 300000000, 24)),
        'DD': dict(read=(300000, 2000000), grid=geom(300000, 2000000, 16)),
    },
}
CHUNK = 20000000
ZONE_DIRS = {'axis0': 0.0, 'axis90': math.pi / 2}
for (u, v) in ((3, 4), (5, 12), (8, 15), (20, 21), (7, 24), (12, 35), (9, 40)):
    ZONE_DIRS['%d+%di' % (u, v)] = math.atan2(v, u); ZONE_DIRS['%d+%di' % (v, u)] = math.atan2(u, v)
for j, ang in enumerate((0.35, 0.55, 0.70, 0.90, 1.05, 1.25, math.pi / 4)):
    ZONE_DIRS['ctl%d' % j] = ang
ZONE_MS = (0, 4)
GROUP_KEYS = {'A': ['A4', 'A8', 'A12', 'A16', 'A20'], 'B': ['B6', 'B12'], 'C23': ['C23'], 'C31': ['C31'], 'C47': ['C47a', 'C47b']}

def windows_for(key, read):
    d = L.SERIES[key][1]; return L.windows(read[0] ** (1 / d), read[1] ** (1 / d))

def chunks(grid, step):
    out = []
    for a, b in zip(grid[:-1], grid[1:]):
        x = a
        while x < b:
            y = min(b, x + step); out.append((x, y)); x = y
    return out

def ser(r):
    return {'sharp': [[z.real, z.imag] for z in r['sharp']], 'smooth': [[[z.real, z.imag] for z in row] for row in r['smooth']],
            'noise_sharp': r['noise_sharp'], 'noise_smooth': list(map(float, r['noise_smooth'])), 'count': r['count']}

def special_in(key, Lo, U):
    sp = L.special_primes(key)
    if key[0] == 'C': sp = {2: sp[2]}            # the odd ramified prime comes from the enumeration (lam = 1)
    return {p: l for p, l in sp.items() if Lo < p <= U}

def add_rec(r1, r2):
    """sum of two serialised records (same windows)."""
    a = lambda x, y: [x[0] + y[0], x[1] + y[1]]
    return {'sharp': [a(x, y) for x, y in zip(r1['sharp'], r2['sharp'])],
            'smooth': [[a(x, y) for x, y in zip(rx, ry)] for rx, ry in zip(r1['smooth'], r2['smooth'])],
            'noise_sharp': r1['noise_sharp'] + r2['noise_sharp'], 'noise_smooth': [x + y for x, y in zip(r1['noise_smooth'], r2['noise_smooth'])],
            'count': r1['count'] + r2['count']}

def do_abc(Lo, U, cfg):
    ok = L.sieve_odd(Lo, U); rec = {}; zones = None
    for grp, keys in GROUP_KEYS.items():
        p, lams, ex = L.lam_series(grp, Lo, U, ok, keys)
        d, Q = 2, L.SERIES[keys[0]][2]
        W = windows_for(keys[0], cfg['read'])
        st = L.seg_stats(p, lams, d, Q, W)
        for k in keys:
            rec[k] = ser(st[k])
            sp = special_in(k, Lo, U)          # p = 2 and the ramified prime (not produced by the odd-norm enumerations)
            if sp:
                s2 = L.seg_stats(np.array(list(sp.keys()), np.int64), {k: np.array(list(sp.values()))}, d, Q, W)[k]
                rec[k] = add_rec(rec[k], ser(s2))
        if grp == 'A':
            zs = L.zone_sums(p.astype(np.float64), ex['a'], ex['b'], list(ZONE_DIRS.values()), ZONE_MS)
            zones = {name: dict({str(m): [zs[(ph, m)].real, zs[(ph, m)].imag] for m in ZONE_MS}, noise=zs[(ph, 'noise')]) for name, ph in ZONE_DIRS.items()}
    return rec, zones

_pari = [None]
def pari():
    if _pari[0] is None:
        import cypari2
        _pari[0] = cypari2.Pari(); _pari[0].allocatemem(900 * 10 ** 6)
    return _pari[0]

def do_d(group, Lo, U, cfg):
    ai = {'D11a': '[0,-1,1,-10,-20]', 'D37a': '[0,0,1,-1,0]'}[group]
    P = pari()
    v = P('my(E=ellinit(%s), Pr=primes([%d,%d])); [Pr, apply(p->ellap(E,p), Pr)]' % (ai, Lo + 1, U))
    p = np.array([int(x) for x in v[0]], np.int64); ap = np.array([int(x) for x in v[1]], np.float64)
    lam = ap / np.sqrt(p); rec = {}
    rec[group] = ser(L.seg_stats(p, {group: lam}, 2, L.SERIES[group][2], windows_for(group, cfg['read']))[group])
    if group == 'D11a':
        lS = lam * lam - 1; lS[p == 11] = 1.0 / 11
        rec['E11a'] = ser(L.seg_stats(p, {'E11a': lS}, 3, 121, windows_for('E11a', cfg['read']))['E11a'])
    return rec

def do_dd(grid, cfg):
    P = pari(); X = grid[-1]
    v = P('my(c=mfcoefs(mfDelta(), %d), Pr=primes([2,%d])); [Pr, apply(p->c[p+1], Pr)]' % (X, X))
    p = np.array([int(x) for x in v[0]], np.int64); tau = np.array([float(int(x)) for x in v[1]])
    lam = tau / p.astype(np.float64) ** 5.5; out = []
    for (a, b) in zip(grid[:-1], grid[1:]):
        m = (p > a) & (p <= b)
        out.append(((a, b), {'DD': ser(L.seg_stats(p[m], {'DD': lam[m]}, 2, 1, windows_for('DD', cfg['read']))['DD'])}))
    return out

def main():
    run = sys.argv[1]; assert run in RUNS
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 540.0
    t0 = time.time(); fn = os.path.join(HERE, 'sums_%s.jsonl' % run)
    done = set()
    if os.path.exists(fn):
        for line in open(fn):
            j = json.loads(line); done.add((j['group'], j['L'], j['U']))
    out = open(fn, 'a'); todo = 0
    # A, B, C
    cfg = RUNS[run]['ABC']
    for (Lo, U) in chunks(cfg['grid'], CHUNK):
        if ('ABC', Lo, U) in done: continue
        if time.time() - t0 > budget: todo += 1; continue
        t1 = time.time(); rec, zones = do_abc(Lo, U, cfg)
        out.write(json.dumps({'group': 'ABC', 'L': Lo, 'U': U, 'res': rec, 'zones': zones}) + '\n'); out.flush()
        print('ABC (%d, %d] %.1fs  total %.0fs' % (Lo, U, time.time() - t1, time.time() - t0), flush=True)
    # elliptic curves (and Sym^2 11a)
    cfg = RUNS[run]['D']
    for group in ('D11a', 'D37a'):
        for (Lo, U) in chunks(cfg['grid'], 2 * 10 ** 7):
            if (group, Lo, U) in done: continue
            if time.time() - t0 > budget: todo += 1; continue
            t1 = time.time(); rec = do_d(group, Lo, U, cfg)
            out.write(json.dumps({'group': group, 'L': Lo, 'U': U, 'res': rec}) + '\n'); out.flush()
            print('%s (%d, %d] %.1fs  total %.0fs' % (group, Lo, U, time.time() - t1, time.time() - t0), flush=True)
    # Delta
    cfg = RUNS[run]['DD']
    if not any(g == 'DD' for (g, a, b) in done):
        if time.time() - t0 > budget: todo += 1
        else:
            for (ab, rec) in do_dd(cfg['grid'], cfg):
                out.write(json.dumps({'group': 'DD', 'L': ab[0], 'U': ab[1], 'res': rec}) + '\n')
            out.flush(); print('DD done total %.0fs' % (time.time() - t0), flush=True)
    print('ALL DONE' if todo == 0 else 'budget spent: %d chunks left, run again' % todo)

if __name__ == '__main__':
    main()

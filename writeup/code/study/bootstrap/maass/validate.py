# Checks every number quoted in maass.md against the data files.  Usage: validate.py   (exit code 1 on any failure)
import json, os, re, sys, math, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
ok = 0; bad = 0
def check(name, got, want, tol):
    global ok, bad
    good = abs(got - want) <= tol
    ok += good; bad += (not good)
    if not good: print('FAIL %s: got %r, quoted %r' % (name, got, want))
def true(name, cond):
    global ok, bad
    ok += bool(cond); bad += (not cond)
    if not cond: print('FAIL %s' % name)
def rows(f): return sorted((json.loads(l) for l in open(os.path.join(D, f))), key=lambda d: (d.get('r', 0), d.get('l1', 0), d.get('l2', 0)))
def at(rs, r): return min(rs, key=lambda d: abs(d['r'] - r))
txt = lambda f: open(os.path.join(D, f)).read()

# --- section 3: characters
t = txt('chars_precise_out.txt')
for kap, X, want in [(1, 8, 2.9999977254), (1, 14, 2.9999999996), (0, 8, 4.9998922240), (0, 14, 4.9999998848), (0, 24, 4.9999999994)]:
    m = re.search(r'kappa %d X %d: smallest conductor ([\d.]+)' % (kap, X), t); check('chars %d %d' % (kap, X), float(m.group(1)), want, 1e-10)
true('chars odd X 24 failed', 'kappa 1 X 24: LP solver failed' in t)
from sympy import kronecker_symbol, primerange
c = txt('chars_out.txt').splitlines()
for kap, Dd, last in [(1, -3, 37), (0, 5, 31)]:
    line = [l for l in c if l.startswith('kappa %d X 40' % kap)][0]
    vals = {int(a): float(b) for a, b in re.findall(r'(\d+):([+-][\d.]+)', line)}
    true('character values D=%d' % Dd, all(abs(vals[p] - int(kronecker_symbol(Dd, p))) < 2e-3 for p in primerange(2, last + 1)))
    if kap == 0: check('chi_5 at 37', vals[37], -0.94, 0.006)

# --- section 4: holomorphic weight
h = rows('hol_X8.jsonl')
for k, want in [(2.0, -2.390), (11.75, -0.0366), (12.0, -0.00012), (12.25, 0.0229), (14.0, 0.165)]: check('hol t*(%g)' % k, at(h, k)['t_lp'], want, 6e-4 if abs(want) > 0.01 else 1e-5)
true('hol positive from 12.25 to 30', all(d['t_lp'] > 0 for d in h if d['r'] >= 12.25) and all(d['t_lp'] < 0 for d in h if d['r'] <= 12.0))
for a, want in zip(at(h, 12.0)['mean_ap'], [-0.5297, 0.5975, 0.6969, -0.4095]): check('tau at support 8', a, want, 6e-5)
tau = {2: -24, 3: 252, 5: 4830, 7: -16744}
for p, want in zip((2, 3, 5, 7), [-0.5303, 0.5987, 0.6912, -0.3765]): check('tau true %d' % p, tau[p]/p**5.5, want, 6e-5)
lv = txt('levels_out.txt')
for X, wants in [(8, [10.906, 4.963, 2.967, 1.951]), (16, [10.936, 4.969, 2.983, 1.958])]:
    for k, want in zip((2, 4, 6, 8), wants):
        m = re.search(r'weight %d X %d N \d+ bad 1: smallest level ([\d.]+)' % (k, X), lv); check('level w%d X%d' % (k, X), float(m.group(1)), want, 6e-4)
m = re.search(r'weight 12 X 16 .*?a_p 2:([+-][\d.]+) 3:([+-][\d.]+) 5:([+-][\d.]+) 7:([+-][\d.]+)', lv)
for got, want in zip(m.groups(), [-0.530, 0.598, 0.691, -0.376]): check('tau at support 16', float(got), want, 6e-4)
th = txt('truth_hol_out.txt').splitlines()
true('truth_hol: 30 lines, all forms pass', len(th) == 30 and all(min(float(x) for x in re.findall(r'(?:f|S2|S3|S4) ([+-][\d.e+-]+)', l)) > -1e-12 for l in th))
true('truth_hol levels', {int(re.match(r'level (\d+)', l).group(1)) for l in th} == {1, 11, 5, 3, 2, 14, 15, 7})

# --- section 5: round 165
key = json.load(open(os.path.join(D, 'lmfdb_key.json'))); true('key: ten forms below 20', sum(1 for f in key if f['r'] < 20) == 10)
tr = txt('truth_out.txt').splitlines(); true('truth: 30 lines', len(tr) == 30)
vals = [[float(x) for x in re.findall(r'(?:f|S2|S3|S4) ([+-][\d.e+-]+)', l)] for l in tr]
true('truth L(f), Sym2 >= -2e-14', min(min(v[0], v[1]) for v in vals) >= -2e-14); true('truth Sym3, Sym4 >= 0.015', min(min(v[2], v[3]) for v in vals) >= 0.015)
S = {(e, X): rows('scan_eps%d_X%d_k4_J61.jsonl' % (e, X)) for e in (0, 1) for X in (8, 16, 32)}
true('grid sizes', [len(S[(1, X)]) for X in (8, 16, 32)] == [381, 381, 191] and [len(S[(0, X)]) for X in (8, 16, 32)] == [381, 381, 191])
for e, r1, wants in [(1, 9.53369526, [-0.128, -0.172, -0.205]), (0, 13.77975135, [-0.211, -0.216, -0.261])]:
    for X, want in zip((8, 16, 32), wants): check('D1 e%d X%d' % (e, X), max(d['t_lp'] for d in S[(e, X)] if d['r'] < r1 - 0.3), want, 6e-4)
o8 = S[(1, 8)]
for r, want in [(9.55, -0.0022), (9.50, -0.013), (9.75, -0.029)]: check('odd X8 t*(%g)' % r, at(o8, r)['t_lp'], want, 6e-4 if abs(want) > 0.005 else 6e-5)
true('odd X8: nothing feasible up to 12.15', all(d['t_lp'] < -0.002 for d in o8 if d['r'] <= 12.151))
e8 = S[(0, 8)]; check('even X8 t*(13.80)', at(e8, 13.80)['t_lp'], -0.008, 6e-4); true('even X8 nothing else below 17.70', all(d['t_lp'] < -0.002 for d in e8 if d['r'] < 17.699))
true('D2 misses', all(-0.0252 - 1e-4 <= x <= -0.0022 + 1e-4 for x in [max(d['t_lp'] for d in S[(e, X)] if abs(d['r'] - r1) <= 0.0501) for e, r1 in ((1, 9.53369526), (0, 13.77975135)) for X in (8, 16, 32)]))
def ap(e, X, r): d = at(S[(e, X)], r); return dict(zip(d['primes'], d['mean_ap']))
for e, r, wants in [(1, 9.5337, [-1.065, -0.464, -0.281, -0.738]), (1, 12.173, [0.277, -1.202, 0.038, 0.436]), (0, 13.7798, [1.542, 0.257, 0.725, -0.227])]:
    a = ap(e, 16, r)
    for p, want in zip((2, 3, 5, 7), wants): check('a_%d at %g' % (p, r), a[p], want, 6e-4)
for r, wants in [(9.53369526, [-1.068, -0.456, -0.291, -0.745]), (12.17300832, [0.289, -1.202, 0.040, 0.448]), (13.77975135, [1.549, 0.247, 0.737, -0.261])]:
    f = min(key, key=lambda f: abs(f['r'] - r))
    for p, want in zip((2, 3, 5, 7), wants): check('true a_%d at %g' % (p, r), f['ap'][str(p)], want, 6e-4)
f3 = min(key, key=lambda f: abs(f['r'] - 14.3585)); a = ap(1, 16, 14.3585); check('third odd form: largest difference', max(abs(a[p] - f3['ap'][str(p)]) for p in (2, 3, 5, 7)), 1.07, 6e-3)
def comps(rs, tau=0.002, step=0.0501):
    out = []
    for d in rs:
        if d['t_lp'] >= -tau:
            if out and d['r'] - out[-1][1] < step: out[-1][1] = d['r']
            else: out.append([d['r'], d['r']])
    return out
true('even X16 components', [[round(a, 2), round(b, 2)] for a, b in comps(S[(0, 16)])] == [[17.7, 17.9], [19.35, 19.4]])
check('even X16 t*(13.80)', at(S[(0, 16)], 13.80)['t_lp'], -0.011, 6e-4)
true('even feasible counts', [sum(d['t_lp'] >= -0.002 for d in S[(0, X)]) for X in (8, 16, 32)] == [29, 7, 3])
true('odd components start at 12.20', all([round(v, 2) for v in comps(S[(1, X)], step=0.1001)[-1]] == [12.2, 20.0] for X in (8, 16, 32)))
true('odd feasible counts', [sum(d['t_lp'] >= -0.002 for d in S[(1, X)]) for X in (8, 16, 32)] == [157, 157, 79])
act = sum(1 for X in (16, 32) for e in (0, 1) for d in S[(e, X)] if min(d['per_form'][2:]) < d['t_lp'] + 1e-4)
true('Sym3/Sym4 active at a handful of points only', act <= 12)
tw = txt('twists_out.txt').splitlines(); tv = [float(re.search(r't\* = ([+-]?[\d.]+)', l).group(1)) for l in tw]
true('twists: unchanged to 5 decimals', len(tv) == 8 and all(abs(tv[i] - tv[i + 4]) < 1e-5 for i in range(4)))
d13 = at(S[(1, 16)], 13.0); a2 = [x for x in d13['atoms'][0] if x[1] > 0.1]
true('two atoms at p = 2 at r = 13', min(x[0] for x in a2) < 1.0 and max(x[0] for x in a2) > 1.6)

# --- section 6: round 166 (genuine angles, incomplete)
pm = rows('pm_eps1_X16_k2.jsonl'); true('pm: 34 points', len(pm) == 34); check('pm seconds', sum(d['secs'] for d in pm), 6735, 1)
lo = [d for d in pm if 11.99 < d['r'] < 13.51]; true('pm 12.0-13.5: 16 points, all infeasible', len(lo) == 16 and all(d['t_upper'] < -0.002 for d in lo))
check('pm upper bounds, smallest', min(d['t_upper'] for d in lo), -0.031, 6e-4); check('pm upper bounds, largest', max(d['t_upper'] for d in lo), -0.0024, 6e-5)
conv = S[(1, 16)]; true('convex feasible on 12.2-13.5', all(d['t_lp'] >= -0.002 for d in conv if 12.199 < d['r'] < 13.501)); check('convex max there', max(d['t_lp'] for d in conv if 12.199 < d['r'] < 13.501), 0.21, 6e-3)
true('pm 16.0-16.6 infeasible (7)', sum(1 for d in pm if 15.99 < d['r'] < 16.61 and d['t_upper'] < -0.002) == 7)
hi = [d for d in pm if d['r'] > 16.65]; true('pm 16.7-17.7 feasible (11)', len(hi) == 11 and all(d['t_upper'] >= -0.002 for d in hi) and abs(hi[0]['r'] - 16.7) < 1e-9 and abs(hi[-1]['r'] - 17.7) < 1e-9)
check('pm a_2 at 16.7', hi[0]['ap'][0], -1.52, 6e-3); check('pm a_3 at 16.7', hi[0]['ap'][1], 0.92, 6e-3)
f = min(key, key=lambda f: abs(f['r'] - 16.6443)); check('true a_2', f['ap']['2'], -1.540, 6e-4); check('true a_3', f['ap']['3'], 0.977, 6e-4)
true('E3 on the part run: 11 against 34', sum(d['t_upper'] >= -0.002 for d in pm) == 11 and sum(1 for d in pm if at(conv, d['r'])['t_lp'] >= -0.002) == 34 - 2)

# --- section 7: round 167 (SL(3,Z))
g3 = json.load(open(os.path.join(D, 'gl3_key.json'))); true('gl3 key: 1050 forms', len(g3['forms']) == 1050)
k3 = sorted(g3['forms'], key=lambda f: f['laplace'])
loc = [json.loads(l) for l in open(os.path.join(D, 'gl3_local_X8.jsonl'))]; co = [json.loads(l) for l in open(os.path.join(D, 'gl3_coarse_X8.jsonl'))]
true('gl3 coarse: 861 points', len(co) == 861)
q = lambda r: r['l1']**2 + r['l1']*r['l2'] + r['l2']**2
low = [r for r in co if q(r) <= 180]; true('G3: 237 points, none feasible', len(low) == 237 and all(r['t_lp'] < -0.002 for r in low))
feas = [r for r in co if r['t_lp'] >= -0.002]; m = min(feas, key=q); check('smallest value on the feasible set', q(m), 222.25, 1e-9); true('at (13.5, 2.5)', (m['l1'], m['l2']) == (13.5, 2.5)); check('lowest form', k3[0]['laplace'] - 1, 239.31, 6e-3)
ex = []; near = []; d2 = []; d3 = []
for f in k3[:5]:
    l1, l2 = f['lambda']; e = [r for r in loc if r['kind'] == 'exact' and abs(r['l1'] - l1) < 1e-6 and abs(r['l2'] - l2) < 1e-6][0]; ex.append(e['t_lp'])
    grid = [r for r in loc if r['kind'] == 'grid' and abs(r['l1'] - round(l1, 1)) < 0.41 and abs(r['l2'] - round(l2, 1)) < 0.41]
    inf = [np.hypot(r['l1'] - l1, r['l2'] - l2) for r in grid if r['t_lp'] < -0.002]; near.append(min(inf) if inf else np.inf)
    d2.append(abs(complex(*e['a_p'][0]) - complex(*f['A_p1']['2']))); d3.append(abs(complex(*e['a_p'][1]) - complex(*f['A_p1']['3'])))
    if abs(l1 - 13.5966) < 1e-3: a2 = complex(*e['a_p'][0]); check('a_2 re', a2.real, 1.040, 6e-4); check('a_2 im', a2.imag, -0.382, 6e-4); check('true re', f['A_p1']['2'][0], 1.048, 6e-4); check('true im', f['A_p1']['2'][1], -0.375, 6e-4)
for got, want in zip(ex, [0.019, 0.167, 0.002, 0.030, 0.261]): check('gl3 t* at the forms', got, want, 6e-4)
check('form 3 t*', ex[2], 0.0016, 6e-5)
for got, want in zip(near, [0.185, 0.348, 0.035, 0.176]): check('gl3 nearest infeasible', got, want, 6e-4)
true('form 5 none', near[4] == np.inf)
for got, want in zip(d2, [0.081, 1.31, 0.011, 0.058, 0.79]): check('gl3 a_2 differences', got, want, 6e-3 if want > 0.5 else 6e-4)
check('gl3 a_3 difference, form 3', d3[2], 0.048, 6e-4)
true('both forms active at 844 points', sum(1 for r in co if max(r['per_form']) < r['t_lp'] + 1e-4) == 844)
true('self-dual line begins at 16.5', min(r['l1'] for r in co if r['l2'] == 0 and r['t_lp'] >= -0.002) == 16.5)

# --- section 8: round 168 (central zero) and its controls
tc = txt('truth_center_out.txt').splitlines(); odd = [l for l in tc if ' odd ' in l]
true('truth with central zero: 7 odd forms >= -2.6e-14', len(odd) == 7 and min(min(float(x) for x in re.findall(r'(?:f|S2) ([+-][\d.e+-]+)', l)) for l in odd) >= -2.6e-14)
C = {X: rows('scan_eps1_X%d_k2_J61_center.jsonl' % X) for X in (8, 16)}
true('168 grids', all(len(C[X]) == 1101 for X in C))
cc = lambda rs: [[round(a, 2), round(b, 2)] for a, b in comps(rs, step=0.0101)]
true('168 pieces X8', cc(C[8]) == [[9.53, 9.54], [12.17, 12.19], [14.28, 14.43], [15.82, 20.0]])
true('168 pieces X16', cc(C[16]) == [[9.54, 9.54], [12.17, 12.19], [14.35, 14.37], [15.93, 20.0]])
odd_true = [f for f in key if f['parity'] == 'odd' and f['r'] < 20]; true('seven odd', len(odd_true) == 7)
feas16 = [d['r'] for d in C[16] if d['t_lp'] >= -0.002]; true('415 feasible', len(feas16) == 415)
true('384 far from every eigenvalue', sum(1 for x in feas16 if min(abs(x - f['r']) for f in odd_true) > 0.03) == 384)
for f, wants in zip(odd_true[:3], [(-1.068, -0.455, -1.068, -0.456), (0.288, -1.201, 0.289, -1.202), (-0.230, 0.694, -0.231, 0.696)]):
    d = at(C[16], f['r']); got = dict(zip(d['primes'], d['mean_ap']))
    check('168 a_2', got[2], wants[0], 6e-4); check('168 a_3', got[3], wants[1], 6e-4); check('168 true a_2', f['ap']['2'], wants[2], 6e-4); check('168 true a_3', f['ap']['3'], wants[3], 6e-4)
for f in odd_true[3:]:
    d = at(C[16], f['r']); got = dict(zip(d['primes'], d['mean_ap'])); true('not returned above 15.9', abs(got[2] - f['ap']['2']) > 0.05)
H = {X: rows('hol_X%d_center.jsonl' % X) for X in (8, 16)}
for X, k, want, tol in [(8, 14.0, -0.470, 6e-4), (8, 17.75, -0.026, 6e-4), (8, 18.0, 0.00001, 6e-6), (8, 18.25, 0.026, 6e-4), (16, 17.75, -0.027, 6e-4), (16, 18.0, -0.0002, 6e-5)]: check('hol center', at(H[X], k)['t_lp'], want, tol)
true('weight 18 is the first', all(d['t_lp'] < -0.002 for X in H for d in H[X] if d['r'] < 17.99) and all(d['t_lp'] >= -0.002 for X in H for d in H[X] if d['r'] > 17.99))
for got, want in zip(at(H[8], 18.0)['mean_ap'][:3], [-1.4585, -0.3771, -1.1750]): check('w18 X8', got, want, 6e-5)
for got, want in zip(at(H[16], 18.0)['mean_ap'][:4], [-1.4575, -0.3783, -1.1697, 0.2032]): check('w18 X16', got, want, 6e-5)
import cypari2
pari = cypari2.Pari(); co = pari.mfcoefs(pari.mfeigenbasis(pari.mfinit([1, 18, 1], 0))[0], 8)
true('weight 18 coefficients', [int(co[p]) for p in (2, 3, 5, 7)] == [-528, -4284, -1025850, 3225992])
for p, want in zip((2, 3, 5, 7), [-1.4584, -0.3770, -1.1745, 0.2115]): check('w18 true', int(co[p])/p**8.5, want, 6e-5)
lc = txt('levels_center_out.txt')
for X, want in [(8, 32.46), (16, 34.76)]: check('level 37 bound', float(re.search(r'weight 2 X %d .*?smallest level ([\d.]+)' % X, lc).group(1)), want, 6e-3)
check('6.0 percent', 100*(1 - 34.76301/37), 6.0, 0.06)
def wN(N):
    mf = pari.mfinit([N, 2, 1], 0)
    if int(pari.mfdim(mf)) == 0: return []
    return [int(x) for orb in pari.mfatkineigenvalues(mf, N) for x in orb]
true('weight 2: root number -1 (W_N = +1) first occurs at level 37 (PARI)', all(all(w == -1 for w in wN(N)) for N in range(11, 37)) and 1 in wN(37))
E16 = rows('scan_eps0_X16_k2_J61.jsonl'); E32 = rows('scan_eps0_X32_k2_J61.jsonl'); true('even fine grid: 701 points', len(E16) == 701)
true('even pieces X16', cc(E16) == [[13.78, 13.78], [17.69, 17.94], [19.32, 19.43]]); true('even pieces X32', cc(E32) == [[13.78, 13.78], [17.7, 17.88], [19.38, 19.43]])
ev = [f for f in key if f['parity'] == 'even' and f['r'] < 20]
for f, wants, twants in zip([ev[0], ev[2]], [(1.549, 0.247, 0.736), (-0.698, 1.555, -0.039)], [(1.549, 0.247, 0.737), (-0.693, 1.562, -0.038)]):
    d = at(E16, f['r']); got = dict(zip(d['primes'], d['mean_ap']))
    for p, w, tw in zip((2, 3, 5), wants, twants): check('even fine a_p', got[p], w, 6e-4); check('even true a_p', f['ap'][str(p)], tw, 6e-4)
d = at(E16, ev[1]['r']); true('17.74 not returned', abs(dict(zip(d['primes'], d['mean_ap']))[2] - ev[1]['ap']['2']) > 0.2)
true('six isolated: the pieces at support 16 hold one eigenvalue each', all(sum(1 for f in ev if a - 0.011 <= f['r'] <= b + 0.011) == 1 for a, b in cc(E16)) and all(sum(1 for f in odd_true if a - 0.011 <= f['r'] <= b + 0.011) == 1 for a, b in cc(C[16])[:3]))

# --- In short, item 3
five = []
for f in [x for x in key if x['r'] < 20]:
    d = at(C[16], f['r']) if f['parity'] == 'odd' else at(E16, f['r']); got = dict(zip(d['primes'], d['mean_ap']))
    five.append(max(abs(got[2] - f['ap']['2']), abs(got[3] - f['ap']['3'])))
true('five of the six isolated forms return a_2, a_3 within 0.01', sorted(five)[:5][-1] <= 0.01 and sorted(five)[5] > 0.2)

# --- section 8b: round 170 (whole-number a_p)
it = txt('integral_out.txt').splitlines(); true('integral: 6 lines', len(it) == 6)
lev = [float(re.search(r'smallest level ([\d.]+)', l).group(1)) for l in it]
for got, want in zip(lev, [10.947, 10.992, 10.995, 36.260, 36.710, 36.906]): check('integral level', got, want, 6e-4)
def top(line, p):
    m = re.search(r' %d:(\S+)' % p, line.split('a_p(weight):')[1]); best = max(re.findall(r'([+-]\d+)S?\(([\d.]+)\)', m.group(1)), key=lambda x: float(x[1])); return int(best[0]), float(best[1])
for line, wants in [(it[2], [(-2, 1.00), (-1, 0.99), (1, 0.98), (-2, 0.98)]), (it[5], [(-2, 1.00), (-3, 1.00), (-2, 0.89), (-1, 0.73)])]:
    for p, (a, w) in zip((2, 3, 5, 7), wants): g = top(line, p); true('integral a_%d' % p, g[0] == a and abs(g[1] - w) < 0.006)
E11 = pari.mfcoefs(pari.mfeigenbasis(pari.mfinit([11, 2, 1], 0))[0], 20); true('11a', [int(E11[p]) for p in (2, 3, 5, 7, 13, 19)] == [-2, -1, 1, -2, 4, 0])
E37 = [pari.mfcoefs(F, 8) for F in pari.mfeigenbasis(pari.mfinit([37, 2, 1], 0))]; true('37a', any([int(c[p]) for p in (2, 3, 5, 7)] == [-2, -3, -2, -1] for c in E37))
true('a_13, a_19 at support 32 are not the curve\'s', top(it[2], 13)[0] in (5, 6) and top(it[2], 19)[0] == -2)
check('0.05 percent', 100*(1 - lev[2]/11), 0.05, 0.006); check('0.25 percent', 100*(1 - lev[5]/37), 0.25, 0.006)
true('I2 fails at support 16 on a_7', top(it[4], 7)[0] in (0, 1))
check('unrestricted bounds', float(re.search(r'weight 2 X 16 .*?smallest level ([\d.]+)', lv).group(1)), 10.936, 6e-4)

# --- section 8c: prescribed values
cl = txt('chars_cond_out.txt').splitlines(); true('18 cases', len(cl) == 18)
ratio = [float(re.search(r'ratio ([\d.]+)', l).group(1)) for l in cl]; cond = [int(re.search(r'conductor (\d+)', l).group(1)) for l in cl]; mt = [int(re.search(r'matched at (\d+) of 12', l).group(1)) for l in cl]
ex = [(c, m) for r_, c, m in zip(ratio, cond, mt) if r_ > 0.99999]
true('8 exact cases', sorted(c for c, m in ex) == sorted([3, 4, 7, 8, 5, 8, 12, 13]) and all(m >= 11 for c, m in ex))
true('17 nearly', any(c == 17 and abs(float(re.search(r'bound ([\d.]+)', l).group(1)) - 16.995) < 6e-4 for c, l in zip(cond, cl)))
true('nine far below', sum(1 for r_ in ratio if r_ < 0.95) == 9)
for c, want in [(19, 4.98), (24, 5.07), (21, 7.28), (11, 6.30), (73, 68.9)]:
    l = [x for x, cc_ in zip(cl, cond) if cc_ == c][0]; check('cond bound %d' % c, float(re.search(r'bound ([\d.]+)', l).group(1)), want, 0.06 if c == 73 else 6e-3)
dl = txt('chars_depth_out.txt').splitlines(); true('29 discriminants up to 44', len(dl) == 29)
depth = {int(re.match(r'D =\s*(-?\d+)', l).group(1)): re.search(r'k = (\w+)', l).group(1) for l in dl}
want = {'0': [-3, 5], '1': [-4, -7, 8], '2': [-8, 12, 13], '3': [-11, -15, 24], '4': [-20, 21], '5': [-24], '6': [17, -23]}
true('depth table', all(sorted(D for D, k in depth.items() if k == kk) == sorted(v_) for kk, v_ in want.items()) and sum(1 for k in depth.values() if k == 'None') == 13)
for D_, want_ in [(-19, 18.91), (28, 27.91), (37, 36.99), (-43, 19.7)]:
    l = [x for x in dl if re.match(r'D =\s*%d:' % D_, x)][0]; check('depth last value %d' % D_, float(re.findall(r'k=6: ([\d.]+)', l)[0]), want_, 0.06 if D_ == -43 else 6e-3)

# --- corrections after the second review
true('max |a_p| of the key forms in truth_out <= 1.74', max(float(re.search(r'max\|a_p\| ([\d.]+)', l).group(1)) for l in tr) <= 1.74)
for r_, want in [(9.54, -0.0013), (12.17, -0.0006), (14.36, -0.0002)]: check('t* at the island grid points', at(C[16], r_)['t_lp'], want, 6e-5)
check('x_7 at support 8, odd', float(re.search(r'kappa 1 X 8:.*? 7:([+-][\d.]+)', txt('chars_precise_out.txt')).group(1)), 0.9185, 6e-5)

# --- section 8d: round 169
p9 = rows('pm_eps1_X16_k2_center.jsonl'); true('169: 25 points', len(p9) == 25); check('169 seconds', sum(d['secs'] for d in p9), 3268, 1)
inf9 = [d for d in p9 if d['t_upper'] < -0.002]; true('169 infeasible 15.9-16.6', [round(d['r'], 1) for d in inf9] == [15.9, 16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6])
check('169 smallest upper', min(d['t_upper'] for d in inf9), -0.032, 6e-4); check('169 largest upper', max(d['t_upper'] for d in inf9), -0.0025, 6e-5)
fe9 = [d for d in p9 if d['t_upper'] >= -0.002]; true('169 feasible 16.7-18.3 (17)', len(fe9) == 17 and abs(fe9[0]['r'] - 16.7) < 1e-9 and abs(fe9[-1]['r'] - 18.3) < 1e-9)
check('169 a_2', fe9[0]['ap'][0], -1.53, 6e-3); check('169 a_3', fe9[0]['ap'][1], 0.95, 6e-3)
true('168 convex feasible from 15.93', abs(min(d['r'] for d in C[16] if d['r'] > 15 and d['t_lp'] >= -0.002) - 15.93) < 1e-9)

# --- section 8e: round 171
mt_ = txt('mixed_truth_out.txt').splitlines(); true('171 truth: 7 odd forms', len(mt_) == 7 and min(float(x) for l in mt_ for x in re.findall(r' ([+-]\d\.\de[+-]\d\d)', l)) >= -2.6e-14)
true('171 truth: up to nine known forms', max(int(re.search(r'(\d+) known forms', l).group(1)) for l in mt_) == 9)
MO = rows('mixed_eps1_X16.jsonl'); ME = rows('mixed_eps0_X16.jsonl'); true('171 grids', len(MO) == 411 and len(ME) == 701)
true('171 odd pieces', cc(MO) == [[16.14, 16.18], [16.65, 16.69], [18.18, 18.23], [19.49, 19.53]]); true('171 even pieces', cc(ME) == [[13.78, 13.78], [17.74, 17.74], [19.43, 19.43]])
true('171 counts', sum(d['t_lp'] >= -0.002 for d in MO) == 21 and sum(d['t_lp'] >= -0.002 for d in ME) == 3 and sum(1 for d in C[16] if d['r'] >= 15.899 and d['t_lp'] >= -0.002) == 408 and sum(d['t_lp'] >= -0.002 for d in E16) == 39)
tab = [(16.1381, MO, 1.156, 1.162, -1.279, -1.282), (16.6443, MO, -1.539, -1.540, 0.973, 0.977), (18.1809, MO, 0.374, 0.374, 0.101, 0.102), (19.4847, MO, -1.681, -1.700, -0.610, -0.615),
       (13.7798, ME, 1.549, 1.549, 0.247, 0.247), (17.7386, ME, -0.767, -0.765, -0.975, -0.978), (19.4235, ME, -0.686, -0.693, 1.563, 1.562)]
worst = 0
for rj, M, a2, t2, a3, t3 in tab:
    f = min(key, key=lambda f: abs(f['r'] - rj)); d = [x for x in M if x['r'] >= rj - 0.005 and x['t_lp'] >= -0.002][0]; got = dict(zip(d['primes'], d['mean_ap']))
    check('171 a_2 %g' % rj, got[2], a2, 6e-4); check('171 a_3 %g' % rj, got[3], a3, 6e-4); check('171 true a_2', f['ap']['2'], t2, 6e-4); check('171 true a_3', f['ap']['3'], t3, 6e-4)
    worst = max(worst, abs(got[2] - f['ap']['2']), abs(got[3] - f['ap']['3']))
true('all ten: a_2, a_3 within 0.02', worst <= 0.02 and max(five[i] for i in (0, 1, 2, 3)) <= 0.02)
allp = cc(C[16])[:3] + cc(MO) + cc(ME)
true('ten pieces, each at most 0.05 long, one eigenvalue each', len(allp) == 10 and all(b - a <= 0.0501 for a, b in allp) and all(sum(1 for f in key if f['r'] < 20 and a - 0.011 <= f['r'] <= b + 0.011) == 1 for a, b in allp))

# --- section 8f: the certificate note (its own validation is theory/validate_exact3.py; here: the quoted numbers are those of the note)
ex3 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'theory', 'exact3.md')).read()
true('exact3: own validation passed', '157 checks, 0 failed' in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'theory', 'validate_exact3_out.txt')).read())
for frag in ['| 8 | [5.136, 5.137]·10⁻⁷ | 3.082·10⁻⁶ |', '| 12 | [7.479, 7.480]·10⁻¹¹ | 4.488·10⁻¹⁰ |', '| 16 | [1.510, 1.511]·10⁻¹⁷ | 9.062·10⁻¹⁷ |', '| 22 | [6.349, 6.350]·10⁻²³ | 3.810·10⁻²² |', 'B = 1 gives the bound 2.740045']:
    true('exact3 table row', frag in ex3)
true('exact3 tail bound quoted', '10⁻²⁹⁷' in ex3 and '600' in ex3)

# --- section 8g: round 172 (the self-contained chain)
fo = json.load(open(os.path.join(D, 'chain_X16_found.json'))); ch = [json.loads(l) for l in open(os.path.join(D, 'chain_X16.jsonl'))]
true('172: 505 points, ten forms', len(ch) == 505 and len(fo) == 10); true('172: 1,700 s', '1700 s' in txt('chain_X16_out.txt'))
tab72 = [(9.53443, 1, 9.53370, -1.0679, -1.0683, -0.4566, -0.4562, 0.0005, 0.0047), (12.17376, 1, 12.17301, 0.2895, 0.2893, -1.2017, -1.2019, 0.0005, 0.0120), (13.77916, 0, 13.77975, 1.5494, 1.5493, 0.2471, 0.2469, 0.0007, 0.0048),
         (14.35814, 1, 14.35851, -0.2313, -0.2309, 0.6960, 0.6956, 0.0004, 0.0049), (16.13660, 1, 16.13807, 1.1613, 1.1619, -1.2816, -1.2820, 0.0022, 0.0050), (16.64339, 1, 16.64426, -1.5397, -1.5402, 0.9769, 0.9775, 0.0011, 0.0114),
         (17.73718, 0, 17.73856, -0.7657, -0.7655, -0.9771, -0.9778, 0.0056, 0.0138), (18.17936, 1, 18.18092, 0.3745, 0.3741, 0.1012, 0.1020, 0.0008, 0.0056), (19.42290, 0, 19.42348, -0.6939, -0.6928, 1.5610, 1.5624, 0.0021, 0.0055),
         (19.48992, 1, 19.48471, -1.6780, -1.7002, -0.6087, -0.6146, 0.0331, 0.0618)]
errs_r = []; errs_a = []
for g, (rh, e, rt, a2, t2, a3, t3, e7, e13) in zip(fo, tab72):
    f = min([x for x in key if (x['parity'] == 'odd') == (g['eps'] == 1)], key=lambda x: abs(x['r'] - g['r']))
    check('172 r^', g['r'], rh, 6e-6); true('172 parity', g['eps'] == e); check('172 true r', f['r'], rt, 6e-6)
    check('172 a_2', g['ap']['2'], a2, 6e-5); check('172 true a_2', f['ap']['2'], t2, 6e-5); check('172 a_3', g['ap']['3'], a3, 6e-5); check('172 true a_3', f['ap']['3'], t3, 6e-5)
    check('172 err p<=7', max(abs(g['ap'][p] - f['ap'][p]) for p in '2357'), e7, 6e-5); check('172 err 11,13', max(abs(g['ap'][p] - f['ap'][p]) for p in ('11', '13')), e13, 6e-5)
    errs_r.append(abs(g['r'] - f['r'])); errs_a.append(max(abs(g['ap']['2'] - f['ap']['2']), abs(g['ap']['3'] - f['ap']['3'])))
check('172 largest r error', max(errs_r), 0.0052, 6e-5); check('172 r error, first nine', max(errs_r[:9]), 0.0016, 6e-5)
check('172 a_2, a_3 first nine', max(errs_a[:9]), 0.0014, 6e-5); check('172 tenth', errs_a[9], 0.022, 6e-4)
fe72 = [d for d in ch if d['t_lp'] >= -0.002]; true('172: 28 feasible points', len(fe72) == 28)
true('172 N3', all(any(f['r'] - 0.02 - 1e-9 <= d['r'] <= f['r'] + 0.15 + 1e-9 for f in key if f['r'] < 20 and (f['parity'] == 'odd') == (d['eps'] == 1)) for d in fe72))
true('172 last piece [19.48, 19.53]', sorted(d['r'] for d in fe72 if d['eps'] == 1 and d['r'] > 19)[0] == 19.48 and sorted(d['r'] for d in fe72 if d['eps'] == 1)[-1] == 19.53)

ex5 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'theory', 'exact5.md')).read()
true('exact5: own validation passed', '137 checks, 0 failed' in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'theory', 'validate_exact5_out.txt')).read())
true('exact5 rows', all(x in ex5 for x in ['2.207·10⁻⁶', '1.148·10⁻⁹', '2.764·10⁻¹⁴']))

# --- END
if __name__ == '__main__':
    print('%d checks pass, %d fail' % (ok, bad)); sys.exit(1 if bad else 0)

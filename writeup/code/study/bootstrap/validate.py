# Checks every number quoted in bootstrap.md against the output files of the scripts in this folder.
import re, os
here = os.path.dirname(os.path.abspath(__file__))
rd = lambda f: open(os.path.join(here, f)).read()
n = 0
def ok(cond, what):
    global n; n += 1
    assert cond, what
# section 1 and 2
t = rd('ef_check_out.txt'); z = float(re.search(r'zeros\s+(\S+)', t).group(1)); r = float(re.search(r'rhs\s+(\S+)', t).group(1))
zs = re.search(r'zeros\s+(\S+)', t).group(1); rs = re.search(r'rhs\s+(\S+)', t).group(1)
from mpmath import mp, mpf
mp.dps = 40
ok(abs(mpf(zs) - mpf(rs))/mpf(zs) < mpf(10)**(-18) and abs(mpf(zs) - mpf(rs)) < mpf(10)**(-27), 'explicit formula to 27 digits')
t = rd('cc_check_out.txt')
ok('mu=3 N=48 even 5.6799e-8' in t and float(re.search(r'mu=3 N=48 even (\S+)', t).group(1)) < 6e-8, 'mu=3')
ok('mu=11 N=48 even 2.9274e-48' in t, 'mu=11')
v = [float(x) for x in re.findall(r'mu=11 N=\d+ even (\S+)', t)]; ok(v[0] > v[1] > v[2] > 2.389e-48, 'still falling, above the published value')
ok(all(float(x) < 0 for x in re.findall(r'one point at \S+ : lam_min (\S+)', t)), 'negative at 1.9999 and 2.0005')
# 3.1
t = rd('scan_out.txt'); a, b = t.split('X 4.0')
ok('empty -0.4156' in a and 'feasible grid points: [2.]' in a, 'X=3 one point')
f3 = re.search(r'feasible \(> -1e-9\): (.*)', a).group(1); ys = [float(x) for x in re.findall(r'np\.float64\(([\d.]+)\)', f3)]
ok(set(ys[0::2]) == {2.0} and min(ys[1::2]) == 2.96, 'X=3 two points')
ok('empty -0.7036' in b and 'feasible grid points: []' in b and '-0.2415' in b, 'X=4 one point')
ok(re.search(r'feasible \(> -1e-9\): \[\(np\.float64\(2\.0\), np\.float64\(3\.0\)\)\]', b) is not None, 'X=4 only (2,3)')
ok(abs(float(re.search(r'3 points: top \[\(np\.float64\((-[\d.]+)\)', b).group(1)) + 0.00718) < 1e-5, 'X=4 three points')
# 3.2
t = rd('pins_out.txt'); blocks = t.split('\nX ')
want = {'3.0': ('8.46e-8', {2: (4, 6)}), '4.0': ('1.11e-12', {2: (8, 11), 3: (7, 8)}),
        '6.0': ('4.21e-22', {2: (18, 20), 3: (16, 18), 5: (9, 10)}),
        '8.0': ('3.3e-29', {2: (25, 27), 3: (23, 25), 5: (18, 20), 7: (10, 11)}),
        '11.0': ('2.15e-43', {2: (38, 42), 3: (37, 40), 5: (33, 34), 7: (27, 28)})}
for blk in blocks:
    X = blk.replace('X ', '').split()[0]
    ev, pins = want[X]
    ok("['%s'" % ev in blk, 'eigenvalue at X=' + X)
    for p, (u, d) in pins.items():
        ok('prime %d: infeasible down to shift +1e-%d, -1e-%d' % (p, u, d) in blk, 'pin %s %d' % (X, p))
# 3.3
t = rd('sequential_out.txt'); lev = {}
for m in re.finditer(r'support X=(\S+) .*\n((?:     .*\n)+)', t):
    lev[m.group(1)] = [(float(l.split()[0]), float(l.split('=')[-1])) for l in m.group(2).strip().split('\n')]
ok(abs(lev['2.5'][0][0] - 2.0027) < 1e-4, 'X=2.5')
ok([d for _, d in lev['3.5']] == [1.18e-6, 0.000266], 'X=3.5')
ok([d for _, d in lev['4.5']] == [1.23e-10, 6.6e-9], 'X=4.5')
ok([d for _, d in lev['6']] == [3.76e-17, 3.69e-15, 3.22e-9], 'X=6')
ok([d for _, d in lev['9']] == [1.01e-18, 2.62e-18, 1.6e-16, 3.36e-12], 'X=9')
ok([d for _, d in lev['10.5']] == [2.98e-17, 8.83e-17, 1.41e-16, 3.97e-15], 'X=10.5')
l = lev['13.5']; ok(all(d < 5e-9 for _, d in l[:4]) and l[4][1] == 9.74e-8 and abs(l[5][0] - 13.002) < 1e-3 and len(l) == 10, 'X=13.5')
ok(all(round(x) == p for (x, _), p in zip(lev['10.5'], (2, 3, 5, 7))), 'primes 2 3 5 7')
res = [float(x) for x in re.findall(r'final lam_min=(\S+)', t)][:7]
ok(min(res[3:]) > -1e-13 and res[2] == 1.85e-11, 'residual eigenvalues')
# 4
t = rd('datum_out.txt'); ok('X=9     lam_min  5.17e-02  points: -1.6630204321' in t, 'q=4 odd')
t = rd('characters_check_out.txt'); ok('chi_-4 X 9 true character lam_min 7.330e-09' in t and '(0.1776, 1.24, -1.0, 1.98, 1.0)' in t, 'characters, signs')
from math import pi, exp
ok(abs(pi*exp(1.0858608) - 9.3) < 0.01 and abs(pi*exp(4.2274535) - 215.4) < 0.1, 'empty-system conductors: pi exp(-psi(3/4)), pi exp(-psi(1/4))')
# 5
t = rd('conductor_out.txt')
def row(k, X): 
    m = re.search(r'kappas %s X=%-4g max lam found \S+\s+-> smallest conductor found (\S+)\s+points: (.*)' % (k, X), t)
    return float(m.group(1)), [float(x) for x in m.group(2).split()]
for X, q in ((4, 2.976788), (6, 2.994999), (8, 2.999540), (10, 3.001629)): ok(row('01', X)[0] == q, 'Gamma_C %d' % X)
for X, q in ((4, 4.965711), (6, 4.979974), (8, 4.991092), (10, 4.997454)): ok(row('00', X)[0] == q, 'Gamma_R^2 %d' % X)
for X in (4, 6, 8): ok(row('0', X)[0] == 1.0, 'zeta %d' % X)
ok([round(x) for x in row('0', 8)[1]] == [2, 3, 5, 7] and max(abs(x - round(x)) for x in row('0', 8)[1]) < 1e-5, 'zeta points')
ok([round(x, 3) for x in row('01', 8)[1]] == [3.002, 3.996, 6.831] and [round(x, 3) for x in row('00', 10)[1]] == [4.018, 4.963], 'points')
qs = [row('01', X)[0] for X in (4, 5, 6, 7, 8)]; ok(qs == sorted(qs) and qs[-1] < 3, 'rising below 3')
qs = [row('00', X)[0] for X in (4, 5, 6, 7, 8, 9, 10)]; ok(qs == sorted(qs) and qs[-1] < 5, 'rising below 5')
t = rd('conductor_true_out.txt')
for name, vals in (('Q(sqrt-3)', (2.998013, 2.999988, 3.0, 3.0)), ('Q(i)     ', (3.960526, 3.999689, 3.99999, 4.0)), ('Q(sqrt 5)', (4.965711, 4.999271, 4.999977, 5.0)), ('Q(sqrt 2)', (7.244072, 7.940275, 7.996105, 7.999969))):
    got = [float(x) for x in re.findall(re.escape(name) + r' d=\d+ X=(?:4|6|8|12)\s+true norms: conductor (\S+)', t)]
    ok(got == list(vals), name)
for name, q in (('Q(sqrt-3)', 3.000014), ('Q(i)     ', 3.000071), ('Q(sqrt 5)', 4.997769), ('Q(sqrt 2)', 5.000372)):
    ok(re.search(re.escape(name) + r' d=\d+ X=12 .*maximised: conductor %.6f' % q, t) is not None, name + ' search')
print('all %d checks passed' % n)

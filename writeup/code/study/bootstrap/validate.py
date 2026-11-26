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
t = rd('conductor_cubic_out.txt')
got = [float(x) for x in re.findall(r'smallest conductor found (\S+)', t)]
ok([round(x, 4) for x in got] == [21.0020, 22.5285, 22.5525, 22.5935, 22.5983, 22.6300] and got == sorted(got) and got[-1] < 23 and got[1] > 22.27, 'cubic')
ok(t.strip().split('\n')[-1].split('points:')[1].split() == ['5.94013', '5.94014', '7.33387'], 'cubic points')
# 5b
from math import log, sqrt
t = rd('edit_cost_out.txt')
rows = {m.group(1).strip(): ([float(x) for x in m.group(2).split()], float(m.group(3)), float(m.group(4))) for m in re.finditer(r'(.{14}) lambda_min at X = 8, 16, 32, 64: (.*?)   limit 2 inf h = (\S+)   conductor exp\(-limit\) = (\S+)', t)}
want = {'remove 2': ([-0.5409, -0.5563, -0.5629, -0.5664], -0.5742, 1.776), 'remove 3': ([-0.6343, -0.7310, -0.7452, -0.7643], -0.8042, 2.235),
        'remove 7': ([-0.4030, -0.7354, -0.7355, -0.8811], -1.0675, 2.908), 'add 6': ([-0.6177, -0.7315, -0.7315, -1.1923], -2.4723, 11.849),
        'move 2 to 2.1': ([-0.5879, -1.0625, -1.3723, -1.8119], -3.8780, 48.330)}
for k, w in want.items():
    ok(rows[k] == w, 'edit ' + k)
    ok(all(a >= b for a, b in zip(rows[k][0], rows[k][0][1:])) and rows[k][0][-1] > rows[k][1], 'falls towards the limit: ' + k)
for p_ in (2, 3, 7): ok(abs(-2*log(p_)/(sqrt(p_) + 1) - want['remove %d' % p_][1]) < 1e-4 and abs(p_**(2/(sqrt(p_) + 1)) - want['remove %d' % p_][2]) < 1e-3, 'closed form remove')
ok(abs(-2*log(6)/(sqrt(6) - 1) - want['add 6'][1]) < 1e-4 and abs(6**(2/(sqrt(6) - 1)) - 11.849) < 2e-3, 'closed form add')
ok(abs(-2*log(2)/(sqrt(2) + 1) - 2*log(2.1)/(sqrt(2.1) - 1) - want['move 2 to 2.1'][1]) < 2e-3, 'prices add')
t = rd('empty_check_out.txt')
xs = [float(x) for x in re.findall(r'positive up to X = (\S+)', t)]
ok(xs == [2.11508, 2.10917, 2.10709] and xs[-1] > 2.1027 and xs[-1] < 2.11, 'empty system threshold')
last = [float(l.split()[-1]) for l in t.strip().split('\n')[3:]]
ok(all(abs(x - y) < 1e-3 for x, y in zip(last, [-0.012, -0.050, -0.086, -0.122, -0.109])), 'next prime forced')
t = rd('relaxed_out.txt')
lam = [float(x) for x in re.findall(r'max lambda_min over real weights = (\S+)', t)]
ok(lam == [7.009e-04, 3.248e-05, 1.808e-05, 1.612e-06] and lam == sorted(lam, reverse=True), 'relaxed optimum')
ws = [eval(l.split('weight):')[1]) for l in t.split('\n') if 'weight):' in l]
tot = lambda w, lo, hi: sum(v for x, v in w if lo <= x <= hi)
ok(abs(tot(ws[1], 1.9, 2.1) - 1.003) < 1e-3 and abs(tot(ws[3], 1.9, 2.1) - 1.0001) < 1e-4 and abs(tot(ws[2], 2.9, 3.1) - 1.197) < 1e-3 and abs(tot(ws[3], 2.9, 3.1) - 1.005) < 1e-3, 'relaxed weights')
ok(all(1.99 < x < 2.04 or 2.99 < x < 3.05 for w in ws for x, _ in w), 'weights sit at 2 and 3')
ok('conductor 0.99999839' in t, 'relaxed conductor')
t = rd('relaxed_fields_out.txt')
def col(k):
    return {float(m.group(1)): (float(m.group(2)), m.group(3)) for m in re.finditer(r'kappas %s X=(\S+)\s+N=\d+: lambda \S+ conductor (\S+) .*?centre\(weight\): (.*)' % k, t)}
z = col('0')
ok(z[3.0][0] == 0.9999677 and z[4.0][0] == 0.9999984 and all(abs(z[X][0] - 1) < 5e-8 for X in (5.0, 6.0, 7.0, 8.0)), 'convex, zeta')
ok(z[8.0][1].split() == ['2.000(1.00)', '3.000(1.00)', '5.000(1.00)', '7.000(1.00)'], 'convex, zeta weights')
c = col('01'); vals = [c[X][0] for X in (4.0, 6.0, 8.0, 10.0, 12.0, 14.0)]
ok([round(x, 5) for x in vals] == [2.97588, 2.98099, 2.98117, 2.98123, 2.98124, 2.98126] and vals == sorted(vals) and vals[-1] < 2.9921, 'convex, Gamma_C')
ok(c[14.0][1].split() == ['3.244(1.65)', '6.259(1.92)', '8.862(0.65)'], 'convex, Gamma_C weights')
c = col('00'); vals = [c[X][0] for X in (4.0, 6.0, 8.0, 10.0, 12.0, 14.0)]
ok([round(x, 5) for x in vals] == [4.96578, 4.97863, 4.98386, 4.98476, 4.98480, 4.98485] and c[14.0][1].split() == ['4.320(1.83)', '9.224(2.30)'], 'convex, Gamma_R^2')
c = col('001'); vals = [c[X][0] for X in (6.0, 8.0, 10.0, 12.0, 14.0)]
ok([round(x, 4) for x in vals] == [22.5288, 22.5523, 22.5664, 22.5852, 22.5953] and c[14.0][1].split() == ['6.206(2.82)'], 'convex, cubic complex')
c = col('000'); vals = [c[X][0] for X in (6.0, 8.0, 10.0, 12.0, 14.0)]
ok([round(x, 4) for x in vals] == [46.5468, 48.0437, 48.1272, 48.1324, 48.1486] and c[14.0][1].split() == ['8.931(4.21)'], 'convex, cubic real')
t = rd('relaxed_k1_out.txt')
lam = [float(x) for x in re.findall(r'lambda (\S+) \(LP', t)]
ok([('%.1e' % x) for x in lam[:3]] == ['3.2e-05', '1.6e-06', '9.9e-08'] and '4.004(0.507)' in t and '4.000(0.500)' in t and '2.000(1.000) 3.000(1.000)' in t, 'k=1 class: same optimum, von Mangoldt weights')
t = rd('relaxed_dual_out.txt')
gv = {float(m.group(1)): float(m.group(2)) for m in re.finditer(r'x=(\S+)\s+g\(log x\)/g\(0\) =\s+(\S+)', t)}
ok(abs(gv[2.504] - 1.586e-3) < 1e-6 and abs(gv[3.502] - 4.499e-5) < 1e-8 and abs(gv[4.503] - 2.488e-6) < 1e-9 and abs(gv[6.0] - 2.264e-7) < 1e-10 and all(abs(gv[x]) < 1e-9 for x in (2.0, 3.0, 4.0, 5.0, 7.0)), 'dual solution: bumps and zeros')
t = rd('magic_check_out.txt')
ours = [float(x) for x in re.search(r'a_1..a_5 \(ours\)  : (.*)', t).group(1).split()]; theirs = [float(x) for x in re.search(r'a_1..a_5 \(paper\) : (.*)', t).group(1).split()]
ok(all(abs(a/b - 1) < 6e-9 for a, b in zip(ours, theirs)), 'magic function: coefficients a_1..a_5')
m = re.search(r'pole term / g\(0\) = (\S+)\s+archimedean term / g\(0\) = (\S+)\s+W\(g\)/g\(0\) = (\S+)', t)
ok(m.group(1) == '1.120099144' and m.group(2) == '-1.120099143' and abs(float(m.group(3))) < 2e-9, 'magic function: W = 0 to 2e-9')
ok(abs(float(re.search(r'min over t in \[0, 60\] \(step 0.1\): (\S+)', t).group(1))) < 1e-11, 'magic function: transform >= 0')
zs = [(float(a), float(b), float(c)) for a, b, c in re.findall(r'zero \S+ : (\S+)\s+midpoints -1, \+1: (\S+) (\S+)', t)]
ok(len(zs) == 10 and max(abs(z[0]) for z in zs) < 2.4e-10 and zs[0][1:] == (1.562e-03, 5.690e-04), 'magic function: transform vanishes at ten zeros')
cmpv = {float(a): (float(b), float(c)) for a, b, c in re.findall(r'x=(\S+)\s+(\S+e[-+]\d+)\s+(\S+e[-+]\d+)\n', t)}
ok(all(abs(cmpv[x][0] - v) < 6e-6*max(1, v/1e-3) for x, v in [(1.02, 0.99637), (1.327, 0.46193), (1.502, 0.18325)]) and abs(cmpv[2.504][0] - 1.561e-3) < 1e-6 and all(abs(cmpv[x][0]/cmpv[x][1] - 1) < 0.02 for x in (1.02, 1.1, 1.201, 1.327, 1.502, 1.748, 2.504)), 'magic function against our dual solution')
tr = [float(x) for x in re.findall(r'X=\S+\s+(\S+e-\d+)\n', t.split('truncated magic function')[1] + '\n')]
ok([('%.1e' % x) for x in tr] == ['1.9e-02', '1.1e-02', '7.6e-03', '4.2e-03', '1.8e-03', '6.8e-04', '1.5e-04', '3.6e-05'] and tr == sorted(tr, reverse=True), 'truncation table')
ok(abs(5*exp(-2*pi*4.83e-4) - 4.9848) < 1e-4 and abs(5*exp(-2*pi*1.64e-5) - 4.9995) < 1e-4 and abs(3*exp(-2*pi*9.99e-4) - 2.9812) < 1e-4 and abs(3*exp(-2*pi*1.50e-4) - 2.9972) < 1e-4, 'Johns Prop 5.11 as conductors')
t = rd('riemann_kernel_cost_out.txt')
ok('W(g)/g(0) = 0.0751971327' in t and 'g(0) = 0.319751812' in t and '= 16.16596749   quadrature: 16.16596749   12 c1 = 193.992' in t and abs(exp(-0.0751971327) - 0.9276) < 5e-5 and abs(2*pi*mpf('1.4291572e-1060')/mpf('8.98e-1060') - 1) < mpf('1e-3'), 'Riemann kernel cost, truncation constant')
t = rd('cert_check_out.txt')
ok('B = 1: g(0) = 0.319751811961621  pole = 0.5' in t and t.count('W/g(0) = 0.075197132697534236812') == 2, 'certificates: B = 1, zero side and prime side')
rows = [l.split('|') for l in t.split('class A certificates')[1].split('class E optima')[0].splitlines() if re.match(r'\s*\d+ \|', l)]
bounds = {4: '6.19402005455e-3', 8: '3.70719233613e-4', 12: '4.64624013887e-5', 16: '3.58920272265e-6', 20: '4.43833859341e-8', 24: '1.70453120117e-9', 28: '1.35254770865e-10', 32: '1.07689720773e-12', 36: '1.12161611521e-16', 40: '6.85355154306e-19', 44: '2.33130821521e-21', 48: '2.37581468402e-26'}
ok(len(rows) == 12 and all(mpf(r[1]) == mpf(r[2]) and mpf(r[3]) < mpf('1e-80') and int(r[5]) == 0 for r in rows), 'certificates class A: zero side = prime side, B without real root')
ok(all(mpf(r[1]) <= mpf(bounds[int(r[0])]) and mpf(bounds[int(r[0])])/mpf(r[1]) - 1 < mpf('1e-10') for r in rows), 'certificates class A: recomputed cost below the certified bound, equal to 10 digits')
import json
balls = {}
for d in range(2, 25, 2):
    c = json.load(open(os.path.join(here, 'theory', 'convex', 'certs', 'certificate_A_d%d.json' % d)))
    m = re.match(r'\[(\S+) \+/- (\S+)\]', c['cost_over_g0']); balls[2*d] = (mpf(m.group(1)), mpf(m.group(2)), int(c['B_roots_in_s>=0']), int(c['intervals_[0,U]']), c['U'])
ok(len(balls) == 12 and all(abs(mpf(r[1])/balls[int(r[0])][0] - 1) < mpf('1e-18') for r in rows) and all(b[0] + b[1] <= mpf(bounds[d]) and mpf(bounds[d])/(b[0] + b[1]) - 1 < mpf('1e-10') and b[2] == 0 for d, b in balls.items()), 'certificates class A: recomputed cost = ball centre; table bound = ball rounded up')
ok(max(b[3] for b in balls.values()) == 3524 and [balls[d][4] for d in range(4, 49, 4)] == ['5']*8 + ['6']*4, 'certificates: interval counts and U')
dc = rd('theory/convex/certs/decay.out')
optA = {int(a): float(b) for a, b in re.findall(r'\n\s*(\d+)\s+\d+ \| (\S+)\s+\|', dc.split('class E')[0])}
optE = {int(a): float(b) for a, b in re.findall(r'\n\s*(\d+)\s+\d+ \| (\S+)\s+\|', dc.split('\nclass E   (deg')[1])}
ok([('%.4e' % optA[d]) for d in range(0, 49, 4)] == ['7.5197e-02', '6.1940e-03', '3.7072e-04', '4.6462e-05', '3.5892e-06', '4.4383e-08', '1.7045e-09', '1.3525e-10', '1.0769e-12', '1.1216e-16', '6.8535e-19', '2.3313e-21', '2.3758e-26'] and all(optA[d] < float(bounds[d]) for d in bounds), 'optimum per degree, class A (below the certified bound)')
rowsE = [l.split('|') for l in t.split('class E optima')[1].split('sampled positivity')[0].splitlines() if re.match(r'\s*\d+ \|', l)]
ok(len(rowsE) == 8 and all(abs(float(r[1])/optE[int(r[0])] - 1) < 1e-12 and mpf(r[1]) == mpf(r[2]) for r in rowsE) and [('%.4e' % optE[d]) for d in range(4, 33, 4)] == ['6.1424e-03', '7.5310e-05', '1.2559e-07', '8.9013e-14', '5.9355e-18', '4.7813e-25', '1.1538e-29', '4.4835e-37'], 'optimum per degree, class E: cost recomputed')
ok('2.072189606' in dc and '2.005724955, 3.083922359' in dc and '2.000066576, 3.001754992, 4.044541305, 5.092273215' in dc and '2, 3, 4, 5, 7, 8, 9, 11, 13, 16.00000048, 17.00000307, 19.00014963' in dc and 'deg 48: n=23 99.9%' in dc and 's = -1.149, a = 1.74, rms 0.688' in dc, 'touch points, cost carrier and decay fit')
pos = re.findall(r'deg\s+(\d+): non-positive values on u = k/100, k <= 380: (\d+); smallest g/g\(0\) = (\S+) at', t)
ok([p[0] for p in pos] == ['8', '24', '48'] and all(p[1] == '0' and float(p[2]) > 0 for p in pos), 'certificates: g > 0 on the sampled grid')
ok('= 2.788432137 ' in t and t.count('4.212556303e-19') == 2 and 'W/(2 pi g(0)) = 1.510725776e-19' in t and 'W/g(0) = 9.492169997e-19' in t and 2.7884318 < 2.788432137 < 2.7884322 and 4.212556303e-19 <= 4.2125564e-19 and 1 - 4.212556303e-19/4.2125564e-19 < 3e-8 and 1.510725776e-19 <= 1.5107258e-19 and optA[40] < 9.492169997e-19, 'Johns Part I Table 1, row J = 10, rebuilt')
ok(abs(2.37581468402e-26 - 2.376e-26) < 1e-29 and 1 - exp(-2.376e-26) < 2.4e-26, 'certified bound as conductor')
print('all %d checks passed' % n)

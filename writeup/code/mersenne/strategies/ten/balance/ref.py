"""Independent reference for prime5.html: every rule, written apart from the page, for the numbers up to X.
usage: ref.py X out.json [mexponly]     (the coin rule is random and left out)"""
import sys, json, math
X = int(sys.argv[1])
MEXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
MEXPONLY = len(sys.argv) > 3 and sys.argv[3] == 'mexponly'
LIM = 100 if MEXPONLY else 2 * X + 10
spf = list(range(LIM + 1))
for i in range(2, int(LIM ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, LIM + 1, i):
            if spf[j] == j:
                spf[j] = i
def isp(n): return n >= 2 and spf[n] == n
def facts(n):
    om = w = 0; sq = True
    while n > 1:
        p = spf[n]; e = 0
        while n % p == 0: n //= p; e += 1
        om += e; w += 1
        if e > 1: sq = False
    return om, w, sq
def dsum(n, b):
    s = 0
    while n: s += n % b; n //= b
    return s
P = [] if MEXPONLY else [n for n in range(2, X + 1) if isp(n)]
def items(rule, src):
    """list of (n, side or None) in placing order; side None = decided by the running state (alt, light)"""
    if src == 'mexp':
        out = []
        for n in MEXP:
            if n > X: break
            if rule in ('mod4', 'mod8') and n == 2: continue
            if rule == 'mod3' and n <= 3: continue
            if rule == 'lastdig' and n in (2, 5): continue
            out.append((n, None))
        return out
    if rule in ('liou', 'mob', 'omega'):
        out = []
        for n in range(1, X + 1):
            om, w, sq = facts(n) if n > 1 else (0, 0, True)
            if rule == 'liou': out.append((n, 1 if om % 2 == 0 else -1))
            elif rule == 'mob':
                if sq: out.append((n, 1 if om % 2 == 0 else -1))
            else: out.append((n, 1 if w % 2 == 0 else -1))
        return out
    out = []
    for i, p in enumerate(P):
        prev = P[i - 1] if i else 0
        if rule in ('alt', 'light', 'binsum', 'decsum'): out.append((p, None))
        elif rule in ('mod4', 'mod8'):
            if p > 2: out.append((p, None))
        elif rule == 'mod3':
            if p > 3: out.append((p, None))
        elif rule == 'lastdig':
            if p not in (2, 5): out.append((p, None))
        elif rule == 'repeat':
            if prev > 5: out.append((p, 1 if p % 10 == prev % 10 else -1))
        elif rule in ('twincous', 'twinsexy', 'sophie'):
            if rule == 'sophie': a, b = isp(2 * p + 1), isp(2 * p - 1)
            else: a, b = isp(p + 2), isp(p + (4 if rule == 'twincous' else 6))
            if a: out.append((p, 1))
            if b: out.append((p, -1))
        elif rule in ('gap24', 'gap26'):
            g = p - prev if prev else 0
            if g == 2: out.append((p, 1))
            elif g == (4 if rule == 'gap24' else 6): out.append((p, -1))
        elif rule == 'updown':
            if i >= 2:
                b = P[i - 1]; a = P[i - 2]
                if b - a != p - b: out.append((b, 1 if p - b > b - a else -1))
    return out
def side_of(rule, n, k, L, R):
    if rule == 'alt': return 1 if k % 2 == 0 else -1
    if rule == 'light': return 1 if L <= R else -1
    if rule == 'mod4': return 1 if n % 4 == 1 else -1
    if rule == 'mod3': return 1 if n % 3 == 1 else -1
    if rule == 'mod8': return 1 if n % 8 == 1 else -1
    if rule == 'lastdig': return 1 if n % 10 in (1, 9) else -1
    if rule == 'binsum': return 1 if dsum(n, 2) % 2 == 0 else -1
    if rule == 'decsum': return 1 if dsum(n, 10) % 2 == 0 else -1
def run(rule, src, wt):
    L = R = 0.0 if wt in ('log', 'recip') else 0
    k = eq = aL = aR = ch = 0; streak = 0; ss = 0; bL = bR = 0; fL = fR = None; eqsteps = []
    for n, s in items(rule, src):
        if s is None: s = side_of(rule, n, k, L, R)
        w = 1 if wt == 'one' else (0.0 if n == 1 else math.log(n)) if wt == 'log' else 1 / n if wt == 'recip' else n
        fac = 3 if rule in ('mod8', 'repeat') and s > 0 else 2 if rule == 'twinsexy' and s > 0 else 1
        w = w * fac
        if s > 0: L += w
        else: R += w
        k += 1
        d = L - R
        if wt in ('log', 'recip'): sg = 0 if abs(d) < 1e-12 * max(1, L + R) else (1 if d > 0 else -1)
        else: sg = 0 if d == 0 else (1 if d > 0 else -1)
        if sg == 0: eq += 1; eqsteps.append(k)
        elif sg > 0: aL += 1
        else: aR += 1
        if sg != 0:
            if sg == ss: streak += 1
            else:
                if ss != 0: ch += 1
                ss = sg; streak = 1
            if sg > 0 and streak > bL: bL = streak
            if sg < 0 and streak > bR: bR = streak
        if sg > 0 and fL is None: fL = [k, n]
        if sg < 0 and fR is None: fR = [k, n]
    return dict(k=k, L=L, R=R, eq=eq, aL=aL, aR=aR, ch=ch, bL=bL, bR=bR, fL=fL, fR=fR, eqsteps=eqsteps[:200])
RULES = ['alt', 'light', 'mod4', 'mod3', 'mod8', 'lastdig', 'repeat', 'binsum', 'decsum', 'twincous', 'twinsexy', 'sophie', 'gap24', 'gap26', 'updown', 'liou', 'mob', 'omega']
MRULES = ['alt', 'light', 'mod4', 'mod3', 'mod8', 'lastdig', 'binsum', 'decsum']
out = {}
for src, rules in ((('mexp', MRULES),) if MEXPONLY else (('primes', RULES), ('mexp', MRULES))):
    for r in rules:
        for wt in ('one', 'val', 'log', 'recip'):
            out['%s|%s|%s' % (src, r, wt)] = run(r, src, wt)
json.dump(out, open(sys.argv[2], 'w'))
print(len(out), 'runs')

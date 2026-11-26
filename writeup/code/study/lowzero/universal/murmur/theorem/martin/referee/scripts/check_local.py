# Lemma F test: q-means at FIXED l of local factors F_v (v not | P) and G_v (v | M), by enumeration of q mod v^J.
import sys; sys.path.insert(0, '/tmp/claude-1000/referee_martin/scripts')
from ref_core import e_v, e_new
def vval(n, v):
    c = 0
    while n % v == 0: n //= v; c += 1
    return c, n
def local(D, v):
    """returns (a, chi) at v for Delta = -D: a = v-adic valuation of conductor f0, chi = (Delta0/v)."""
    Delta = -D
    w, u = vval(Delta, v)
    if v != 2:
        if w % 2 == 1: return (w-1)//2, 0
        return w//2, (1 if pow(u % v, (v-1)//2, v) == 1 else -1)
    if w % 2 == 0:
        if u % 4 == 1: return w//2, (1 if u % 8 == 1 else -1)
        return (w-2)//2, 0
    return (w-3)//2, 0
def F(D, v):
    a, chi = local(D, v); return sum(v**-j for j in range(a)) + v**-a/(1-chi/v)
def G(D, v):
    a, chi = local(D, v); return v**-a*(chi-1)/(1-chi/v)
J = {2: 15, 3: 9, 5: 6, 7: 5}
worst = 0; n = 0
for l in [1201, 3001, 6007, 9001, 7, 11]:
    for Q0 in (1, 3):
        for s in range(0, 7):
            if s == 0 and False: pass
            for v in (2, 3, 5, 7):
                if l % v == 0: continue
                mod = v**J[v]
                tot = 0.0; cnt = 0
                for q in range(1, mod):
                    if q % v == 0: continue
                    D = Q0*q*(4*l - s*s*Q0*q)
                    if D == 0: continue
                    D = abs(D) if D < 0 else D   # only the v-adic class matters; sign handled below
                    Dtrue = Q0*q*(4*l - s*s*Q0*q)
                    # v-adic data depend only on Dtrue mod squares of units; use |.| times sign as integer
                    Dv = Dtrue
                    if v == 3 and Q0 == 1:   # 3 | M: G_3
                        tot += G(Dv, v)
                    elif v == 3 and Q0 == 3:  # 3 | Q0: expect exactly 1
                        tot += F(Dv, v)
                    else:
                        tot += F(Dv, v)
                    cnt += 1
                mean = tot/cnt
                if v == 3 and Q0 == 1: pred = e_new(3, s)
                elif v == 3 and Q0 == 3: pred = 1.0
                else: pred = e_v(v, s)
                dev = abs(mean - pred); worst = max(worst, dev); n += 1
                if dev > 5e-3: print('DEV', l, Q0, s, v, mean, pred)
print('Lemma F / local means over q at fixed l: %d cases, max |mean - closed form| = %.2e' % (n, worst))

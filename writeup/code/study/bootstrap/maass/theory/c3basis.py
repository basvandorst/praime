# Basis data for the optimisation: moments m_j = g_j(0) (j <= 2*DMAX) and g_k(log n) for prime powers n <= NMAX, k <= DMAX.
import sys; sys.dont_write_bytecode = True
# Output: c3basis_D<DMAX>.json (mid-points as decimal strings and radii).  Usage: c3basis.py DMAX NMAX PREC
import sys, json, time
sys.path.insert(0, '/home/bas/praime/writeup/code/study/bootstrap/maass/theory')
from flint import arb, ctx
import c3core as C
DMAX, NMAX, PREC = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ctx.prec = PREC
t0 = time.time()
def is_pp(n):
    for p in range(2, n + 1):
        if n % p == 0:
            while n % p == 0: n //= p
            return n == 1
    return False
pps = [n for n in range(2, NMAX + 1) if is_pp(n)]
mom = C.g_basis(1, 2*DMAX)
print('moments done %.1f s; max rel radius %.2e' % (time.time() - t0, max(float(m.rad()/abs(m.mid())) for m in mom)), flush=True)
vals = {}
worst = 0.0
for n in pps:
    v = C.g_basis(n, DMAX)
    vals[n] = v
    worst = max(worst, max(float(x.rad()) for x in v))
print('values done %.1f s; max abs radius %.2e' % (time.time() - t0, worst), flush=True)
dig = int(PREC*0.30103) + 5
out = {'DMAX': DMAX, 'NMAX': NMAX, 'PREC': PREC,
       'moments': [[m.mid().str(dig, radius=False), float(m.rad())] for m in mom],
       'values': {str(n): [[x.mid().str(dig, radius=False), float(x.rad())] for x in vals[n]] for n in pps}}
json.dump(out, open('/home/bas/praime/writeup/code/study/bootstrap/maass/theory/c3basis_D%d.json' % DMAX, 'w'))
print('written; time %.1f s' % (time.time() - t0))

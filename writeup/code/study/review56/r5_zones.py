"""r5_zones.py (review, own code): the axis zones.  Model (own derivation) for the half zone around the positive real axis
(points a + bi, b >= 1, weight w(b/sqrt a) = exp(-pi u^2/4), e(-|z|), e^{i m arg z}), density rho per lattice point:
  Z = rho [ (2/3)(r1^1.5 - r0^1.5) F_half  +  s i m (r1 - r0) G  -  (m^2/2) 2 (sqrt r1 - sqrt r0) H2  +  i (pi/4) 2 (sqrt r1 - sqrt r0) H4 ]
      - (1/2) C0 (r1 - r0)
  F_half = int_0^inf e^{-pi(1/4+i)u^2} du, G = int_0^inf u e^{..} du, H2 = int u^2 e^{..}, H4 = int u^4 e^{..};
  s = +1 for the zone at angle 0 (e^{i m theta}, theta ~ u/sqrt a), s = -1 for the zone at 90 degrees (mirror: e^{-i m theta});
  the last term: the line b = 0 (Euler-Maclaurin, Sum_{b>=1} = int_0^inf - g(0)/2), C0 = density on the line (lattice 1;
  primes: inert Gaussian primes q, weight ln q^2, density 1 per integer, as in the symmetric Hardy-Littlewood model).
Primes: rho = 4/pi, C0 = 1.  Lattice: rho = 1, C0 = 1."""
import os, json, math, cmath
import numpy as np
SD = os.path.dirname(os.path.abspath(__file__))+'/'
PS = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/resonance/primeside/'
al = math.pi * (0.25 + 1j)
Fh = 0.5 * cmath.sqrt(math.pi / al); G = 1 / (2 * al)
H2 = math.sqrt(math.pi) / 4 * al ** -1.5; H4 = 3 * math.sqrt(math.pi) / 8 * al ** -2.5
# numerical check of the closed forms
u = np.linspace(0, 30, 300001); f = np.exp(-al * u * u)
num = [np.trapezoid(f, u), np.trapezoid(u * f, u), np.trapezoid(u * u * f, u), np.trapezoid(u ** 4 * f, u)]
assert all(abs(a - b) < 1e-8 for a, b in zip(num, [Fh, G, H2, H4])), (num, [Fh, G, H2, H4])

def model(r0, r1, m, s, rho, C0=1.0, terms=('main', 'first', 'second', 'quartic', 'boundary')):
    z = 0j
    if 'main' in terms: z += rho * (2 / 3) * (r1 ** 1.5 - r0 ** 1.5) * Fh
    if 'first' in terms: z += rho * s * 1j * m * (r1 - r0) * G
    if 'second' in terms: z += -rho * (m * m / 2) * 2 * (math.sqrt(r1) - math.sqrt(r0)) * H2
    if 'quartic' in terms: z += rho * 1j * (math.pi / 4) * 2 * (math.sqrt(r1) - math.sqrt(r0)) * H4
    if 'boundary' in terms: z += -0.5 * C0 * (r1 - r0)
    return z

def C(v): return complex(v[0], v[1])
def show(tag, read, mod, sd):
    d = read - mod; print('  %-34s reading %12.0f %+12.0fi  model %12.0f %+12.0fi  dev %8.0f %+8.0fi = %5.2f sd' % (tag, read.real, read.imag, mod.real, mod.imag, d.real, d.imag, abs(d) / sd))

# 1. lattice (theorem check) on (3e7, 1.93263832e8]
A = json.load(open(SD + '/r2_A.json')); r0, r1 = math.sqrt(3e7), math.sqrt(A['X1'])
print('LATTICE, N in (3e7, %.4g], density 1 (no noise; differences are model error)' % A['X1'])
for key, m, s in (('ax0_m0', 0, 1), ('ax0_m4', 4, 1), ('ax90_m4', 4, -1)):
    z = C(A['zones'][key])
    for terms in (('main',), ('main', 'boundary'), ('main', 'first', 'boundary'), ('main', 'first', 'second', 'quartic', 'boundary')):
        mm = model(r0, r1, m, s, 1.0, 1.0, terms); d = z - mm
        print('  %-8s %-40s reading %11.1f %+11.1fi model %11.1f %+11.1fi  diff %8.1f %+8.1fi' % (key, '+'.join(terms), z.real, z.imag, mm.real, mm.imag, d.real, d.imag))

# 2. primes, fresh range: registered model vs corrected model
pred = json.load(open(PS + 'predictions_primeside.json'))
for run, (x0, x1), fn in (('fresh', (3e7, 1e9), 'readings_fresh.json'), ('dry', (6e5, 2e7), 'readings_dry.json')):
    rd = json.load(open(PS + fn)); r0, r1 = math.sqrt(x0), math.sqrt(x1)
    print('PRIMES %s, r in (%.1f, %.1f]' % (run, r0, r1))
    for name, s in (('axis0', 1), ('axis90', -1)):
        z = rd['zones'][name]; sd = z['sd']
        for mm_, key in ((0, 'm0'), (4, 'm4')):
            read = C(z[key]); reg = C(z['pred_' + key])
            pure = model(r0, r1, mm_, s, 4 / math.pi, 1.0, ('main',))
            neigh = reg - pure      # what the registered model adds beyond the bare axis (neighbouring rays)
            corr = model(r0, r1, mm_, s, 4 / math.pi, 1.0) + neigh
            nob = model(r0, r1, mm_, s, 4 / math.pi, 1.0, ('main', 'first', 'second', 'quartic')) + neigh
            print(' %s m=%d  (registered model - bare axis = %.0f %+.0fi)' % (name, mm_, neigh.real, neigh.imag))
            show('registered model', read, reg, sd); show('+ first/second order, no boundary', read, nob, sd); show('+ all terms (corrected)', read, corr, sd)

# 3. primes on (3e7, 1.93263832e8] from my own sums
zs = {'ax0_m0': 0j, 'ax0_m4': 0j, 'ax90_m4': 0j, 'noise': 0.0}
for t in ('_lo',):
    for ch in json.load(open(SD + '/r1_out%s.json' % t))['chunks']:
        for k in zs: zs[k] += (C(ch['zones'][k]) if k != 'noise' else ch['zones'][k])
r0, r1 = math.sqrt(3e7), math.sqrt(193263832); sd = math.sqrt(zs['noise'] / 2)
print('PRIMES own sums, r in (%.1f, %.1f], sd %.0f' % (r0, r1, sd))
for key, m, s in (('ax0_m0', 0, 1), ('ax0_m4', 4, 1), ('ax90_m4', 4, -1)):
    show(key + ' bare axis', zs[key], model(r0, r1, m, s, 4 / math.pi, 1.0, ('main',)), sd)
    show(key + ' corrected', zs[key], model(r0, r1, m, s, 4 / math.pi, 1.0), sd)

"""r6_subreading.py (review, own code): own smooth readings for A4, B6, C23, C47a on (3e7, 193263832] with an own Hann window
on v = sqrt p in [sqrt 3e7, sqrt 193263832], own main term, own prime-power term (lambda-Lambda rule), own c (r3_out.json);
errors: random model and the rms of the 8 control sums; plus the same readings over the full fresh range assembled from my
chunk sums (registered window) for comparison."""
import os, json, math, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/')
import r1_sums as R
SD = R.SD; TP = 2 * math.pi
r3 = json.load(open(SD + '/r3_out.json'))
cs = {k: r3['accel'][k]['c'] for k in ('A4', 'B6', 'C23', 'C47a')}
Qs = {'A4': 4, 'B6': 3, 'C23': 23, 'C47a': 47}
ks = {'A4': 5, 'B6': 7, 'C23': 1, 'C47a': 1}
V0, V1 = math.sqrt(3e7), math.sqrt(193263832)

def K_of(key):
    k = ks[key]; Q = Qs[key]; ph = (2 * k - 1) * math.pi / 4
    return cs[key] * complex(math.cos(ph), -math.sin(ph)) * (2 * math.sqrt(2) / 3) * Q ** -0.25
def MT(key, v0, v1, n=200001):
    k = ks[key]; Q = Qs[key]; kap = 1 / 16 - (k - 1) ** 2 / 4
    v = np.linspace(v0, v1, n); f = R.np.where(True, 1, 1) * np.sin(np.pi * (v - v0) / (v1 - v0)) ** 2 * (1 - 1j * kap / (TP * v / math.sqrt(Q)))
    h = v[1] - v[0]; return 1.5 * K_of(key) * h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())

# small primes and lambda(p) for the prime-power term (p <= sqrt(193263832) = 13902)
PM = 14000
fl = R.prime_flags(0, PM + 1); pr = 1 + np.flatnonzero(fl); pr = pr[pr >= 2]
lam = {k: {} for k in Qs}
p, a, b = R.gauss_primes(0, PM + 1, fl)
for pp, l in zip(p, 2 * np.real(((a + 1j * b) / np.sqrt(p)) ** 4)): lam['A4'][int(pp)] = l
p, a, b = R.eis_primes(0, PM + 1, fl)
for pp, l in zip(p, 2 * np.real(((a + 1j * math.sqrt(3) * b) / np.sqrt(p)) ** 6)): lam['B6'][int(pp)] = l
odd = pr[pr > 2]; lg = R.leg(odd, 23); sp = R.frob_split(odd)
for pp, g, s in zip(odd, lg, sp): lam['C23'][int(pp)] = (2.0 if s else -1.0) if g == 1 else 0.0
l47 = R.c47_lam(odd)
for pp, l in zip(odd, l47): lam['C47a'][int(pp)] = float(l)
lam['A4'][2] = -1.0; lam['B6'][3] = -1.0; lam['C23'][2] = -1.0; lam['C47a'][2] = R.C5a; lam['C23'][23] = 1.0; lam['C47a'][47] = 1.0
def chi(key, p):
    if key == 'A4': return 0 if p == 2 else (1 if p % 4 == 1 else -1)
    if key == 'B6': return 0 if p == 3 else (1 if p % 3 == 1 else -1)
    q = Qs[key]
    if p == q: return 0
    if p == 2: return 1
    return 1 if pow(p, (q - 1) // 2, q) == 1 else -1
def PP(key, x0, x1, v0, v1):
    Q = Qs[key]; tot = 0j
    for p in pr:
        p = int(p); l1 = lam[key].get(p, 0.0); ch = chi(key, p); L = [1.0, l1]
        k = 1
        while p ** (k + 1) <= x1:
            k += 1
            L.append(l1 * L[-1] - ch * L[-2] if ch != 0 else l1 ** k)
            n = p ** k
            if n <= x0: continue
            v = math.sqrt(n); W = math.sin(math.pi * (v - v0) / (v1 - v0)) ** 2 if v0 < v < v1 else 0.0
            tot += L[k] * math.log(p) * n ** -0.25 * W * complex(math.cos(TP * 2 * math.sqrt(n / Q)), -math.sin(TP * 2 * math.sqrt(n / Q)))
    return tot

ch = json.load(open(SD + '/r1_out_lo.json'))['chunks']
print('own window [sqrt 3e7, sqrt 193263832], own MT, own PP, own c:')
for key in Qs:
    S = sum(complex(*c[key]['mine']) for c in ch); nz = sum(c[key]['noise_mine'] for c in ch)
    ctrl = np.array([sum(complex(*c[key]['ctrl_mine'][j]) for c in ch) for j in range(8)])
    M = MT(key, V0, V1); P = PP(key, 3e7, 193263832, V0, V1)
    rho = (S + P) / M; e_rand = math.sqrt(nz / 2) / abs(M); e_ctrl = math.sqrt(np.mean(np.abs(ctrl) ** 2) / 2) / abs(M)
    print('  %-5s c = %.6f  reading %.4f %+.4fi  err random %.4f controls %.4f   (PP/MT %.4f; constant-1 reading %.3f; lattice reading %.3f)' % (
        key, cs[key], rho.real, rho.imag, e_rand, e_ctrl, abs(P / M), 1 / cs[key],
        {'A4': 4 / math.pi, 'B6': 3 * math.sqrt(3) / math.pi, 'C23': math.sqrt(23) / (3 * math.pi), 'C47a': math.sqrt(47) / (5 * math.pi)}[key] / cs[key]))

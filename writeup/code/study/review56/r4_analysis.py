"""r4_analysis.py (review, own code): errors (sandwich over chunks, delete-group jackknife), alternatives' chi^2,
sieved-lattice readings, and the Pythagorean pooled amplitude with mirror pairs counted once."""
import os, json, math, cmath
import numpy as np
SD = os.path.dirname(os.path.abspath(__file__))+'/'
PS = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/resonance/primeside/'
rd = json.load(open(PS + 'readings_fresh.json')); pred = json.load(open(PS + 'predictions_primeside.json'))
const = json.load(open(PS + 'constants.json')); r3 = json.load(open(SD + '/r3_out.json'))
recs = [json.loads(l) for l in open(PS + 'sums_fresh.jsonl')]
def C(v): return complex(v[0], v[1])
TP = 2 * math.pi
keys = list(rd['series'].keys())

def hann(v, v0, v1):
    u = (v - v0) / (v1 - v0); return np.where((u > 0) & (u < 1), np.sin(np.pi * u) ** 2, 0.0)
def mt_chunk(ent, x0, x1, W, n=4001):
    """own main term of the smooth sum over p in (x0, x1]: ((d+1)/2) K int_{x0^(1/d)}^{x1^(1/d)} W(v) (1 - i kappa/t) dv,
    t = 2 pi v / Q^(1/d).  (Derivation: dP = ((d+1)/2d) K x^((1-d)/2d) dx, weight x^(-1/2+1/2d) W, v = x^(1/d).)"""
    d, Q, K, kap = ent['d'], ent['Q'], C(ent['K']), ent['kappa']
    v = np.linspace(x0 ** (1 / d), x1 ** (1 / d), n); f = hann(v, *W) * (1 - 1j * kap / (TP * v / Q ** (1 / d)))
    h = v[1] - v[0]; return (d + 1) / 2 * K * h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())

# ---------------------------------------------------------------- 1. per-series errors
print('1. errors per series (relative to |MT|): registered (controls) / random model / sandwich over chunks / delete-8-groups jackknife')
G = 8; info = {}
for k in keys:
    ent = pred['series'][k]; fr = ent['runs']['fresh']; x0, x1 = fr['read']; W = fr['windows'][0]
    ch = sorted([(r['L'], r['U'], C(r['res'][k]['smooth'][0][0]), r['res'][k]['noise_smooth'][0]) for r in recs if k in r['res'] and r['L'] >= x0 and r['U'] <= x1])
    S = np.array([c[2] for c in ch]); MTc = np.array([mt_chunk(ent, c[0], c[1], W) for c in ch]); nz = np.array([c[3] for c in ch])
    MT = C(fr['MT_smooth'][0]); PP = C(fr['PP_smooth'][0])
    assert abs(MTc.sum() - MT) / abs(MT) < 1e-5, (k, MTc.sum(), MT)
    rho = (S.sum() + PP) / MT
    assert abs(rho - C(rd['series'][k]['smooth'])) < 1e-9
    res = S + PP * MTc / MT - rho * MTc; n = len(S)
    sand = math.sqrt(n / (n - 1) * np.sum(np.abs(res) ** 2) / 2) / abs(MT)
    # delete-one-group jackknife with G contiguous groups of equal noise share
    cum = np.cumsum(nz) / nz.sum(); grp = np.minimum((cum * G - 1e-12).astype(int), G - 1)
    rj = []
    for g in range(G):
        m = grp != g; rj.append((S[m].sum() + PP * MTc[m].sum() / MT) / MTc[m].sum())
    rj = np.array(rj); jk = math.sqrt((G - 1) / G * np.sum(np.abs(rj - rj.mean()) ** 2) / 2)
    info[k] = dict(rho=rho, reg=rd['series'][k]['smooth_err'], rand=rd['series'][k]['smooth_err_random_model'], sand=sand, jk=jk, n=n, grp=grp, S=S, MTc=MTc, MT=MT, PP=PP)
    print('  %-5s n=%2d  reading %.4f %+.4fi   reg %.4f  rand %.4f  sandwich %.4f  jackknife8 %.4f   (reg/rand %.2f, sand/rand %.2f)' % (
        k, n, rho.real, rho.imag, info[k]['reg'], info[k]['rand'], sand, jk, info[k]['reg'] / info[k]['rand'], sand / info[k]['rand']))
rr = np.array([info[k]['reg'] / info[k]['rand'] for k in keys]); sr = np.array([info[k]['sand'] / info[k]['rand'] for k in keys])
print('  mean ratio registered/random %.3f (sd of mean %.3f); sandwich/random %.3f (%.3f)' % (rr.mean(), rr.std(ddof=1) / math.sqrt(len(rr)), sr.mean(), sr.std(ddof=1) / math.sqrt(len(sr))))

def pool(errkey, comp='re'):
    x = np.array([info[k]['rho'].real if comp == 're' else info[k]['rho'].imag for k in keys]); s = np.array([info[k][errkey] for k in keys])
    w = 1 / s ** 2; m = np.sum(w * x) / np.sum(w); e = 1 / math.sqrt(np.sum(w)); chi = np.sum(((x - m) / s) ** 2)
    return m, e, chi
for ek in ('reg', 'rand', 'sand', 'jk'):
    m, e, chi = pool(ek); mi, ei, chii = pool(ek, 'im')
    chi_law = sum(abs(info[k]['rho'] - 1) ** 2 / info[k][ek] ** 2 for k in keys)
    print('  pooled with %-5s errors: slope %.4f +- %.4f  chi2(Re about mean) %.1f/14   Im mean %+.4f +- %.4f chi2(Im) %.1f/14   chi2(law, 30) %.1f' % (ek, m, e, chi, mi, ei, chii, chi_law))
# pooled delete-group jackknife (same groups removed from every series at once)
pj = []
w = {k: 1 / info[k]['reg'] ** 2 for k in keys}
for g in range(G):
    vals = []
    for k in keys:
        I = info[k]; m = I['grp'] != g
        vals.append(((I['S'][m].sum() + I['PP'] * I['MTc'][m].sum() / I['MT']) / I['MTc'][m].sum()).real)
    pj.append(sum(w[k] * v for k, v in zip(keys, vals)) / sum(w.values()))
pj = np.array(pj); print('  pooled slope, delete-1-of-8-groups jackknife (registered weights): err %.4f (registered 0.0055)' % math.sqrt((G - 1) / G * np.sum((pj - pj.mean()) ** 2)))

# ---------------------------------------------------------------- 2. alternatives
print('\n2. alternatives: chi^2 over the 30 smooth readings (registered errors | sandwich errors)')
L1chi = {'A': math.pi / 4, 'B': math.pi / (3 * math.sqrt(3)), 'C23': 3 * math.pi / math.sqrt(23), 'C31': 3 * math.pi / math.sqrt(31), 'C47': 5 * math.pi / math.sqrt(47)}
def l1chi(k): return L1chi[k[0]] if k[0] in 'AB' else (L1chi[k[:3]] if k[0] == 'C' else None)
def phase_flip(k):
    e = pred['series'][k]; ph = e['phi0'] + math.pi / 4; return cmath.exp(2j * ph)
c = {k: const[k]['c_exact'] for k in keys}
alts = {
    'law': {k: 1.0 for k in keys},
    'constant 1 (registered alt.)': {k: 1 / c[k] for k in keys},
    'no effect': {k: 0.0 for k in keys},
    'lattice baseline = 1/L(1,chi_D) only (CM/dihedral; constant 1 for D,E,DD)': {k: ((1 / l1chi(k)) / c[k] if l1chi(k) else 1 / c[k]) for k in keys},
    'lattice baseline (CM/dihedral; law for D,E,DD)': {k: ((1 / l1chi(k)) / c[k] if l1chi(k) else 1.0) for k in keys},
    '1/L(1,Ad) (B dropped)': {k: 1 / const[k]['B'] if 'B' in const[k] else None for k in keys},
    'direct Euler product p<=100': {k: (r3['direct'][k]['100'] / c[k] if k in r3['direct'] else r3['direct100_nonCM'][k] / c[k]) for k in keys},
    'sieved lattice P=10 (ray form; law for D,E,DD)': {k: (r3['sieved'][k]['P=10'] / c[k] if k in r3['sieved'] else 1.0) for k in keys},
    'sieved lattice P=100 (ray form; law for D,E,DD)': {k: (r3['sieved'][k]['P=100'] / c[k] if k in r3['sieved'] else 1.0) for k in keys},
    'sieved lattice P=1000 (ray form; law for D,E,DD)': {k: (r3['sieved'][k]['P=1000'] / c[k] if k in r3['sieved'] else 1.0) for k in keys},
    'phase conjugated: e^{+i(phi0+pi/4)}': {k: phase_flip(k) for k in keys},
    'without root number (D37a -> -1)': {k: (-1.0 if k == 'D37a' else 1.0) for k in keys},
}
for nm, a in alts.items():
    if any(v is None for v in a.values()): a = {k: (v if v is not None else 1.0) for k, v in a.items()}
    chi_r = sum(abs(info[k]['rho'] - a[k]) ** 2 / info[k]['reg'] ** 2 for k in keys)
    chi_s = sum(abs(info[k]['rho'] - a[k]) ** 2 / info[k]['sand'] ** 2 for k in keys)
    cm = [k for k in keys if k[0] in 'ABC']
    chi_cm = sum(abs(info[k]['rho'] - a[k]) ** 2 / info[k]['reg'] ** 2 for k in cm)
    print('  %-75s chi2 %9.1f | %9.1f   (CM/dihedral 22 readings: %8.1f)' % (nm, chi_r, chi_s, chi_cm))
print('  predicted readings under the lattice baseline and sieved lattices (CM/dihedral):')
for k in keys:
    if k[0] in 'ABC':
        print('   %-5s reading %.4f +- %.4f   lattice %.3f  P=10 %.3f  P=100 %.3f  P=1000 %.4f  direct p<=100 %.3f' % (k, info[k]['rho'].real, info[k]['reg'],
              (1 / l1chi(k)) / c[k], r3['sieved'][k]['P=10'] / c[k], r3['sieved'][k]['P=100'] / c[k], r3['sieved'][k]['P=1000'] / c[k], r3['direct'][k]['100'] / c[k]))
# what the data measure: the arithmetic factor c * L(1,chi_D) per CM/dihedral series, measured vs predicted
print('  arithmetic factor beyond the lattice, c L(1,chi_D): predicted vs measured (reading x predicted)')
for k in keys:
    if k[0] in 'ABC':
        f = c[k] * l1chi(k); print('   %-5s predicted %.3f  measured %.3f +- %.3f' % (k, f, f * info[k]['rho'].real, f * info[k]['reg']))

# ---------------------------------------------------------------- 3. Pythagorean pooled amplitude
print('\n3. Pythagorean windows (m = 0): mirror pairs')
zp = pred['zones']['fresh']; Z = rd['zones']
names = rd['zones_pooled_pyth']['windows']
num = den = 0.0; seen = set()
for n in names:
    z = C(Z[n]['m0']); p = C(zp[n]['0']); sd = Z[n]['sd']
    u, v = n.split('+'); mir = v[:-1] + '+' + u + 'i'
    same = abs(z - C(Z[mir]['m0'])) < 1e-6 * abs(z)
    amp = (z * p.conjugate()).real / abs(p) ** 2; e = sd / abs(p)
    print('  %-7s amplitude %.4f +- %.4f (%+.1f sd)   identical to mirror %s: %s' % (n, amp, e, (amp - 1) / e, mir, same))
    key = tuple(sorted([n, mir]))
    if key in seen: continue
    seen.add(key); num += (z * p.conjugate()).real / sd ** 2; den += abs(p) ** 2 / sd ** 2
print('  pooled over the 6 distinct windows: %.4f +- %.4f   (registered, 12 windows: %.4f +- %.4f)' % (num / den, 1 / math.sqrt(den), rd['zones_pooled_pyth']['amplitude'], rd['zones_pooled_pyth']['err']))

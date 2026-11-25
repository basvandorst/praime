"""scan.py: chirp spectrum P(alpha, beta) = |sum_p a(p) ln p W(p^beta) e(-alpha p^beta)|^2 / sum_p (a(p) ln p W)^2 on the
DISCOVERY range p in [1e5, 3e7) (beta = 1: p in [1e5, 1e7)), W = Hann window in u = p^beta.
Labels from labels_3e7.npz (scan_labels.py); controls: rnd (random +-1) and NSCR scrambles (random permutations of each label).
Evaluation on the uniform grid alpha_j = j * dalpha, dalpha = 1/(OS (u1 - u0)) (Rayleigh resolution / OS), alpha <= A_beta,
by binning in u (bin width h, pi A h <= 1) with a Taylor expansion of e(-alpha delta) to order NT (error < 1e-9 relative).
beta = 1: exact real FFT on the integers.  Output: scan_peaks.json (all local maxima with P > PMIN) and scan_max.json."""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import spectro_lib as S

Z = np.load(os.path.join(S.HERE, 'labels_3e7.npz'))
p_all = Z['p']; m = p_all >= 100000; P = p_all[m].astype(np.float64)
LABELS = ['l11', 'l37', 'c23', 'a4', 'rnd', 'prod', 'sym2']
LAB = {k: Z[k][m] for k in LABELS}
LNP = np.log(P)
BETAS = [(1, 5), (1, 4), (1, 3), (2, 5), (1, 2), (2, 3), (3, 4), (1, 1)]
AMAX = {(1, 5): 5.0, (1, 4): 5.0, (1, 3): 5.0, (2, 5): 5.0, (1, 2): 5.0, (2, 3): 1.2, (3, 4): 0.6, (1, 1): 0.5}
OS = 4; NT = 12; PMIN = 12.0; NSCR = 3
rng = np.random.default_rng(777)

def spectrum(u, w, u0, u1, A):
    """S(alpha_j), alpha_j = j * dalpha, j = 0..J (complex), by binning + Taylor; returns alphas, S."""
    dal = 1.0 / (OS * (u1 - u0)); h = 1.0 / (math.pi * A)
    M = int(math.ceil(1.0 / (dal * h))); M = 1 << int(math.ceil(math.log2(M)))
    h = 1.0 / (M * dal)                                   # exact relation M dal h = 1
    J = int(A / dal)
    k = np.rint((u - u0) / h).astype(np.int64); d = (u - u0) - k * h
    assert k.max() < M
    out = np.zeros(J + 1, complex); al = np.arange(J + 1) * dal
    dpow = np.ones_like(d)
    for t in range(NT + 1):
        g = np.bincount(k, weights=w * dpow, minlength=M)
        F = np.fft.fft(g)[:J + 1]                          # sum_k g_k e(-j k / M)
        out += ((-2j * math.pi * al) ** t / math.factorial(t)) * F
        dpow = dpow * d
    out *= np.exp(-2j * math.pi * al * u0)
    return al, out

def spectrum_beta1(pp, w, A):
    n0 = int(pp.min()); n1 = int(pp.max()) + 1
    g = np.zeros(n1 - n0); np.add.at(g, (pp - n0).astype(np.int64), w)
    F = np.fft.rfft(g); al = np.arange(len(F)) / len(g)
    k = al <= A
    return al[k], F[k] * np.exp(-2j * math.pi * al[k] * n0)

def peaks(al, Pw, thr):
    i = np.nonzero((Pw[1:-1] > Pw[:-2]) & (Pw[1:-1] >= Pw[2:]) & (Pw[1:-1] > thr))[0] + 1
    return i

t0 = time.time(); allpk = []; maxes = {}
for (a, b) in BETAS:
    beta = a / b; A = AMAX[(a, b)]
    if (a, b) == (1, 1):
        mk = P < 1e7; pp = P[mk]; u = pp; u0, u1 = 1e5, 1e7
    else:
        mk = slice(None); pp = P; u = P ** beta; u0, u1 = 1e5 ** beta, 3e7 ** beta
    W = S.hann(u, u0, u1) * (LNP[mk] if (a, b) != (1, 1) else np.log(pp))
    for lab in LABELS + ['scr_%s_%d' % (l, j) for l in ('l11', 'c23', 'a4') for j in range(NSCR)]:
        if lab.startswith('scr_'):
            base = lab.split('_')[1]; lv = LAB[base][mk].copy(); rng.shuffle(lv)
        else:
            lv = LAB[lab][mk]
        w = lv * W; nz = float(np.dot(w, w))
        if (a, b) == (1, 1): al, Sv = spectrum_beta1(pp, w, A)
        else: al, Sv = spectrum(u, w, u0, u1, A)
        Pw = np.abs(Sv) ** 2 / nz
        lo = int(0.5 * OS) + 1                               # skip alpha ~ 0 (the DC lobe: label mean)
        Pw[:lo] = 0.0
        ii = peaks(al, Pw, PMIN)
        maxes['%s|%d/%d' % (lab, a, b)] = dict(max=float(Pw.max()), at=float(al[int(np.argmax(Pw))]), ngrid=len(Pw), dal=float(al[1] - al[0]))
        for i in ii:
            allpk.append(dict(label=lab, beta='%d/%d' % (a, b), alpha=float(al[i]), P=float(Pw[i]), S=[float(Sv[i].real), float(Sv[i].imag)], noise=nz))
        print('%-10s beta=%d/%d  grid %d (dalpha %.3g)  max P %.1f at alpha %.6f  peaks>%g: %d   %.0fs' % (lab, a, b, len(Pw), al[1] - al[0], Pw.max(), al[int(np.argmax(Pw))], PMIN, len(ii), time.time() - t0), flush=True)
json.dump(allpk, open(os.path.join(S.HERE, 'scan_peaks.json'), 'w'), indent=0)
json.dump(maxes, open(os.path.join(S.HERE, 'scan_max.json'), 'w'), indent=1)
print('done %.0fs' % (time.time() - t0))

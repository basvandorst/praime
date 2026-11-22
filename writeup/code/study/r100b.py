"""Round 100 addendum (after the thirty-fourth review): the deficit of the swing against GUE from the arithmetic
correction to the GUE form factor (Berry 1988; Bogomolny-Keating 1996; diagonal approximation, heuristic). With
L(D) = 1/(4 + D^2) = (1/4) int e^{-2|tau|} e^{i D tau} d tau, the expected Q per zero differs from GUE by the form factor
at tau < T = ln(gamma/2pi): GUE has the ramp tau/T, the zeros have peaks Lambda(n)^2/n delta(tau - ln n)/T, so
  Q_zeros - Q_GUE ~ sum_{gamma>0} a_gamma^2 [sum_{ln n < T} Lambda(n)^2/n^3 - int_0^T tau e^{-2 tau} d tau]/T, a = |c|^2
(times 2 for both signs of gamma, which the factor-2 coherence of the RMS formula already counts once more).
usage: r100b.py"""
import os, json, math
import numpy as np
import sympy

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, '..', '..', 'data', 'study')
Z = np.array([float(l) for l in open(os.path.join(HERE, '..', '..', 'data', 'riemann', 'zeros_10000_flint.txt')) if l.strip()])
pp = []
for p in sympy.primerange(2, 3000):
    v = p
    while v < 3000:
        pp.append((v, math.log(p))); v *= p
pp.sort()
r100 = json.load(open(os.path.join(OUT, 'r100.json')))['summary']
out = {}
for d in (0.01, 0.005, 0.003):
    rho = 0.5 + 1j * Z; a = np.abs(((1 + d) ** rho - 1) / rho) ** 2
    T = np.log(Z / (2 * math.pi))
    prim = np.array([sum(L * L / n ** 3 for n, L in pp if math.log(n) < t) for t in T])
    ramp = 0.25 - (2 * T + 1) * np.exp(-2 * T) / 4                       # int_0^T tau e^{-2 tau} d tau
    dQ_half = float(np.sum(a * a * (prim - ramp) / T))                    # gamma > 0 only
    rms = r100[str(d)]['actual']; Q = (rms * d / 2) ** 2 / 2
    for factor in (1, 2, 4):
        out.setdefault(str(d), {})['dRMS_factor_%d' % factor] = (2 / d) * math.sqrt(2) * factor * dQ_half / (2 * math.sqrt(Q))
    out[str(d)]['observed_vs_r100_sets'] = rms - r100[str(d)]['gue_mean']
    print(d, {k: round(v, 5) for k, v in out[str(d)].items()})
json.dump(out, open(os.path.join(OUT, 'r100b.json'), 'w'), indent=1)

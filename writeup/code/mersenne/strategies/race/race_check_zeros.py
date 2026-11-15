"""Independent recheck of race_theory.py.

(a) L(1, chi_d) and L'(1, chi_d) from the Dirichlet series sum chi(n)/n and -sum chi(n) ln n / n
    (numpy, n <= 2*10^6, partial sums averaged over the last full period), Kronecker symbol via sympy.jacobi_symbol;
    then V(chi_d) = ln(|d|/pi) + psi((1+a)/2) + 2 L'/L with psi(1) = -gamma, psi(1/2) = -gamma - 2 ln 2.
(b) For d = -3, -4, 5, -19: the zeros themselves. Lambda(1/2 + it) = (q/pi)^((s+a)/2) Gamma((s+a)/2) L(s, chi) is real
    for real primitive chi; its sign changes for 0 < t <= 60 are the zeros; V = sum 1/(1/4 + gamma^2) over the found
    zeros (both signs) plus the tail (1/pi)(ln(qT/2pi) + 1)/T from the zero density (1/2pi) ln(qt/2pi).
"""
import json
import math
import numpy as np
import mpmath as mp
import sympy

th = json.load(open("race_" + "theory.json"))
EG = 0.57721566490153286


def kron(d, n):
    if n % 2 == 0:
        if d % 2 == 0:
            return 0
        e = 0
        while n % 2 == 0:
            n //= 2
            e += 1
        s2 = 1 if d % 8 in (1, 7) else -1
        return s2 ** e * (int(sympy.jacobi_symbol(d % n, n)) if n > 1 else 1)
    return int(sympy.jacobi_symbol(d % n, n)) if n > 1 else 1


M = 2_000_000
n = np.arange(1, M + 1, dtype=np.float64)
worst = 0.0
print("   d   L(1) series   L(1) theory   L'(1) series   L'(1) theory   V series   V theory")
for ds, row in th["V"].items():
    d = int(ds)
    q = abs(d)
    per = np.array([kron(d, k) if k > 0 else 0 for k in range(q)], dtype=np.float64)
    chi = per[(np.arange(1, M + 1)) % q]
    s1 = np.cumsum(chi / n)
    s2 = np.cumsum(-chi * np.log(n) / n)
    L = s1[-q:].mean()
    Lp = s2[-q:].mean()
    a = 0 if d > 0 else 1
    psi = -EG if a == 1 else -EG - 2 * math.log(2)
    V = math.log(q / math.pi) + psi + 2 * Lp / L
    worst = max(worst, abs(V - row["V"]))
    print(f"{d:4d}  {L:11.6f}  {row['L1']:11.6f}  {Lp:12.6f}  {row['Lp1']:12.6f}  {V:9.5f}  {row['V']:9.5f}")
print(f"max |V(series) - V(theory)| = {worst:.2e}")

mp.mp.dps = 25
T = 60.0
for d in (-3, -4, 5, -19):
    q = abs(d)
    a = 0 if d > 0 else 1
    chi = [kron(d, k) if k > 0 else 0 for k in range(q)]

    def Lam(t):
        s = mp.mpf(0.5) + 1j * mp.mpf(t)
        v = (mp.mpf(q) / mp.pi) ** ((s + a) / 2) * mp.gamma((s + a) / 2) * mp.dirichlet(s, chi)
        return v

    ts = np.arange(0.05, T + 1e-9, 0.05)
    vals = []
    maximag = 0.0
    for t in ts:
        v = Lam(t)
        maximag = max(maximag, abs(float(v.imag)) / max(abs(float(v.real)), 1e-300))
        vals.append(float(v.real))
    zeros = []
    for i in range(len(ts) - 1):
        if vals[i] == 0 or vals[i] * vals[i + 1] < 0:
            g = mp.findroot(lambda t: Lam(t).real, (ts[i], ts[i + 1]), solver="anderson")
            zeros.append(float(g))
    S = 2 * sum(1 / (0.25 + g * g) for g in zeros)
    tail = (math.log(q * T / (2 * math.pi)) + 1) / (math.pi * T)
    NT = T / (2 * math.pi) * math.log(q * T / (2 * math.pi * math.e))
    vt = [r["V"] for k, r in th["V"].items() if int(k) == d][0]
    print(f"d = {d:3d}: zeros in (0, {T:.0f}]: {len(zeros)} (density formula ~{NT:.1f}); first three {', '.join(f'{g:.4f}' for g in zeros[:3])}; "
          f"sum {S:.5f} + tail {tail:.5f} = {S + tail:.5f}  vs theory {vt:.5f}; max |Im/Re| of Lambda on the line {maximag:.1e}; "
          f"first zero pair contributes {2 / (0.25 + zeros[0] ** 2):.4f}")

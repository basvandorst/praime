"""Synthetic L'-values on B's zero heights (no L' value of B is used): y = ln L + 0.1 + sum_{p<=13, p not | D} [A_p + eps_p] + N(0, 0.45),
so the planted per-prime effects are the LDM model's (mean and variance) and the within-function coefficients are 1, 1, 0."""
import numpy as np, math, os, sys, importlib.util
spec = importlib.util.spec_from_file_location('lt', sys.argv[3]); lt = importlib.util.module_from_spec(spec); spec.loader.exec_module(lt)
z = np.load(sys.argv[1]); out = sys.argv[2]; os.makedirs(out, exist_ok=True); rng = np.random.default_rng(5)
for k in z.files:
    D = int(k[1:]); q = abs(D); g = np.sort(z[k]); b = g[(g >= 100) & (g <= 600)]; L = np.log(q * b / (2 * math.pi)); y = np.log(L) + 0.1
    for p in (2, 3, 5, 7, 11, 13):
        c = lt.kron(D, p)
        if c == 0: continue
        w = c * p ** -0.5 * np.exp(-1j * b * math.log(p)); y = y - np.log(np.abs(1 - w)) - 2 * math.log(p) / L * np.real(w / (1 - w))
    y = y + rng.normal(0, 0.45 ** 0.5, len(b))
    np.savez(os.path.join(out, k + '.npz'), D=D, gamma=b, logabs=y, arg=0 * b, absL=1e-15 * np.exp(y), seconds=0.0)

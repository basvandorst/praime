"""Priority 5: SL(3,Z), X = 8. (1) smallest eigenvalue of the Hermitian form L(F) and of L(Ad F) with the true
coefficients of (13.5966, 4.7647) from gl3_key.json (and with conjugated coefficients); (2) t* there:
one probability measure per prime on the classes (phi1, phi2, phi3 = -phi1-phi2) of SU(3), solved by column
generation (barrier SDP on a finite atom set; new atoms at the exact maximisers of the dual function), with the
dual bound UB = <Z, A> + sum_p max_phi Phi_p(phi) as a certificate.
L(F): Gamma_R(s+i l1) Gamma_R(s+i l2) Gamma_R(s-i(l1+l2)), c(p^m) = tau_m log p, tau_m = sum_j alpha_j^m,
alpha_j roots of x^3 - a_p x^2 + conj(a_p) x - 1 (key convention); conductor e^{-t}.
Ad: Gamma_R(s)^2 prod_d Gamma_R(s +- i d), c(p^m) = (|tau_m|^2 - 1) log p; conductor e^{-2t}."""
import sys, json
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review2")
import numpy as np
from scipy.optimize import minimize
from weil2 import Engine, prime_powers

key = json.load(open("/home/bas/praime/writeup/data/study/maass/gl3_key.json"))
X = 8.0
N = int(sys.argv[1]) if len(sys.argv) > 1 else 24
form = [f for f in key["forms"] if abs(f["lambda"][0] - 13.5966) < 1e-3][0]
l1, l2 = form["lambda"]
ap = {int(p): complex(*v) for p, v in form["A_p1"].items()}


def R(M):
    return np.block([[M.real, -M.imag], [M.imag, M.real]])


def taus(phi1, phi2, m):
    return np.exp(1j * m * phi1) + np.exp(1j * m * phi2) + np.exp(-1j * m * (phi1 + phi2))


eng = Engine(X, N, nu=800)
pps = prime_powers(X)
kF = [1j * l1, 1j * l2, -1j * (l1 + l2)]
ds = [l1 - l2, 2 * l1 + l2, l1 + 2 * l2]
kA = [0, 0] + [s * 1j * d for d in ds for s in (1, -1)]
AF = eng.arch_herm(kF)
AA = eng.arch_sym(kA).astype(complex)
P = {}; Q = {}; S = {}
for p, m in pps:
    C = eng.C(m * np.log(p)); w = p ** (-m / 2) * np.log(p)
    P[(p, m)] = -w * (C.T + C)                 # coefficient of Re tau_m
    Q[(p, m)] = -1j * w * (C.T - C)            # coefficient of Im tau_m
    S[(p, m)] = -w * (C.T + C)                 # coefficient of (|tau_m|^2 - 1) in Ad


def mats(local):
    """local[(p,m)] = (E tau_m, E |tau_m|^2)."""
    MF = AF.copy(); MA = AA.copy()
    for pm in pps:
        tm, t2 = local[pm]
        MF = MF + tm.real * P[pm] + tm.imag * Q[pm]
        MA = MA + (t2 - 1) * S[pm]
    return MF, MA


def true_local(a):
    loc = {}
    for p, m in pps:
        al = np.roots([1, -a[p], np.conj(a[p]), -1])
        tm = np.sum(al ** m)
        loc[(p, m)] = (tm, abs(tm) ** 2)
    return loc


print("form (%.5f, %.5f), X = 8, N = %d" % (l1, l2, N))
for lab, a in (("true", ap), ("conjugated", {p: np.conj(v) for p, v in ap.items()})):
    MF, MA = mats(true_local(a))
    print("  %-10s: lambda_min L(F) = %+.3e   L(Ad F)/2 = %+.3e   (|alpha|-1 max %.1e)" % (
        lab, np.linalg.eigvalsh(MF)[0], np.linalg.eigvalsh(MA)[0] / 2,
        max(abs(abs(r) - 1) for p in a for r in np.roots([1, -a[p], np.conj(a[p]), -1]))), flush=True)
MF, MA = mats(true_local(ap))
t_true = min(np.linalg.eigvalsh(MF)[0], np.linalg.eigvalsh(MA)[0] / 2)

# ---- column generation for t*
primes = sorted({p for p, m in pps})
ms = {p: [m for q, m in pps if q == p] for p in primes}


def atom_cols(p, phi):
    out = {}
    for m in ms[p]:
        tm = taus(phi[0], phi[1], m)
        out[m] = (tm, abs(tm) ** 2)
    return out


def solve_atoms(atoms, gap_end=1e-9):
    """barrier: maximise t, M_F(w) - t I >= 0, M_A(w) - 2t I >= 0, w_p in simplex over atoms[p]."""
    # variables: free weights w_{p,j}, j < n_p - 1 ; last weight = 1 - sum ; t
    idx = []; base = {}
    for p in primes:
        for j in range(len(atoms[p]) - 1):
            idx.append((p, j))
    nv = len(idx) + 1
    # matrices for each atom
    def blockmats(p, phi):
        col = atom_cols(p, phi)
        BF = sum(col[m][0].real * P[(p, m)] + col[m][0].imag * Q[(p, m)] for m in ms[p])
        BA = sum((col[m][1] - 1) * S[(p, m)] for m in ms[p])
        return R(BF), R(BA)
    C0F = R(AF); C0A = R(AA)
    n2 = C0F.shape[0]
    BFs = np.zeros((nv, n2, n2)); BAs = np.zeros((nv, n2, n2))
    for p in primes:
        last = blockmats(p, atoms[p][-1])
        C0F = C0F + last[0]; C0A = C0A + last[1]
        for j in range(len(atoms[p]) - 1):
            bf, ba = blockmats(p, atoms[p][j])
            i = idx.index((p, j))
            BFs[i] = bf - last[0]; BAs[i] = ba - last[1]
    BFs[-1] = -np.eye(n2); BAs[-1] = -2 * np.eye(n2)
    # linear constraints G x + h >= 0: w >= 0 and 1 - sum >= 0
    G = []; h = []
    for i, (p, j) in enumerate(idx):
        g = np.zeros(nv); g[i] = 1; G.append(g); h.append(0.0)
    for p in primes:
        g = np.zeros(nv)
        for i, (q, j) in enumerate(idx):
            if q == p:
                g[i] = -1
        G.append(g); h.append(1.0)
    G = np.array(G); h = np.array(h)
    blocks = [(C0F, BFs), (C0A, BAs)]
    x = np.zeros(nv)
    for i, (p, j) in enumerate(idx):
        x[i] = 1.0 / len(atoms[p])
    Fx = lambda b, x: b[0] + np.tensordot(x, b[1], axes=(0, 0))
    x[-1] = min(np.linalg.eigvalsh(Fx((b[0], b[1][:-1]), x[:-1]))[0] / k for b, k in zip(blocks, (1, 2))) - 1
    ntot = 2 * n2 + len(h); mu = 1.0; good = (x.copy(), mu)
    while True:
        ok = False
        for it in range(60):
            g = np.zeros(nv); g[-1] = -1 / mu; H = np.zeros((nv, nv))
            for b in blocks:
                Fi = np.linalg.inv(Fx(b, x)); Fi = 0.5 * (Fi + Fi.T)
                W = np.matmul(Fi[None], b[1])                     # F^-1 A_i
                g -= np.einsum("iaa->i", W)
                nvv = W.shape[0]
                H += W.reshape(nvv, -1) @ np.transpose(W, (0, 2, 1)).reshape(nvv, -1).T
            s = G @ x + h
            g -= G.T @ (1 / s); H += (G / s[:, None] ** 2).T @ G
            Dg = np.sqrt(np.diag(H))
            if not np.all(np.isfinite(H)) or np.min(Dg) <= 0:
                print("   bad Hessian: min diag", np.min(np.diag(H)), "finite", np.all(np.isfinite(H)), "min s", np.min(s), flush=True)
                break
            Hs = H / np.outer(Dg, Dg)
            try:
                d = -np.linalg.solve(Hs, g / Dg) / Dg
            except np.linalg.LinAlgError:
                d = -np.linalg.lstsq(Hs, g / Dg, rcond=1e-14)[0] / Dg
            lam = np.sqrt(max(-g @ d, 0))
            if lam < 1e-7:
                ok = True; break
            step = 1.0 if lam < 0.2 else 1 / (1 + lam)
            xn = x + step * d
            try:
                for b in blocks:
                    np.linalg.cholesky(Fx(b, xn))
                assert np.all(G @ xn + h > 0)
            except (np.linalg.LinAlgError, AssertionError):
                break
            x = xn
        if not ok:
            break
        good = (x.copy(), mu)
        if mu * ntot < gap_end:
            break
        mu /= 4
    xg, mug = good
    ZF = mug * np.linalg.inv(Fx(blocks[0], xg)); ZA = mug * np.linalg.inv(Fx(blocks[1], xg))
    sc = np.trace(ZF) + 2 * np.trace(ZA); ZF /= sc; ZA /= sc
    # weights
    W = {}
    for p in primes:
        ws = [x[idx.index((p, j))] for j in range(len(atoms[p]) - 1)]
        W[p] = ws + [1 - sum(ws)]
    # dual function pieces
    coef = {}
    for p in primes:
        for m in ms[p]:
            coef[(p, m)] = (np.sum(ZF * R(P[(p, m)])), np.sum(ZF * R(Q[(p, m)])), np.sum(ZA * R(S[(p, m)])))
    base = np.sum(ZF * R(AF)) + np.sum(ZA * R(AA))
    return x[-1], W, coef, base


def Phi(p, coef, phi1, phi2):
    v = 0.0
    for m in ms[p]:
        a, b, c = coef[(p, m)]
        tm = taus(phi1, phi2, m)
        v = v + a * tm.real + b * tm.imag + c * (np.abs(tm) ** 2 - 1)
    return v


g1 = np.linspace(0, 2 * np.pi, 241)[:-1]
PH1, PH2 = np.meshgrid(g1, g1)
atoms = {p: [(a, b) for a in np.linspace(0, 2 * np.pi, 7)[:-1] for b in np.linspace(0, 2 * np.pi, 7)[:-1]] for p in primes}
for rd in range(14):
    t, W, coef, base = solve_atoms(atoms)
    UB = base; new = 0
    for p in primes:
        V = Phi(p, coef, PH1, PH2)
        i = np.unravel_index(np.argmax(V), V.shape)
        res = minimize(lambda z: -Phi(p, coef, z[0], z[1]), [PH1[i], PH2[i]], method="Nelder-Mead",
                       options=dict(xatol=1e-12, fatol=1e-15, maxiter=4000))
        UB += max(-res.fun, V[i])
        atoms[p].append(tuple(res.x)); new += 1
        # also add the next best grid maxima far from the first
        flat = np.argsort(V.ravel())[::-1][:400:40]
        for f in flat[1:4]:
            atoms[p].append((PH1.ravel()[f], PH2.ravel()[f]))
    a2 = sum(w * taus(*phi, 1) for w, phi in zip(W[2], atoms[2][:len(W[2])]))
    print("round %d: t* in [%+.6f, %+.6f]  (truth %+.3e); mean a_2 = %.4f%+.4fi (true %.4f%+.4fi); atoms %d"
          % (rd, t, UB, t_true, a2.real, a2.imag, ap[2].real, ap[2].imag, sum(len(v) for v in atoms.values())), flush=True)
    if UB - t < 2e-6:
        break

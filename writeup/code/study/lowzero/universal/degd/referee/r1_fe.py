"""Referee check 1: functional equation L(s) = eps X(s) Ltilde(1-s) with
X(s) = Q^{1/2-s} prod Gamma_R(1-s+conj mu_j)/Gamma_R(s+mu_j), independent of the note's code.
Cases: complex Dirichlet character mod 7 (order 6, odd) and mod 7 (order 3, even);
the shifted L_y(s) = L(s+iy, chi) (complex mu = delta + iy, root number eps*Q^{-iy});
Sym^2 of 11a1 (degree 3, PARI lfunsympow)."""
import cypari2, mpmath as mp
P = cypari2.Pari(); P.allocatemem(400 * 10**6); P.set_real_precision(38)
mp.mp.dps = 30

def GR(z): return mp.power(mp.pi, -z / 2) * mp.gamma(z / 2)
def X(s, Q, mu, conjmu=True):
    v = mp.power(Q, mp.mpf(1) / 2 - s)
    for m in mu:
        mb = mp.conj(m) if conjmu else m
        v *= GR(1 - s + mb) / GR(s + m)
    return v
def c(x): return mp.mpc(complex(x))
def lf(L, s): return c(P.lfun(L, P(f"{mp.re(s)}+I*({mp.im(s)})")))

pts = [mp.mpc(0.3, 7.1), mp.mpc(0.81, 23.4), mp.mpc(-0.2, 3.3), mp.mpc(0.5, 41.7), mp.mpc(1.3, -12.2)]
print('PARI', P.version())
for name, chi, chib in (('chi7 order6 (odd)', '[1]', '[5]'), ('chi7 order3 (even)', '[2]', '[4]')):
    G = P('znstar(7,1)')
    L = P.lfuncreate(P(f'[znstar(7,1),{chi}]')); Lb = P.lfuncreate(P(f'[znstar(7,1),{chib}]'))
    prm = P.lfunparams(L); Q = int(prm[0]); k = int(prm[1]); Vga = [int(v) for v in prm[2]]
    eps = c(P.lfunrootres(L)[2])
    mu = [mp.mpf(v) for v in Vga]
    errs, wrong = [], []
    for s in pts:
        lhs = lf(L, s); rhs = eps * X(s, Q, mu) * lf(Lb, 1 - s)
        errs.append(abs(lhs - rhs) / abs(lhs)); wrong.append(abs(lhs - (-eps) * X(s, Q, mu) * lf(Lb, 1 - s)) / abs(lhs))
        # using L itself (not the dual) on the right must fail for a complex character
        wrong.append(abs(lhs - eps * X(s, Q, mu) * lf(L, 1 - s)) / abs(lhs))
    print(f'{name}: Q={Q} k={k} Vga={Vga} eps={complex(eps):.12f} |eps|={float(abs(eps)):.15f}; max rel err {float(max(errs)):.2e}; '
          f'min rel err of wrong variants (-eps, L instead of Ltilde) {float(min(wrong)):.3f}')
    # shifted: L_y(s) = L(s+iy), mu_y = mu + iy, Ltilde_y(w) = Lbar(w - iy), eps_y = eps Q^{-iy}
    if Vga == [1]:
        for y in (2.7, -5.3):
            muy = [mp.mpc(Vga[0], y)]
            epsy = eps * mp.power(Q, mp.mpc(0, -y))
            e1, e2 = [], []
            for s in pts:
                lhs = lf(L, s + mp.mpc(0, y)); lt = lf(Lb, 1 - s - mp.mpc(0, y))
                e1.append(abs(lhs - epsy * X(s, Q, muy) * lt) / abs(lhs))
                e2.append(abs(lhs - epsy * X(s, Q, muy, conjmu=False) * lt) / abs(lhs))
            print(f'   shifted y={y}: mu = {complex(muy[0])}: FE with conj(mu) in X: max rel err {float(max(e1)):.2e}; '
                  f'with mu (no conj): min rel err {float(min(e2)):.3f}')
# Sym^2 11a
E = P.ellinit([0, -1, 1, -10, -20]); L = P.lfunsympow(E, 2)
prm = P.lfunparams(L); Q = int(prm[0]); k = int(prm[1]); Vga = [int(v) for v in prm[2]]
sh = mp.mpf(k - 1) / 2
mu = [mp.mpf(v) + sh for v in Vga]
eps = c(P.lfunrootres(L)[2])
errs, wrong = [], []
for s in pts:
    lhs = lf(L, s + sh); rhs = eps * X(s, Q, mu) * lf(L, 1 - s + sh)
    errs.append(abs(lhs - rhs) / abs(lhs)); wrong.append(abs(lhs + rhs) / abs(lhs))
print(f'Sym2 11a: Q={Q} k={k} Vga={Vga} -> analytic mu={[float(m) for m in mu]} eps={complex(eps)}; max rel err {float(max(errs)):.2e}; flipped eps {float(min(wrong)):.3f}')
print('lfuncheckfeq Sym2 11a:', P.lfuncheckfeq(L))

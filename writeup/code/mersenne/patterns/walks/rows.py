"""Class-row layouts: omega (row = Omega(n)), sieverows (row = smallest prime factor), ruler (row = v2(n)).
Exact row counts against theory, page check of the row index of every n."""
import json, math
import numpy as np
from scipy.special import gamma as Gamma
from common import series, load_pos, N, VIS

S = series(N)
Om, spf, isp, primes = S["Om"][: N + 1], S["spf"][: N + 1], S["isp"], S["primes"]
out = {}
n = np.arange(N + 1)

# ---------------- omega ----------------
K = int(Om.max())
cnt = np.bincount(Om[1:], minlength=K + 1)
js = json.load(open(f"{VIS}/data/omega.json"))
fx, fy = load_pos("omega")
row_page = np.rint((fy[1:] - 2) * K / 1400).astype(int)
out["omega_page_rows_match"] = bool(np.all(row_page == Om[1:]))
out["omega_counts_match_page"] = cnt.tolist() == js["extra"]["rowsCount"]
x = float(N); ll = math.log(math.log(x)); lx = math.log(x)
P = primes[primes < 10 ** 6]
def G(z):
    pr = np.prod((1 - z / P) ** -1 * (1 - 1 / P) ** z)
    return pr / Gamma(z + 1)
C2 = np.prod(1 - 1 / (P[1:].astype(float) - 1) ** 2)
Cn = 0.25 / C2
rows = []
for k in range(K + 1):
    landau = x * ll ** (k - 1) / (math.factorial(k - 1) * lx) if k >= 1 else 1
    ss = G((k - 1) / ll) * landau if 1 <= k <= 6 else None
    y = x / 2 ** k
    nic = Cn * y * math.log(y) if y > 1 else None
    rows.append(dict(k=k, count=int(cnt[k]), landau=round(landau, 1),
                     sathe_selberg=None if ss is None else round(ss, 1),
                     nicolas=None if nic is None else round(nic, 1)))
out["omega_rows"] = rows
out["nicolas_C"] = Cn
# exact identity for large k: N_k(x) = #{odd m : Omega(m) <= k, m * 2^(k - Omega(m)) <= x}
odd = np.arange(1, N + 1, 2)
ident = []
for k in range(K + 1):
    ok = (Om[odd] <= k) & (odd * (2.0 ** (k - Om[odd].astype(float))) <= x)
    ident.append(int(np.sum(ok)))
out["omega_identity_matches"] = ident == cnt.tolist()
out["omega_ratio_next"] = [round(cnt[k + 1] / cnt[k], 4) for k in range(1, K)]

# ---------------- sieverows ----------------
PL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
row = np.full(N + 1, len(PL), dtype=int)
for i, p in enumerate(PL):
    row[spf == p] = i
row[isp & (n > 31)] = len(PL) + 1
row[1] = len(PL)
cnt_s = np.bincount(row[1:], minlength=len(PL) + 2)
js = json.load(open(f"{VIS}/data/sieverows.json"))
fx, fy = load_pos("sieverows")
row_page = np.rint((fy[1:] - 2) * (len(PL) + 1) / 1400).astype(int)
out["sieve_page_rows_match"] = bool(np.all(row_page == row[1:]))
out["sieve_counts_match_page"] = cnt_s.tolist() == js["extra"]["rowsCount"]
srows = []
dens_prev = 1.0
for i, p in enumerate(PL):
    d = dens_prev / p
    dens_prev *= (1 - 1 / p)
    # exact: Legendre phi(N//p, primes < p) by inclusion-exclusion over subsets of smaller primes
    small = PL[:i]
    def phi(m, ps):
        if not ps:
            return m
        return phi(m, ps[1:]) - phi(m // ps[0], ps[1:])
    ex = phi(N // p, small)
    srows.append(dict(p=p, count=int(cnt_s[i]), legendre_phi=ex, density_pred=round(N * d, 1),
                      ratio=round(cnt_s[i] / (N * d), 5)))
rough = phi(N, PL)  # numbers <= N with no prime factor <= 31 (includes 1)
srows.append(dict(row="1 + composites with spf>=37", count=int(cnt_s[len(PL)]),
                  pred_exact=rough - (len(primes) - len(PL)), density_pred=round(N * dens_prev - (len(primes) - 11), 1)))
srows.append(dict(row="primes > 31", count=int(cnt_s[len(PL) + 1]), pi_N_minus_11=len(primes) - 11))
out["sieve_rows"] = srows
out["mertens_product_31"] = dens_prev
out["mertens_3rd_theorem_31"] = math.exp(-0.5772156649) / math.log(31)

# ---------------- ruler ----------------
v = np.zeros(N + 1, dtype=int)
m = n.copy()
for k in range(1, 20):
    v[(n % (2 ** k)) == 0] = k
cnt_r = np.bincount(v[1:])
js = json.load(open(f"{VIS}/data/ruler.json"))
fx, fy = load_pos("ruler")
Kr = int(v[1:].max())
row_page = np.rint((fy[1:] - 2) * Kr / 1400).astype(int)
out["ruler_page_rows_match"] = bool(np.all(row_page == v[1:]))
out["ruler_counts_match_page"] = cnt_r.tolist() == js["extra"]["rowsCount"]
out["ruler_exact_formula"] = all(cnt_r[k] == N // 2 ** k - N // 2 ** (k + 1) for k in range(Kr + 1))
out["ruler_rows"] = [[k, int(cnt_r[k]), N / 2 ** (k + 1)] for k in range(Kr + 1)]
# x-position in row 0 is linear in n: x = 2 + 700 + ((n-1)/2 + 0.5 - cnt0/2) * 1400/cnt0
odds = np.arange(1, N + 1, 2)
xpred = 2 + 700 + ((odds - 1) / 2 + 0.5 - cnt_r[0] / 2) * 1400 / cnt_r[0]
out["ruler_row0_x_linear_maxdev"] = float(np.max(np.abs(fx[odds] - xpred)))

json.dump(out, open("rows_out.json", "w"), indent=1)
for k, val in out.items():
    if isinstance(val, list) and val and isinstance(val[0], dict):
        print(k); [print("  ", r) for r in val]
    else:
        print(k, ":", val)

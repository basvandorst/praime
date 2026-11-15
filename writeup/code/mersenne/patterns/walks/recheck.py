"""Independent recheck in plain Python integers (no numpy), with sympy spot checks.
Recomputes the key numbers of walks_out.json, turtle_out.json, rows_out.json, rect_out.txt,
exponents_out.json and compares."""
import json, math, random
import sympy

N = 500000
spf = [0] * (N + 1)
primes = []
for i in range(2, N + 1):
    if spf[i] == 0:
        spf[i] = i; primes.append(i)
    for p in primes:
        if p > spf[i] or i * p > N:
            break
        spf[i * p] = p
Om = [0] * (N + 1); mu = [0] * (N + 1); mu[1] = 1
for n in range(2, N + 1):
    p = spf[n]; q = n // p
    Om[n] = Om[q] + 1
    mu[n] = 0 if q % p == 0 else -mu[q]
ok = []
# sympy spot checks of the sieve
random.seed(1)
for n in random.sample(range(2, N + 1), 3000):
    ok.append(sympy.primeomega(n) == Om[n] and sympy.mobius(n) == mu[n])
ok.append(sympy.primepi(N) == len(primes))
ok.append(sum(sympy.mobius(n) for n in range(1, 10001)) == -23)

M = [0] * (N + 1); L = [0] * (N + 1); R = [0] * (N + 1)
pset = set(primes)
for n in range(1, N + 1):
    M[n] = M[n - 1] + mu[n]
    L[n] = L[n - 1] + (1 if Om[n] % 2 == 0 else -1)
    R[n] = R[n - 1] + ((1 if n % 4 == 3 else -1) if (n in pset and n > 2) else 0)
W = json.load(open("walks_out.json"))
res = {}
res["M(N)"] = M[N] == W["mertens_end"]
res["L(N)"] = L[N] == W["liouville_end"]
res["R(N)"] = R[N] == W["race4_end"]
res["M zeros <=N"] = sum(1 for n in range(1, N + 1) if M[n] == 0) == W["mertens_stats_N"]["zeros"]
mx = max(range(1000, N + 1), key=lambda n: abs(M[n]) / math.sqrt(n))
res["max|M|/sqrt n, n>=1000"] = (mx, round(abs(M[mx]) / math.sqrt(mx), 6)) == (W["mertens_stats_N"]["argmax"], round(W["mertens_stats_N"]["max_abs"], 6))
res["L zeros"] = [n for n in range(1, N + 1) if L[n] == 0] == [2, 4, 6, 10, 16, 26, 40, 96, 586]
res["L<=0 for 2..N"] = max(L[2:]) == 0
lm = min(range(1000, N + 1), key=lambda n: L[n] / math.sqrt(n))
res["min L/sqrt n"] = lm == W["liouville_stats_N"]["argmax"]
res["R<0 only at 26861,26862"] = [n for n in range(1, N + 1) if R[n] < 0] == [26861, 26862]
res["R=0 count"] = sum(1 for n in range(1, N + 1) if R[n] == 0) == W["race4_stats_N"]["zeros"]
# turtle
x = y = 0; d = 0; DX = [1, 0, -1, 0]; DY = [0, -1, 0, 1]
TX = [0] * (N + 1); TY = [0] * (N + 1)
for n in range(2, N + 1):
    x += DX[d]; y += DY[d]; TX[n] = x; TY[n] = y
    if n in pset:
        d = (d + 1) % 4
T = json.load(open("turtle_out.json"))
res["turtle end"] = [TX[N], TY[N]] == T["end"]
gaps = [primes[i + 1] - primes[i] for i in range(len(primes) - 1)]
rect = sum(1 for k in range(len(gaps) - 3) if gaps[k] == gaps[k + 2] and gaps[k + 1] == gaps[k + 3])
res["rect loops 255 (gaps only)"] = rect == 255
# run lengths = 1 + gaps + tail: real path uses run sequence starting with the 1-step run
runs = [1] + gaps + [N - primes[-1]]
rect_runs = sum(1 for k in range(len(runs) - 3) if runs[k] == runs[k + 2] and runs[k + 1] == runs[k + 3])
res["rect loops on runs == turtle_out"] = rect_runs == T["real"]["rect_loops"]
pts = set(zip(TX[1:], TY[1:]))
res["distinct points"] = len(pts) == T["real"]["distinct_points"]
# rows
Rw = json.load(open("rows_out.json"))
cnt = [0] * 19
for n in range(1, N + 1):
    cnt[Om[n]] += 1
res["omega counts"] = cnt == [r["count"] for r in Rw["omega_rows"]]
PL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
sc = [sum(1 for n in range(2, N + 1) if spf[n] == p) for p in PL]
res["sieve counts"] = sc == [r["count"] for r in Rw["sieve_rows"][:11]]
ru = [0] * 19
for n in range(1, N + 1):
    k = (n & -n).bit_length() - 1; ru[k] += 1
res["ruler counts"] = ru == [r[1] for r in Rw["ruler_rows"]]
# exponents
E = json.load(open("exponents_out.json"))["exponent_values"]
pi = {}
c = 0
for n in range(1, N + 1):
    if n in pset:
        c += 1
    pi[n] = c
res["exponent M,L,R,X,Y,pi mod 4"] = all(
    (M[int(p)], L[int(p)], R[int(p)], TX[int(p)], TY[int(p)], pi[int(p)] % 4) ==
    (v["M"], v["L"], v["R"], v["X"], v["Y"], v["pi_mod4"]) for p, v in E.items())
res["sympy spot checks"] = all(ok)
for k, v in res.items():
    print(f"{k:40s} {v}")
print("ALL OK" if all(res.values()) else "MISMATCH")

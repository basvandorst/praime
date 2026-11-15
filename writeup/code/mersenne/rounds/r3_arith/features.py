"""Exotic arithmetic / combinatorial features of a prime k. Each group(k) -> {name: float array}.
k: sorted int64 array of odd primes (a contiguous run of primes, as the harness and final exam use).
NaN = undefined for that k."""
import os
import sys
import numpy as np
import gmpy2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/bas/praime/writeup/code/mersenne/structure")
LN2 = np.log(2.0)
PHI = (1 + 5 ** 0.5) / 2
_pari = None


def pari():
    global _pari
    if _pari is None:
        import cypari2
        _pari = cypari2.Pari()
        _pari.allocatemem(300 * 10 ** 6)
    return _pari


def gp_load(k):
    g = pari()
    g(f"P=primes([{int(k[0])},{int(k[-1])}])")
    assert int(g("#P")) == len(k), "prime range mismatch"
    return g


def gpvec(expr, ncol=None):
    s = str(pari()(expr)).replace("[", "").replace("]", "").replace(";", ",").replace(" E", "E")
    a = np.array(s.split(","), dtype=np.float64)
    return a if ncol is None else a.reshape(-1, ncol)


class F32Dict(dict):
    """feature dict that stores float32 (keeps a group under the memory cap)"""
    def __setitem__(self, key, val):
        super().__setitem__(key, np.asarray(val, dtype=np.float32))


def popc(x):
    return np.bitwise_count(np.asarray(x, dtype=np.uint64)).astype(np.float64)


def bitlen(x):
    return np.floor(np.log2(np.asarray(x, dtype=np.float64))) + 1


def isprime_arr(a):
    return np.fromiter((gmpy2.is_prime(int(x)) for x in a), dtype=bool, count=len(a)).astype(float)


# ----------------------------------------------------------------------------- representations
def digits_matrix(k, b):
    L = int(np.ceil(np.log(float(k.max()) + 1) / np.log(b))) + 1
    D = np.zeros((len(k), L), dtype=np.int64)
    x = k.copy()
    n = np.zeros(len(k), dtype=np.int64)
    for i in range(L):
        D[:, i] = x % b
        n += (x > 0)
        x //= b
    return D, n          # D[:, i] = digit of b^i; n = number of digits


def g_rep(k):
    out = {}
    N = len(k)
    # Zeckendorf
    F = [1, 2]
    while F[-1] < 4 * 10 ** 8:
        F.append(F[-1] + F[-2])
    rem = k.copy()
    mask = np.zeros(N, dtype=np.uint64)
    for i in range(len(F) - 1, -1, -1):
        t = rem >= F[i]
        rem = np.where(t, rem - F[i], rem)
        mask |= (t.astype(np.uint64) << np.uint64(i))
    top = np.floor(np.log2(mask.astype(np.float64)))
    out["zeck_terms_frac"] = popc(mask) / (top + 1)
    out["zeck_uses_1"] = (mask & np.uint64(1)).astype(float)
    low = np.zeros(N)
    gap = np.zeros(N)
    last = np.full(N, -1.0)
    for i in range(len(F)):
        b = ((mask >> np.uint64(i)) & np.uint64(1)).astype(bool)
        low = np.where(b & (last < 0), i, low)
        gap = np.where(b & (last >= 0), np.maximum(gap, i - last), gap)
        last = np.where(b, i, last)
    out["zeck_lowest"] = np.minimum(low, 12)
    out["zeck_maxgap"] = np.minimum(gap, 12)
    out["zeck_terms_odd"] = popc(mask) % 2
    # Zeckendorf palindrome
    pal = np.ones(N, dtype=bool)
    for i in range(len(F)):
        j = top - i
        valid = i < j
        bi = (mask >> np.uint64(i)) & np.uint64(1)
        bj = (mask >> np.clip(j, 0, 63).astype(np.uint64)) & np.uint64(1)
        pal &= ~valid | (bi == bj)
    out["zeck_palindrome"] = pal.astype(float)
    # factorial base
    x = k.copy()
    s = np.zeros(N); z = np.zeros(N); mx = np.zeros(N); ln = np.zeros(N); smax = np.zeros(N)
    r = 2
    while (x > 0).any():
        d = x % r
        act = x > 0
        s += d; z += (d == 0) & act; mx += (d == r - 1) & act; ln += act; smax += (r - 1) * act
        x //= r
        r += 1
    out["fact_digitsum_frac"] = s / smax
    out["fact_zero_frac"] = z / ln
    out["fact_maxdigit_frac"] = mx / ln
    # balanced ternary
    x = k.copy(); nz = np.zeros(N); neg = np.zeros(N); ssum = np.zeros(N); ln = np.zeros(N)
    while (x > 0).any():
        d = x % 3
        act = x > 0
        dig = np.where(d == 2, -1, d)
        nz += (dig != 0) & act; neg += (dig == -1) & act; ssum += dig * act; ln += act
        x = (x - dig) // 3
    out["bt_nonzero_frac"] = nz / ln
    out["bt_neg_share"] = neg / nz
    out["bt_digitsum_z"] = ssum / np.sqrt(ln)
    # base -2
    x = k.copy(); pc = np.zeros(N); ln = np.zeros(N); runs = np.zeros(N); prev = np.full(N, -1)
    while (x != 0).any():
        act = x != 0
        d = x & 1
        pc += d * act; ln += act
        runs += act & (d != prev)
        prev = np.where(act, d, prev)
        x = -((x - d) // 2)
    out["negabin_pop_frac"] = pc / ln
    out["negabin_runs_frac"] = runs / ln
    out["negabin_len_excess"] = ln - bitlen(k)
    # binary-derived (not weight / runs of k itself)
    L = bitlen(k)
    out["naf_weight_frac"] = popc(k ^ (3 * k)) / L
    out["pop_k2_frac"] = popc(k * k) / bitlen(k * k)
    for m in (3, 5, 7, 11):
        out[f"pop_{m}k_frac"] = popc(m * k) / bitlen(m * k)
    rs = popc(k & (k >> 1))
    out["rudin_shapiro_odd"] = rs % 2
    out["pairs11_frac"] = rs / L
    # Baum-Sweet: no block of zeros of odd length
    x = k.copy(); ok = np.ones(N, dtype=bool); run = np.zeros(N, dtype=np.int64)
    while (x > 0).any():
        act = x > 0
        b = x & 1
        endrun = act & (b == 1)
        ok &= ~(endrun & (run % 2 == 1))
        run = np.where(act & (b == 0), run + 1, 0)
        x >>= 1
    out["baum_sweet"] = ok.astype(float)
    # Stern fusc
    x = k.copy(); a = np.ones(N, dtype=np.int64); bb = np.zeros(N, dtype=np.int64)
    while (x > 0).any():
        odd = (x & 1) == 1
        act = x > 0
        bb = np.where(act & odd, bb + a, bb)
        a = np.where(act & ~odd, a + bb, a)
        x >>= 1
    out["fusc_log_frac"] = np.log(bb) / np.log(k)
    out["fusc_prime"] = isprime_arr(bb)
    # reversals
    rv = np.zeros(N, dtype=np.int64); x = k.copy()
    while (x > 0).any():
        rv = np.where(x > 0, (rv << 1) | (x & 1), rv); x >>= 1
    out["bin_reversal_prime"] = isprime_arr(rv)
    for b in (3, 10):
        D, n = digits_matrix(k, b)
        rv = np.zeros(N, dtype=np.int64)
        for i in range(D.shape[1]):
            rv = np.where(i < n, rv * b + D[:, i], rv)
        out[f"rev{b}_prime"] = isprime_arr(rv)
    # decimal games
    D, n = digits_matrix(k, 10)
    dsum = D.sum(1)
    out["dec_digitsum_prime"] = isprime_arr(dsum)
    x = k.copy(); steps = np.zeros(N)
    for _ in range(20):
        act = x >= 10
        Dx, nx = digits_matrix(np.maximum(x, 1), 10)
        prod = np.where(np.arange(Dx.shape[1])[None, :] < nx[:, None], Dx, 1).prod(1)
        x = np.where(act, prod, x); steps += act
    out["mult_persistence"] = steps
    out["mult_digital_root"] = x.astype(float)
    x = k.copy(); st = np.full(N, np.nan); happy = np.zeros(N)
    for i in range(30):
        Dx, _ = digits_matrix(np.maximum(x, 1), 10)
        x = (Dx ** 2).sum(1)
        hit1 = (x == 1) & np.isnan(st); hit4 = (x == 4) & np.isnan(st)
        happy[hit1] = 1
        st[hit1 | hit4] = i + 1
    out["happy"] = happy
    out["happy_steps"] = st
    # Kaprekar routine with the number of digits of k
    kap = np.full(N, 30.0)
    for lo in range(0, N, 50000):
        x = k[lo:lo + 50000].copy(); nd = n[lo:lo + 50000]
        hist = [x.copy()]
        res = np.full(len(x), 30.0)
        for it in range(1, 30):
            Dx = np.zeros((len(x), int(nd.max())), dtype=np.int64)
            y = x.copy()
            for i in range(Dx.shape[1]):
                Dx[:, i] = y % 10; y //= 10
            Dx = np.where(np.arange(Dx.shape[1])[None, :] < nd[:, None], Dx, -1)
            srt = np.sort(Dx, 1)                      # ascending, pads (-1) first
            pw = 10 ** np.arange(Dx.shape[1])[::-1]
            valid = srt >= 0
            asc = np.zeros(len(x), dtype=np.int64); desc = np.zeros(len(x), dtype=np.int64)
            for i in range(Dx.shape[1]):
                c = srt[:, i]; v = c >= 0
                asc = np.where(v, asc * 10 + c, asc)
            for i in range(Dx.shape[1] - 1, -1, -1):
                c = srt[:, i]; v = c >= 0
                desc = np.where(v, desc * 10 + c, desc)
            x = desc - asc
            seen = np.zeros(len(x), dtype=bool)
            for h in hist:
                seen |= (h == x)
            res = np.where(seen & (res == 30), it, res)
            hist.append(x.copy())
        kap[lo:lo + 50000] = res
    out["kaprekar_steps"] = kap
    return {nm: np.asarray(v, dtype=np.float32) for nm, v in out.items()}


def g_bases(k):
    out = {}
    N = len(k)
    # bases 3..64: normalized digit sum; 3..16: palindrome, longest run, distinct digits
    for b in range(3, 65):
        D, n = digits_matrix(k, b)
        s = D.sum(1)
        mu = n * (b - 1) / 2.0
        sd = np.sqrt(n * (b * b - 1) / 12.0)
        out[f"b{b}_digitsum_z"] = (s - mu) / sd
        if b <= 16:
            idx = np.arange(N)
            pal = np.ones(N, dtype=bool)
            for i in range(D.shape[1]):
                j = n - 1 - i
                valid = i < j
                pal &= ~valid | (D[:, i] == D[idx, np.clip(j, 0, D.shape[1] - 1)])
            out[f"b{b}_palindrome"] = pal.astype(float)
            run = np.ones(N); best = np.ones(N)
            for i in range(1, D.shape[1]):
                same = (D[:, i] == D[:, i - 1]) & (i < n)
                run = np.where(same, run + 1, 1); best = np.maximum(best, run)
            out[f"b{b}_maxrun"] = best
            present = np.zeros((N, b), dtype=bool)
            for i in range(D.shape[1]):
                present[idx, D[:, i]] |= (i < n)
            out[f"b{b}_distinct_frac"] = (present.sum(1) / np.minimum(n, b)).astype(np.float32)
        out[f"b{b}_digitsum_z"] = out[f"b{b}_digitsum_z"].astype(np.float32)
    return out


# ----------------------------------------------------------------------------- Beatty / fractional parts
def g_beatty(k):
    out = {}
    kf = k.astype(np.float64)
    consts = {"invphi": 1 / PHI, "sqrt2": 2 ** 0.5, "sqrt3": 3 ** 0.5, "sqrt5": 5 ** 0.5, "pi": np.pi,
              "e": np.e, "gamma": 0.5772156649015329, "log2_3": np.log2(3), "zeta3": 1.2020569031595942,
              "ln2": np.log(2), "cbrt2": 2 ** (1 / 3), "feigenbaum": 4.669201609102990, "khinchin": 2.685452001065306}
    for nm, a in consts.items():
        # exact split: k*a = k*hi + k*lo to keep ~1e-9 accuracy
        hi = np.float64(np.float32(a)); lo = a - hi
        f = np.mod(np.mod(kf * hi, 1.0) + kf * lo, 1.0)
        out[f"frac_k_{nm}"] = f
        out[f"dist_k_{nm}"] = np.minimum(f, 1 - f)
        if nm in ("invphi", "sqrt2", "pi", "e"):
            with np.errstate(divide="ignore"):
                y = 1 / f
                a1 = np.floor(y)
                a2 = np.floor(1 / (y - a1))
            out[f"cf1_k_{nm}"] = np.minimum(a1, 50)
            out[f"cf2_k_{nm}"] = np.minimum(a2, 50)
    out["frac_sqrt_k"] = np.mod(np.sqrt(kf), 1)
    out["frac_sqrt_2k"] = np.mod(np.sqrt(2 * kf), 1)
    out["frac_cbrt_k"] = np.mod(np.cbrt(kf), 1)
    out["frac_k_pow_quarter"] = np.mod(kf ** 0.25, 1)
    out["frac_k_pow_2_3"] = np.mod(kf ** (2 / 3), 1)
    out["frac_ln_k"] = np.mod(np.log(kf), 1)
    out["frac_lnln_k"] = np.mod(np.log(np.log(kf)), 1)
    out["frac_log10_k"] = np.mod(np.log10(kf), 1)
    out["frac_log3_k"] = np.mod(np.log(kf) / np.log(3), 1)
    out["frac_logphi_k"] = np.mod(np.log(kf) / np.log(PHI), 1)
    out["frac_k_over_lnk"] = np.mod(kf / np.log(kf), 1)
    out["frac_li_ratio"] = np.mod(np.sqrt(kf) * np.log(kf), 1)
    return out


# ----------------------------------------------------------------------------- continued fraction of sqrt(k)
def g_cf(k):
    N = len(k)
    a0 = np.floor(np.sqrt(k.astype(np.float64))).astype(np.int64)
    a0 = np.where((a0 + 1) ** 2 <= k, a0 + 1, a0); a0 = np.where(a0 ** 2 > k, a0 - 1, a0)
    m = np.zeros(N, dtype=np.int64); d = np.ones(N, dtype=np.int64); a = a0.copy()
    L = np.zeros(N); mx = np.zeros(N); slog = np.zeros(N); ones = np.zeros(N)
    act = np.arange(N)
    while len(act):
        mm, dd, aa, kk, a00 = m[act], d[act], a[act], k[act], a0[act]
        mm = dd * aa - mm
        dd = (kk - mm * mm) // dd
        aa = (a00 + mm) // dd
        m[act], d[act], a[act] = mm, dd, aa
        L[act] += 1
        fin = aa == 2 * a00
        inner = ~fin
        mx[act[inner]] = np.maximum(mx[act[inner]], aa[inner])
        slog[act] += np.log(aa)
        ones[act] += aa == 1
        act = act[~fin]
    kf = k.astype(np.float64)
    return {"cfsqrt_period_norm": L / np.sqrt(kf),
            "cfsqrt_period_mod4": L % 4,
            "cfsqrt_maxq_norm": mx / np.sqrt(kf),
            "cfsqrt_meanlogq": slog / L,
            "cfsqrt_ones_frac": ones / L,
            "cfsqrt_reg_norm": slog / np.sqrt(kf)}


# ----------------------------------------------------------------------------- PARI: forms x^2 + n y^2
IDONEAL = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25, 28, 30, 33, 37, 40, 42, 45, 48, 57,
           58, 60, 70, 72, 78, 85, 88, 93, 102, 105, 112, 120}


def g_forms(k):
    g = gp_load(k)
    out = F32Dict()
    for n in range(5, 65):
        if n in IDONEAL:
            continue
        v = gpvec(f"vector(#P,i,my(p=P[i]);if(p<={n}||kronecker(-4*{n},p)!=1,-1,#qfbsolve(Qfb(1,0,{n}),p)>0))")
        out[f"form_x2_{n}y2"] = np.where(v < 0, np.nan, v)
    for n in (1, 2, 3):
        ab = gpvec(f"vector(#P,i,my(p=P[i],s);if(p<=3||kronecker(-4*{n},p)!=1,[0,0],s=qfbsolve(Qfb(1,0,{n}),p);[abs(s[1]),abs(s[2])]))", 2)
        a, b = ab[:, 0], ab[:, 1]
        ok = a > 0
        if n == 1:   # a odd, b even
            a, b = np.where(a % 2 == 1, a, b), np.where(a % 2 == 1, b, a)
        ang = np.arctan2(np.sqrt(n) * b, a) / (np.pi / 2)
        out[f"angle_x2_{n}y2"] = np.where(ok, ang, np.nan)
        for q in (3, 5):
            out[f"x2_{n}y2_{q}_div_x"] = np.where(ok, (a % q == 0).astype(float), np.nan)
            out[f"x2_{n}y2_{q}_div_y"] = np.where(ok, (b % q == 0).astype(float), np.nan)
        if n == 1:
            out["gauss_a_1mod4"] = np.where(ok, (a % 4 == 1).astype(float), np.nan)
            out["gauss_apb_prime"] = np.where(ok, isprime_arr(a + b), np.nan)
            out["gauss_amb_prime"] = np.where(ok, isprime_arr(np.abs(a - b)), np.nan)
    return out


# ----------------------------------------------------------------------------- PARI: orders, primitive roots, Lucas sequences
GP_ORD = r"""
ordm(M,N,p)=my(a=N,f=factor(N)[,1]);for(j=1,#f,while(a%f[j]==0&&M^(a/f[j])==matid(2),a/=f[j]));a;
issc(X)=X[1,2]==0&&X[2,1]==0&&X[1,1]==X[2,2];
ords(M,N,p)=my(a=N,f=factor(N)[,1]);for(j=1,#f,while(a%f[j]==0&&issc(M^(a/f[j])),a/=f[j]));a;
lucasinfo(M,D,p)=my(N=p-kronecker(D,p),a=ords(Mod(M,p),N,p),c=lift((Mod(M,p)^a)[1,1]));[N/a,valuation(a,2)-valuation(N,2),znorder(Mod(c,p))];
"""


def g_orders(k):
    g = gp_load(k)
    for line in GP_ORD.strip().split("\n"):
        g(line.rstrip(";"))
    out = F32Dict()
    kf = k.astype(np.float64)
    big = k > 11
    for b in (2, 3, 5, 6, 7, 10):
        o = gpvec(f"vector(#P,i,my(p=P[i]);if(p<=11,0,znorder(Mod({b},p))))")
        idx = np.where(big, (k - 1) / np.maximum(o, 1), np.nan)
        out[f"idx{b}"] = np.minimum(idx, 64)
        out[f"idx{b}_log_frac"] = np.log(idx) / np.log(kf)
        out[f"idx{b}_is1"] = (idx == 1).astype(float)
        odd = idx.copy()
        okk = ~np.isnan(odd)
        oi = odd[okk].astype(np.int64)
        oi = oi // (oi & -oi)
        tmp = np.full(len(k), np.nan); tmp[okk] = (oi == 1)
        out[f"idx{b}_oddpart_is1"] = tmp
        v2o = np.full(len(k), np.nan)
        oo = o[okk].astype(np.int64); km = (k[okk] - 1)
        v2o[okk] = np.log2((km & -km)) - np.log2((oo & -oo))
        out[f"ord{b}_v2deficit"] = v2o
        for r in (3, 4, 5, 7, 8, 9):
            cond = big & ((k - 1) % r == 0)
            out[f"{b}_is_{r}th_power"] = np.where(cond, (np.nan_to_num(idx) % r == 0).astype(float), np.nan)
    out["least_primroot"] = gpvec("vector(#P,i,my(p=P[i],g=2);while(znorder(Mod(g,p))<p-1,g++);g)")
    out["least_prime_primroot"] = gpvec("vector(#P,i,my(p=P[i],g=2);while(znorder(Mod(g,p))<p-1,g=nextprime(g+1));g)")
    out["least_qnr"] = np.minimum(gpvec("vector(#P,i,my(p=P[i],t=2);while(kronecker(t,p)==1,t++);t)"), 30)
    out["least_noncubic"] = gpvec("vector(#P,i,my(p=P[i],t=2);if(p%3!=1||p<11,0,while(Mod(t,p)^((p-1)/3)==1,t++);t))")
    out["least_noncubic"] = np.where(out["least_noncubic"] == 0, np.nan, np.minimum(out["least_noncubic"], 30))
    out["qr_count_below50"] = gpvec("vector(#P,i,my(p=P[i]);sum(t=1,49,kronecker(t,p)==1))")
    out["legendre_sum_primes50"] = gpvec("vector(#P,i,my(p=P[i]);sum(j=1,15,kronecker(prime(j),p)))")
    out["qr_run_from1"] = np.minimum(gpvec("vector(#P,i,my(p=P[i],t=1);while(kronecker(t+1,p)==1,t++);t)"), 30)
    for nm, M, D in (("fib", "[1,1;1,0]", 5), ("pell", "[2,1;1,0]", 8), ("ll", "[2,3;1,2]", 12)):
        r = gpvec(f"vector(#P,i,my(p=P[i]);if(p<=11,[0,0,0],lucasinfo({M},{D},p)))", 3)
        ok = big
        out[f"{nm}_index"] = np.where(ok, np.minimum(r[:, 0], 64), np.nan)
        out[f"{nm}_index_is1"] = np.where(ok, (r[:, 0] == 1).astype(float), np.nan)
        out[f"{nm}_v2deficit"] = np.where(ok, r[:, 1], np.nan)
        out[f"{nm}_pisano_ratio"] = np.where(ok, r[:, 2], np.nan)
    # Lucas-Lehmer map s -> s^2 - 2 from 4 mod k: cycle length = ord_m(2), m = odd part of ord(2+sqrt3)
    cyc = gpvec("vector(#P,i,my(p=P[i],N,a,m);if(p<=11,1,N=p-kronecker(12,p);a=ordm(Mod([2,3;1,2],p),N,p);m=a>>valuation(a,2);if(m==1,1,znorder(Mod(2,m)))))")
    out["ll_cycle_log_frac"] = np.where(big, np.log(cyc) / np.log(kf), np.nan)
    return out


# ----------------------------------------------------------------------------- PARI: class numbers
def g_classno(k):
    g = gp_load(k)
    out = F32Dict()
    kf = k.astype(np.float64)
    big = k > 11
    for nm, Dexpr in (("h1", "if(p%4==3,-p,-4*p)"), ("h8", "-8*p"), ("h3", "if(p%4==1,-3*p,-12*p)")):
        hv = gpvec(f"vector(#P,i,my(p=P[i]);if(p<=11,0,qfbclassno({Dexpr})))")
        absD = np.where(nm == "h1", np.where(k % 4 == 3, kf, 4 * kf), 0)
        absD = {"h1": np.where(k % 4 == 3, kf, 4 * kf), "h8": 8 * kf, "h3": np.where(k % 4 == 1, 3 * kf, 12 * kf)}[nm]
        out[f"{nm}_L1"] = np.where(big, np.pi * hv / np.sqrt(absD), np.nan)
        hi = hv.astype(np.int64)
        for q in (3, 5, 7):
            out[f"{nm}_div{q}"] = np.where(big, (hi % q == 0).astype(float), np.nan)
        hs = np.maximum(hi, 1)
        out[f"{nm}_v2"] = np.where(big, np.minimum(np.log2(hs & -hs), 6), np.nan)
        out[f"{nm}_oddpart_prime"] = np.where(big, isprime_arr(hs // (hs & -hs)), np.nan)
    return out


# ----------------------------------------------------------------------------- PARI: elliptic curves
CURVES = {"11a": "[0,-1,1,-10,-20]", "37a": "[0,0,1,-1,0]", "389a": "[0,1,1,-2,0]", "5077a": "[0,0,1,-7,6]",
          "14a": "[1,0,1,4,-6]", "15a": "[1,1,1,-10,-10]", "x3px1": "[0,0,0,1,1]", "x3m2": "[0,0,0,0,-2]"}


def g_ell(k):
    g = gp_load(k)
    out = F32Dict()
    kf = k.astype(np.float64)
    big = k > 31
    for nm, c in CURVES.items():
        ap = gpvec(f"E=ellinit({c});vector(#P,i,my(p=P[i]);if(p<=31,0,ellap(E,p)))")
        out[f"ell{nm}_satotate"] = np.where(big, ap / (2 * np.sqrt(kf)), np.nan)
        out[f"ell{nm}_ap0"] = np.where(big, (ap == 0).astype(float), np.nan)
        out[f"ell{nm}_ap_div3"] = np.where(big, (ap % 3 == 0).astype(float), np.nan)
        out[f"ell{nm}_ap_pos"] = np.where(big, (ap > 0).astype(float), np.nan)
        npts = (k + 1 - ap).astype(np.int64)
        out[f"ell{nm}_order_prime"] = np.where(big, isprime_arr(npts), np.nan)
    return out


# ----------------------------------------------------------------------------- PARI: neighbours
NEIGH = {"km1": "p-1", "kp1": "p+1", "2km1": "2*p-1", "2kp1": "2*p+1", "4kp1": "4*p+1", "4km1": "4*p-1",
         "6km1": "6*p-1", "6kp1": "6*p+1", "km2": "p-2", "kp2": "p+2", "kp6": "p+6", "k2p1h": "(p^2+1)/2",
         "k2pkp1": "p^2+p+1", "k2mkp1": "p^2-p+1"}


def g_neigh(k):
    g = gp_load(k)
    g("st(n)=my(f=factor(n),pr=f[,1],m=#pr);if(n==1,return([1,1,0,0,1.,1.,0]));[pr[m],if(m>1,pr[m-1],1),m,vecsum(f[,2]),sigma(f)/n*1.,eulerphi(f)/n*1.,valuation(n,2)]")
    g("isprac(n)=my(f=factor(n),s=1);if(n%2,return(n==1));for(i=1,#f~,if(f[i,1]>s+1,return(0));s*=sigma(f[i,1]^f[i,2]));1")
    out = F32Dict()
    for nm, ex in NEIGH.items():
        r = gpvec(f"vector(#P,i,my(p=P[i]);st({ex}))", 7)
        n = {"km1": k - 1, "kp1": k + 1, "2km1": 2 * k - 1, "2kp1": 2 * k + 1, "4kp1": 4 * k + 1, "4km1": 4 * k - 1,
             "6km1": 6 * k - 1, "6kp1": 6 * k + 1, "km2": k - 2, "kp2": k + 2, "kp6": k + 6,
             "k2p1h": (k * k + 1) // 2, "k2pkp1": k * k + k + 1, "k2mkp1": k * k - k + 1}[nm].astype(np.float64)
        ln = np.log(n)
        out[f"{nm}_lpf_frac"] = np.log(r[:, 0]) / ln
        out[f"{nm}_lpf2_frac"] = np.log(r[:, 1]) / ln
        has2 = (r[:, 6] > 0).astype(float)
        out[f"{nm}_omega_odd"] = np.minimum(r[:, 2] - has2, 6)
        out[f"{nm}_nonsqfree_odd"] = np.minimum(r[:, 3] - r[:, 6] - (r[:, 2] - has2), 4)
        out[f"{nm}_sigma_ratio"] = r[:, 4]
        out[f"{nm}_phi_ratio"] = r[:, 5]
        out[f"{nm}_oddpart_prime"] = ((r[:, 2] - has2 == 1) & (r[:, 3] - r[:, 6] == 1)).astype(float)
        if nm in ("km1", "kp1"):
            out[f"{nm}_v2"] = np.minimum(r[:, 6], 8)
    out["km1_practical"] = gpvec("vector(#P,i,isprac(P[i]-1))")
    out["kp1_practical"] = gpvec("vector(#P,i,isprac(P[i]+1))")
    al = gpvec("vector(#P,i,my(n=P[i]-1,s1=sigma(n)-n,s2=sigma(s1)-s1);[isprime(s1),s2*1./s1])", 2)
    out["aliquot_km1_s1_prime"] = al[:, 0]
    out["aliquot_km1_s2_ratio"] = al[:, 1]
    return out


# ----------------------------------------------------------------------------- dynamics / sequences / gaps
_ALLP = None


def all_primes():
    global _ALLP
    if _ALLP is None:
        path = os.path.join(HERE, "primes_B.npy")
        if os.path.exists(path):
            _ALLP = np.load(path)
        else:
            from common import segmented_primes
            _ALLP = segmented_primes(83880077 + 2000).astype(np.int64)
            np.save(path, _ALLP)
    return _ALLP


def collatz(n):
    n = n.astype(np.int64).copy()
    N = len(n)
    steps = np.zeros(N); odd = np.zeros(N); mx = n.astype(np.float64).copy()
    act = np.flatnonzero(n > 1)
    x = n[act]
    while len(act):
        o = (x & 1) == 1
        x = np.where(o, 3 * x + 1, x >> 1)
        steps[act] += 1; odd[act] += o
        mx[act] = np.maximum(mx[act], x)
        keep = x > 1
        act, x = act[keep], x[keep]
    return steps, odd, mx


def kolakoski(n):
    Kk = np.array([1, 2, 2], dtype=np.int8)
    while len(Kk) < n + 2:
        vals = (1 + (np.arange(len(Kk)) % 2)).astype(np.int8)
        Kk = np.repeat(vals, Kk)
    return Kk[:n + 2]


def g_dyn(k):
    out = {}
    kf = k.astype(np.float64)
    for nm, start in (("k", k), ("2kp1", 2 * k + 1), ("km1", k - 1)):
        s, o, mx = collatz(start)
        ls = np.log(start.astype(np.float64))
        out[f"collatz_{nm}_steps_norm"] = s / ls
        out[f"collatz_{nm}_oddfrac"] = o / s
        out[f"collatz_{nm}_peak"] = np.log(mx) / ls
    Kk = kolakoski(int(k.max()) + 1)          # Kk[0] = K(1)
    out["kolakoski_term"] = Kk[k - 1].astype(float)
    c2k = np.zeros(len(k), dtype=np.int64)       # number of 2s among K(1..k), chunked
    run, CH = 0, 10 ** 7
    for a in range(0, len(Kk), CH):
        cs = np.cumsum(Kk[a:a + CH] == 2, dtype=np.int64) + run
        sel = (k - 1 >= a) & (k - 1 < a + CH)
        c2k[sel] = cs[k[sel] - 1 - a]
        run = int(cs[-1])
    c2 = None
    disc = (2 * c2k - k) / np.sqrt(kf)      # (#2 - #1)/sqrt(k)
    out["kolakoski_disc"] = disc
    del Kk, c2
    G = np.array([0, 1, 2, 2], dtype=np.int64)
    while G.sum() < 2 * k.max():
        m = len(G)
        G = np.concatenate([[0], np.repeat(np.arange(1, m), G[1:m])])[: 4 * m]
    cs = np.cumsum(G)
    Gk = np.searchsorted(cs, k)                 # smallest m with G(1)+..+G(m) >= k  ->  Golomb(k) = m
    out["golomb_ratio"] = Gk / (PHI ** (2 - PHI) * kf ** (PHI - 1))
    ap = all_primes()
    pos = np.searchsorted(ap, k)
    assert (ap[pos] == k).all()
    n = pos + 1                                  # k = p_n
    out["index_prime"] = isprime_arr(n)
    out["index_parity"] = (n % 2).astype(float)
    out["index_mod3"] = (n % 3).astype(float)
    out["index_digitsum_z"] = 0.0 * kf
    D, nd = digits_matrix(n.astype(np.int64), 10)
    out["index_digitsum_z"] = (D.sum(1) - nd * 4.5) / np.sqrt(nd * 8.25)
    out["index_popcount_frac"] = popc(n) / bitlen(n)
    gm = (k - ap[pos - 1]).astype(np.float64); gp = (ap[pos + 1] - k).astype(np.float64)
    out["gap_before_norm"] = gm / np.log(kf)
    out["gap_after_norm"] = gp / np.log(kf)
    out["gap_balanced"] = (gm == gp).astype(float)
    out["gap_before_mod6"] = gm % 6
    out["gap_after_mod6"] = gp % 6
    out["gap_after_gt_before"] = (gp > gm).astype(float)
    gm2 = (ap[pos - 1] - ap[pos - 2]).astype(np.float64)
    out["good_prime"] = (kf * kf > ap[pos - 1].astype(np.float64) * ap[pos + 1]).astype(float)
    out["prev_gap_pair_equal"] = (gm2 == gm).astype(float)
    out["primes_in_next_100_z"] = (np.searchsorted(ap, k + 100) - pos - 1 - 100 / np.log(kf)) / np.sqrt(100 / np.log(kf))
    s3 = ap[:-2] + ap[1:-1] + ap[2:]
    out["sum_3_consecutive_primes"] = np.isin(k, s3).astype(float)
    # Fermat quotients q_b = (b^(k-1)-1)/k mod k, as a fraction of k
    for b in (2, 3, 5):
        out[f"fermat_quotient_{b}"] = np.fromiter(
            ((int(gmpy2.powmod(b, x - 1, x * x)) - 1) // x / x if x > b else np.nan for x in k.tolist()),
            dtype=np.float64, count=len(k))
    # Euler quotient of the LL seed: ((2+sqrt3)^(k-(3/k)) - 1)/k: skip; Wilson-type sign: ((k-1)/2)! for k=3 mod 4
    return out


# ----------------------------------------------------------------------------- sampled: real quadratic fields
def g_realquad(kvals):
    """kvals: any int array (not necessarily contiguous). Returns h+ and regulator features."""
    g = pari()
    g("Q=[" + ",".join(str(int(x)) for x in kvals) + "]")
    r = gpvec("vector(#Q,i,my(p=Q[i],D=if(p%4==1,p,4*p),c=quadclassunit(D));[c.no,c.reg*1.])", 2)
    D = np.where(kvals % 4 == 1, kvals, 4 * kvals).astype(np.float64)
    return {"realquad_h_gt1": (r[:, 0] > 1).astype(float),
            "realquad_reg_norm": r[:, 1] / np.sqrt(D),
            "realquad_hR_norm": r[:, 0] * r[:, 1] / np.sqrt(D)}


GROUPS = {"rep": g_rep, "bases": g_bases, "beatty": g_beatty, "cf": g_cf, "forms": g_forms, "orders": g_orders,
          "classno": g_classno, "ell": g_ell, "neigh": g_neigh, "dyn": g_dyn}

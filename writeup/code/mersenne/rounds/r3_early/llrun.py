"""Stream statistics of the Lucas-Lehmer run s_0 = 4, s_i = s_{i-1}^2 - 2 and of the
PRP base-3 run x_0 = 3, x_i = x_{i-1}^2, both mod M = 2^k - 1, for every survivor k.
Nothing but running sums and checkpoint snapshots is kept.

usage: llrun.py survivors.txt kmin kmax out.pkl
"""
import sys, math, pickle, time
import gmpy2
from gmpy2 import mpz, popcount, gcd, jacobi
from multiprocessing import Pool

PR = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
PPROD = 1
for p in PR:
    PPROD *= p
W32 = (1 << 32) - 1
F61 = (1 << 61) - 1
JAC_MAX = 512           # Jacobi symbol of every s_i only up to this i (cost)


def horizons(k):
    return {'j32': 32, 'j64': 64, 'j128': 128, 'j256': 256,
            'k/16': k // 16, 'k/8': k // 8, 'k/4': k // 4, 'k/2': k // 2}


def run_seq(k, kind):
    M = (mpz(1) << k) - 1
    sk = math.sqrt(k) / 2.0
    half = k / 2.0
    H = horizons(k)
    at = {}
    for name, j in H.items():
        at.setdefault(j, []).append(name)
    L = k - 2 if kind == 'LL' else k
    v = mpz(4) if kind == 'LL' else mpz(3)
    prev = v
    i0 = None
    # running sums (from i0 on)
    n = 0
    S = dict(hw=0.0, hw2=0.0, run=0.0, run2=0.0, lag=0.0, top=0.0, bot=0.0, r32=0.0,
             jac=0.0, njac=0)
    maxabs = 0.0
    prevhwz = None
    hist = {p: [0] * p for p in PR}
    tort = v
    prod = mpz(1)
    brent_first = None      # (i_checkpoint, factor)
    fp = set()
    repeat_at = None
    snaps = {}
    late = {}
    for i in range(1, L + 1):
        if kind == 'LL':
            x = v * v - 2
            if x < 0:
                x += M
        else:
            x = v * v
        if x >= M:
            if i0 is None:
                i0 = i
            x = (x & M) + (x >> k)
            if x >= M:
                x -= M
        prev, v = v, x
        # Brent cycle signal: prod of (v - s_{2^t}) for 2^t < i <= 2^(t+1)
        d = v - tort
        if d < 0:
            d += M
        y = prod * d
        y = (y & M) + (y >> k)
        if y >= M:
            y -= M
        prod = y
        if i & (i - 1) == 0:
            tort = v
        h = int(v % F61)
        if repeat_at is None:
            if h in fp:
                repeat_at = i
            fp.add(h)
        if i0 is not None:
            hw = popcount(v)
            tr = popcount(v ^ (v >> 1))
            hz = (hw - half) / sk
            rz = (tr - half) / sk
            n += 1
            S['hw'] += hz; S['hw2'] += hz * hz
            S['run'] += rz; S['run2'] += rz * rz
            if prevhwz is not None:
                S['lag'] += hz * prevhwz
            prevhwz = hz
            if abs(hz) > maxabs:
                maxabs = abs(hz)
            S['top'] += int(v >> (k - 32)) / 4294967296.0
            S['bot'] += int(v & W32) / 4294967296.0
            S['r32'] += int(v % W32) / W32
            R = int(v % PPROD)
            for p in PR:
                hist[p][R % p] += 1
            if i <= JAC_MAX:
                S['jac'] += jacobi(v, M); S['njac'] += 1
        if i in at:
            f = {}
            f['i0'] = i0
            hw = popcount(v)
            f['v_hwz'] = (hw - half) / sk
            f['v_runz'] = (popcount(v ^ (v >> 1)) - half) / sk
            f['v_dist_pm2m'] = min(hw, k - hw)
            f['v_top'] = int(v >> (k - 32)) / 4294967296.0
            f['v_bot'] = int(v & W32) / 4294967296.0
            f['v_r32'] = int(v % W32) / W32
            R = int(v % PPROD)
            for p in PR:
                f['v_mod%d' % p] = (R % p) / p
            f['v_jac'] = jacobi(v, M)
            if kind == 'LL':
                f['chk_jac_s+2'] = jacobi(v + 2, M)
                f['chk_jac_s-2'] = jacobi((v - 2) % M, M)
            pr = (v * prev) % M
            f['v_prodprev_hwz'] = (popcount(pr) - half) / sk
            f['v_diffprev_hwz'] = (popcount((v - prev) % M) - half) / sk
            g = gcd(prod, M)
            f['v_brent'] = 1 if 1 < g < M else (2 if g == M else 0)
            if f['v_brent'] == 1 and brent_first is None:
                brent_first = (i, int(g))
            f['v_repeat'] = 1 if repeat_at is not None else 0
            # aggregates
            f['a_n'] = n
            nn = max(n, 1)   # before the first reduction (n = 0) all sums are 0
            f['a_hw_meanz'] = S['hw'] / math.sqrt(nn)
            f['a_hw_sd'] = math.sqrt(max(S['hw2'] / nn - (S['hw'] / nn) ** 2, 0))
            f['a_hw_maxabs'] = maxabs
            f['a_run_meanz'] = S['run'] / math.sqrt(nn)
            f['a_run_sd'] = math.sqrt(max(S['run2'] / nn - (S['run'] / nn) ** 2, 0))
            f['a_hw_lag1z'] = S['lag'] / math.sqrt(max(n - 1, 1))
            f['a_top_z'] = (S['top'] / nn - 0.5) * math.sqrt(12 * n)
            f['a_bot_z'] = (S['bot'] / nn - 0.5) * math.sqrt(12 * n)
            f['a_r32_z'] = (S['r32'] / nn - 0.5) * math.sqrt(12 * n)
            chi = 0.0
            df = 0
            for p in PR:
                e = nn / p
                c0 = hist[p][0]
                f['a_zero%d_z' % p] = (c0 - e) / math.sqrt(e * (1 - 1 / p)) if n else 0.0
                chi += sum((c - e) ** 2 for c in hist[p]) / e if n else 0.0
                df += p - 1
            f['a_chi2_z'] = (chi - df) / math.sqrt(2 * df) if n else 0.0
            f['a_jac_z'] = S['jac'] / math.sqrt(S['njac']) if S['njac'] else 0.0
            for name in at[i]:
                snaps[name] = f
        if i >= L - 5:
            late[i] = popcount(v)
    res = {'i0': i0, 'snaps': snaps, 'late_hw': late, 'repeat_at': repeat_at,
           'brent_first': brent_first, 'final_brent_gcd': int(gcd(prod, M)) if 1 < gcd(prod, M) < M else 1}
    if kind == 'LL':
        # here v = s_{k-2}, prev = s_{k-3}
        r = pow(mpz(2), (k + 1) // 2, M)
        res['s_k2_zero'] = (v == 0)
        res['s_k3_sign'] = 1 if prev == r else (-1 if prev == M - r else 0)
    else:
        # v = x_k, prev = x_{k-1}
        res['x_k1_is_m3'] = (prev == M - 3)
        res['x_k_is_9'] = (v == 9)
    return res


def work(k):
    t = time.process_time()
    out = {'k': k, 'LL': run_seq(k, 'LL'), 'PRP': run_seq(k, 'PRP')}
    out['cpu'] = time.process_time() - t
    return out


if __name__ == '__main__':
    fn, kmin, kmax, outfn = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    ks = [int(x) for x in open(fn) if kmin <= int(x) <= kmax]
    ks.sort(reverse=True)
    t0 = time.time()
    with Pool(2) as pool:
        res = list(pool.imap_unordered(work, ks, chunksize=1))
    res.sort(key=lambda r: r['k'])
    pickle.dump(res, open(outfn, 'wb'))
    print('k', kmin, kmax, 'runs', len(res), 'wall %.1f s' % (time.time() - t0),
          'cpu %.1f s' % sum(r['cpu'] for r in res))

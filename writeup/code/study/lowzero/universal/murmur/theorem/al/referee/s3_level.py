"""Referee steps 3-4: N = 2q, k = 2. t_q, t_N from class numbers (Martin Prop 5.5, validated in s1);
t_1, t_2 from PARI's trace form of the new space (independent of the formula). Window means vs densities."""
import sys, math, json, time
import numpy as np
import cypari2
from sympy import primerange, factorint
pari = cypari2.Pari(); pari.allocatemem(700 * 10 ** 6)

C0 = float(np.exp(np.sum(np.log1p(-1.0 / ((np.array(list(primerange(2, 10 ** 6)), float) - 1) ** 2 * (np.array(list(primerange(2, 10 ** 6)), float) + 1))))))

def e_gen(s, l):
    return l * l / (l * l - 1) if s % l == 0 else 1 - 1 / ((l - 1) ** 2 * (l + 1))

def nup(s):
    v = C0
    for l in factorint(s):
        v *= e_gen(0, l) / e_gen(1, l)
    return v

def Dq(y):
    v = math.pi * math.sqrt(2 * y) - 24 * y; s = 1
    while s * s < 8 * y:
        v += 12 / math.pi * (2 if s % 2 else 0.5) * nup(s) * math.sqrt(8 * y - s * s); s += 1
    return v

def DN(y):
    v = -3 * math.pi * math.sqrt(y) + 24 * y; r = 1
    while r * r < 4 * y:
        v -= 24 / math.pi * (1.5 if r % 2 else 0.75) * nup(r) * math.sqrt(4 * y - r * r); r += 1
    return v

hcl = lambda n: float(pari.qfbhclassno(n))

def kron2(D):  # Kronecker (D/2), D odd
    return 1 if D % 8 in (1, 7) else -1

WIN = [(0.25 * i, 0.25 * (i + 1)) for i in range(16)]

def wmean(ps, N, X):
    y = ps / N; w = np.log(ps); out = []
    for a, b in WIN:
        m = (y > a) & (y <= b)
        out.append(float(np.sum(w[m] * X[m]) / np.sum(w[m])))
    return np.array(out)

def level(q, full=True):
    N = 2 * q
    ps = np.array([p for p in primerange(3, 4 * N + 1) if p != q], dtype=np.int64)
    dim = int(pari('mfdim([%d,2],0)' % N))
    t0 = time.time()
    tf = pari('mfcoefs(mftraceform([%d,2],0),%d)' % (N, 8 * N + 8))
    t1 = np.array([float(tf[int(p)]) for p in ps]); t2 = np.array([-float(tf[2 * int(p)]) for p in ps])
    trW = {Q: float(pari('my(mf=mfinit([%d,2],0),A=mfatkininit(mf,%d)); trace(A[2]/A[3])' % (N, Q))) for Q in (2, q, N)} if N < 3000 else None
    out = {'N': N, 'dim': dim, 'ps': ps, 't1': t1, 't2': t2, 'time_tf': time.time() - t0}
    if not full:
        return out
    tq = np.zeros(len(ps)); tN = np.zeros(len(ps))
    smean = {s: [0.0, 0] for s in range(4)}
    for i, p in enumerate(ps.tolist()):
        tot = 0.0; s = -int(math.isqrt(4 * p // q)) - 1
        for s in range(s, -s + 1):
            if s * s * q < 4 * p:
                A = 4 * p * q - s * s * q * q; Delta = -A
                if A % 4 == 0:
                    xh = 2 * hcl(A // 4) - hcl(A)
                else:
                    xh = (kron2(Delta) - 1) * hcl(A)
                tot += xh
                if s >= 0 and s < 4:
                    smean[s][0] += math.log(p) * math.pi * xh / math.sqrt(A); smean[s][1] += math.log(p)
            s += 1
        tq[i] = -tot / 2 - (p + 1)
        tot = 0.0; r = -int(math.isqrt(4 * p // N)) - 1
        for r in range(r, -r + 1):
            if r * r * N < 4 * p:
                tot += hcl(4 * p * N - r * r * N * N)
            r += 1
        tN[i] = -tot / 2 + (p + 1)
    out.update({'tq': tq, 'tN': tN, 'smean': {s: v[0] / v[1] for s, v in smean.items()}})
    return out

if __name__ == '__main__':
    t0 = time.time()
    res = {}
    for q in (10103, 10111, 9973):
        L = level(q)
        N = L['N']; ps = L['ps']; dim = L['dim']
        y = ps / N
        mq = wmean(ps, N, L['tq'] / dim); mN = wmean(ps, N, L['tN'] / dim)
        dq = wmean(ps, N, np.array([Dq(v) for v in y])); dn = wmean(ps, N, np.array([DN(v) for v in y]))
        m1 = wmean(ps, N, L['t1'] / dim); m2 = wmean(ps, N, L['t2'] / dim)
        res[N] = dict(mq=mq.tolist(), mN=mN.tolist(), dq=dq.tolist(), dn=dn.tolist(), m1=m1.tolist(), m2=m2.tolist(), dim=dim,
                      smean=L['smean'])
        # class-number means per s for Q = q vs E^(q,2)(s)
        Eq = {0: math.pi ** 2 / 6 * (1 - 1 / q ** 2) * (-0.5)}
        for s in (1, 2, 3):
            Eq[s] = nup(s) / e_gen(s, 2) / e_gen(s, q) * ((-4 / 3) if s % 2 else (-2 / 3))
        print('N = %d dim %d primes %d (trace form %.0fs, total %.0fs)' % (N, dim, len(ps), L['time_tf'], time.time() - t0), flush=True)
        print('   class-number means pi*xi*H/sqrt(D), Q = q: ' + '  '.join('s=%d %.4f (E %.4f)' % (s, L['smean'][s], Eq[s]) for s in range(4)))
    json.dump(res, open('s3_out.json', 'w'))
    Ns = list(res)
    for key, mk in (('mq', 'dq'), ('mN', 'dn')):
        R = np.array([np.array(res[N][key]) - np.array(res[N][mk]) for N in Ns])
        sig = R.std(axis=0, ddof=1)
        print('%s vs %s: residual window means per level (rows) ' % (key, mk))
        for N, r in zip(Ns, R):
            print('   %d: ' % N + ' '.join('%6.3f' % v for v in r))
        print('   mean residual / (sigma/sqrt3): ' + ' '.join('%5.1f' % v for v in R.mean(axis=0) / (sig / math.sqrt(3))))
        print('   chi2 of the mean residual (16 windows) = %.1f; rms residual %.3f; rms density %.3f' % (
            float(np.sum((R.mean(axis=0) / (sig / math.sqrt(3))) ** 2)), float(np.sqrt(np.mean(R ** 2))), float(np.sqrt(np.mean(np.array([res[N][mk] for N in Ns]) ** 2)))))
    N0 = Ns[0]
    print('data  t_q/dim  (N=%d): ' % N0 + ' '.join('%6.3f' % v for v in res[N0]['mq']))
    print('model D_q         : ' + ' '.join('%6.3f' % v for v in res[N0]['dq']))
    print('data  t_N/dim     : ' + ' '.join('%6.3f' % v for v in res[N0]['mN']))
    print('model D_N         : ' + ' '.join('%6.3f' % v for v in res[N0]['dn']))
    for N in Ns:
        print('t1/dim N=%d: rms %.4f  ' % (N, float(np.sqrt(np.mean(np.array(res[N]['m1']) ** 2)))) + ' '.join('%6.3f' % v for v in res[N]['m1']))
        print('t2/dim N=%d: rms %.4f  ' % (N, float(np.sqrt(np.mean(np.array(res[N]['m2']) ** 2)))) + ' '.join('%6.3f' % v for v in res[N]['m2']))
    print('total %.0fs' % (time.time() - t0))

"""Round 7 (preregister/round7.txt): the gap sequence as a hidden Markov chain on the residues of p mod Q.
usage: r7.py build Q | predict Q"""
import sys, os, math, json, pickle, time
from math import gcd
import numpy as np
import theory2 as T2
import theory_r3 as R3

OUT = R3.OUT
HM = 300
CACHE7 = os.environ.get('STUDY_CACHE7', '/tmp/study_r7.pkl')


def units(Q):
    return [a for a in range(1, Q) if gcd(a, Q) == 1] if Q > 1 else [0]


def build(Q):
    cache = pickle.load(open(CACHE7, 'rb')) if os.path.exists(CACHE7) else {}
    t0 = time.time()
    todo = [(Q, a, h) for h in range(2, HM + 1, 2) for a in units(Q) if (Q, a, h) not in cache]
    for key in todo:
        _, a, h = key
        if gcd(a + h, Q) != 1:
            cache[key] = (0.0, None); continue
        o = R3.order_for(h)
        pt = T2.GPattern([0, h], o, T2.S_cond(a, Q))
        cache[key] = (pt.SB, np.array([pt.log_avoid(l, o) for l in R3.LAMG]) if pt.SB > 0 else None)
        if time.time() - t0 > 450:
            break
    pickle.dump(cache, open(CACHE7, 'wb'))
    done = sum((Q, a, h) in cache for h in range(2, HM + 1, 2) for a in units(Q))
    print('Q', Q, 'cached', done, 'of', len(units(Q)) * (HM // 2), flush=True)


def trans(cache, Q, L, W):
    U = units(Q); idx = {a: i for i, a in enumerate(U)}
    S0 = {a: float(T2.S_cond(a, Q)(np.array([[0]]))[0]) for a in U}
    u = np.linspace(math.log(L), math.log(L + W), 401); w = np.exp(u) / u
    T = np.zeros((20, len(U), len(U)))
    for a in U:
        rows = []
        for h in range(2, HM + 1, 2):
            SB, la = cache[(Q, a, h)]
            if SB == 0:
                continue
            f = SB / S0[a] / u * np.exp(np.interp(u, R3.LAMG, la))
            rows.append((h, float(np.trapezoid(w * f, u) / np.trapezoid(w, u))))
        tot = sum(v for _, v in rows)
        for h, v in rows:
            T[min(h // 2, 20) - 1, idx[a], idx[(a + h) % Q]] += v / tot
    return T


def H(p):
    p = np.asarray(p).ravel(); p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def predict(Q):
    cache = pickle.load(open(CACHE7, 'rb'))
    res = json.load(open(os.path.join(OUT, 'r7.json'))) if os.path.exists(os.path.join(OUT, 'r7.json')) else {}
    meas = {}
    for l in open(os.path.join(OUT, 'results_r6.jsonl')):
        r = json.loads(l)
        if r['block'] == -1:
            meas[r['window']] = r['stats']
    for win, (L, W) in R3.WIN.items():
        T = trans(cache, Q, L, W)
        M = T.sum(0)
        ev, V = np.linalg.eig(M.T)
        pi = np.real(V[:, np.argmin(abs(ev - 1))]); pi /= pi.sum()
        P1 = np.array([pi @ T[c] @ np.ones(len(pi)) for c in range(20)])
        P2 = np.array([[pi @ T[a] @ T[b] @ np.ones(len(pi)) for b in range(20)] for a in range(20)])
        A = np.array([pi @ T[a] for a in range(20)])                       # (20, states)
        AB = np.einsum('as,bst->abt', A, T)                                 # (20, 20, states)
        P3 = np.einsum('abs,cst->abc', AB, T)                              # t summed: (20, 20, 20)
        h1 = H(P2) - H(P1)
        h2 = H(P3) - H(P2)
        key = '%s_Q%d' % (win, Q)
        res[key] = {'H1_pred': h1, 'H2_pred': h2, 'H1_meas': meas[win]['H1_mm'], 'H2_meas': meas[win]['H2_mm'],
                    'stationary_max_dev': float(np.max(abs(pi - 1 / len(pi))))}
        print(key, {k: round(v, 4) for k, v in res[key].items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'r7.json'), 'w'), indent=1)


if __name__ == '__main__':
    {'build': build, 'predict': predict}[sys.argv[1]](int(sys.argv[2]))

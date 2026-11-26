"""Blind test: murmurations on Atkin-Lehner eigenspaces at one level N = p1 q (p1 = 2, 3), weight 2 (al_theorem.md).

Usage (from this folder):
  BUDGET=540 python -I -B al_blind_test.py count blind   # resumable: repeat until it prints 'all levels done'
  python -I -B al_blind_test.py read blind               # decisions D1-D3 (refuses unless every level is done)
  ... count dry / read dry                               # the same on the six discovery levels (in-sample dry run)
The predictions file must have the registered sha256, else the script stops.

Data per level (own code; Martin, arXiv 2409.02338, trace formula for T_p W_Q on S_2^new(N), N squarefree, validated against
PARI in al_check.py C1-C3): for every prime p <= 4N, p !| N, and Q in {1, p1, q, N},
  t_Q(p) = -(1/2) sum_{s in Z, s^2 Q < 4p} xi_{s^2Q^2 - 4Qp}(N/Q) H(4Qp - s^2Q^2) + mu(N/Q)(p + 1),
  xi_D(l) = l H(|D|/l^2)/H(|D|) - 1 if l^2 | D, else (D/l) - 1 (multiplicative in N/Q);
eigenspace mean E^eps(p) = (t_1 + e1 t_p1 + eq t_q + e1 eq t_N) / (4 dim^eps), eps = (e1, eq) = signs at (p1, q);
antisymmetric parts (t_q + e1 t_N)/dim; symmetric parts (t_1 + e1 t_p1)/dim. Window means with weights log p over the
16 windows (0, 0.25], ..., (3.75, 4] of y = p/N.
"""
import sys, os, json, time, math, hashlib
import cypari2

HERE = os.path.dirname(os.path.abspath(__file__))
SHA = {'blind': 'b9aa428ce991e062cacfda8510e2be48216c77f1259c1c194eccfb082deccc91', 'dry': '6b0f27e7b171a46bffd48eed5d0b2fb92aa9a4e1f8a144cc6e177e5147053eba'}
FILES = {'blind': 'predictions_blind.json', 'dry': 'predictions_dry.json'}
EPS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]


def load(which):
    txt = open(os.path.join(HERE, FILES[which]), 'rb').read()
    h = hashlib.sha256(txt).hexdigest()
    if h != SHA[which]:
        sys.exit('REFUSED: %s has sha256 %s, registered %s' % (FILES[which], h, SHA[which]))
    return json.loads(txt)


def primes_upto(n):
    sv = bytearray([1]) * (n + 1); sv[0] = sv[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = bytearray(len(sv[i * i::i]))
    return [i for i in range(n + 1) if sv[i]]


def kron_l(D, l, sq):  # Kronecker symbol (D/l) for l = 2 or an odd prime l (sq: set of squares mod l)
    if l == 2:
        if D % 2 == 0:
            return 0
        return 1 if D % 8 in (1, 7) else -1
    r = D % l
    return 0 if r == 0 else (1 if r in sq else -1)


def level_data(N, p1, q, pari):
    hcl = pari('D -> qfbhclassno(D)')
    nmax = 16 * p1 * N + 8
    Ht = [-1.0 / 12] + [float(x) for x in pari('vector(%d,n,qfbhclassno(n))' % nmax)]
    sqq = set(x * x % q for x in range(1, q)); sq1 = set(x * x % p1 for x in range(1, p1))
    mu_q = -1; mu_p1 = -1
    dim = int(pari('mfdim([%d,2],0)' % N))
    out = []
    for p in primes_upto(4 * N):
        if N % p == 0:
            continue
        t = {}
        # Q = 1 and Q = p1: discriminants up to 16 p1 N (table); xi at q is (D/q) - 1 since q^2 > |D|
        for Q0 in (1, p1):
            tot = 0.0; s = 0
            while s * s * Q0 < 4 * p:
                D = s * s * Q0 * Q0 - 4 * Q0 * p; A = -D
                if Q0 == 1:
                    if A % (p1 * p1) == 0:
                        v = p1 * Ht[A // (p1 * p1)] - Ht[A]
                    else:
                        v = (kron_l(D, p1, sq1) - 1) * Ht[A]
                else:
                    v = Ht[A]
                v *= kron_l(D, q, sqq) - 1
                tot += 0.5 * v if s == 0 else v
                s += 1
            t[Q0] = -tot + (p + 1) * (mu_p1 * mu_q if Q0 == 1 else mu_q)
        for Q in (q, N):
            M = N // Q; tot = 0.0; s = 0
            while s * s * Q < 4 * p:
                A = 4 * p * Q - s * s * Q * Q
                if M == 1:
                    v = float(hcl(A))
                elif A % (M * M) == 0:
                    v = M * float(hcl(A // (M * M))) - float(hcl(A))
                else:
                    v = (kron_l(-A, M, sq1) - 1) * float(hcl(A))
                tot += 0.5 * v if s == 0 else v
                s += 1
            t[Q] = -tot + (p + 1) * (1 if M == 1 else mu_p1)
        out.append((p, t))
    return out, dim


def count(which):
    doc = load(which)
    pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
    state_f = os.path.join(HERE, 'state_%s.json' % which)
    state = json.load(open(state_f)) if os.path.exists(state_f) else {}
    budget = float(os.environ.get('BUDGET', '540')); t0 = time.time()
    edges = doc['meta']['windows']
    todo = [str(N) for T in sorted(doc['levels']) for N in doc['levels'][T] if str(N) not in state]
    for key in todo:
        if time.time() - t0 > budget:
            break
        N = int(key); P = doc['pred'][key]; p1, q = P['p1'], P['q']; tl = time.time()
        data, dim = level_data(N, p1, q, pari)
        if dim != P['dim']:
            sys.exit('dimension mismatch at N = %d' % N)
        de = P['dimeps']
        nq = 4 + 2 + 2
        num = [[0.0] * len(edges) for _ in range(nq)]; den = [0.0] * len(edges)
        for p, t in data:
            y = p / N; w = math.log(p)
            vals = [(t[1] + e1 * t[p1] + eq * t[q] + e1 * eq * t[N]) / (4 * de[k]) for k, (e1, eq) in enumerate(EPS)]
            vals += [(t[q] + e1 * t[N]) / dim for e1 in (1, -1)] + [(t[1] + e1 * t[p1]) / dim for e1 in (1, -1)]
            for i, (a, b) in enumerate(edges):
                if a < y <= b:
                    for j in range(nq):
                        num[j][i] += w * vals[j]
                    den[i] += w
                    break
        mean = [[num[j][i] / den[i] for i in range(len(edges))] for j in range(nq)]
        state[key] = {'data': mean[0:4], 'anti': mean[4:6], 'sym': mean[6:8], 'dim': dim, 'sec': time.time() - tl}
        with open(state_f + '.tmp', 'w') as fh:
            json.dump(state, fh)
        os.replace(state_f + '.tmp', state_f)
        print('N = %d done (%.1f s)' % (N, time.time() - tl), flush=True)
    left = sum(1 for T in doc['levels'] for N in doc['levels'][T] if str(N) not in state)
    print('all levels done' if left == 0 else '%d levels left: run again' % left)


def read(which):
    doc = load(which)
    state_f = os.path.join(HERE, 'state_%s.json' % which)
    state = json.load(open(state_f)) if os.path.exists(state_f) else {}
    missing = [N for T in doc['levels'] for N in doc['levels'][T] if str(N) not in state]
    if missing:
        sys.exit('REFUSED: %d levels not counted yet' % len(missing))
    dec = doc['meta']['decision']; allpass = True
    for T in sorted(doc['levels']):
        sig = doc['meta']['sigma'][T]; siga = doc['meta']['sigma_anti'][T]; W = len(sig[0])
        c2 = {'own': 0.0, 'wn': 0.0, 'e0': 0.0, 'e1': 0.0, 'e2': 0.0, 'e3': 0.0, 'aown': 0.0, 'awn': 0.0}
        nt = 0; nta = 0; rms = 0.0; srms = 0.0; ns = 0
        for N in doc['levels'][T]:
            d = state[str(N)]; P = doc['pred'][str(N)]
            for k in range(4):
                for w in range(W):
                    x = d['data'][k][w]; nt += 1
                    c2['own'] += ((x - P['own'][k][w]) / sig[k][w]) ** 2
                    c2['wn'] += ((x - P['wn'][k][w]) / sig[k][w]) ** 2
                    rms += (x - P['own'][k][w]) ** 2
                    for j in range(4):
                        if j != k:
                            c2['e%d' % j] += ((x - P['alt'][k][j][w]) / sig[k][w]) ** 2
            for i in range(2):
                for w in range(W):
                    x = d['anti'][i][w]; nta += 1
                    c2['aown'] += ((x - P['anti_model'][i][w]) / siga[i][w]) ** 2
                    c2['awn'] += ((x - P['anti_wn'][i][w]) / siga[i][w]) ** 2
                    srms += d['sym'][i][w] ** 2; ns += 1
        alts = ['wn', 'e0', 'e1', 'e2', 'e3']
        nl = len(doc['levels'][T])
        d1 = c2['own'] <= dec['D1_own_chi2_per_term_max'] * nt
        d2 = all(c2[m] - c2['own'] >= dec['D2_alt_delta_per_level_min'] * nl and c2[m] > dec['D2_alt_chi2_per_term_min'] * nt for m in alts)
        d3 = c2['aown'] <= dec['D3_anti_own_chi2_per_term_max'] * nta and c2['awn'] - c2['aown'] >= dec['D3_anti_wn_delta_per_level_min'] * nl
        allpass &= d1 and d2 and d3
        print('type %s (%d levels): eigenspaces (%d terms) chi2 own %.1f [D1 %s]; W_N-only %.1f; other eigenspaces %s [D2 %s]' % (
            T, len(doc['levels'][T]), nt, c2['own'], 'PASS' if d1 else 'FAIL', c2['wn'],
            ', '.join('%s %.0f' % (EPS[j], c2['e%d' % j]) for j in range(4)), 'PASS' if d2 else 'FAIL'))
        print('          antisymmetric parts (%d terms): chi2 own %.1f, W_N-only %.1f [D3 %s]; rms(data - own) %.3f; '
              'reported: rms of the symmetric parts %.3f' % (nta, c2['aown'], c2['awn'], 'PASS' if d3 else 'FAIL',
                                                            math.sqrt(rms / nt), math.sqrt(srms / ns)))
    print('VERDICT: %s' % ('PASS (D1, D2, D3 for both types)' if allpass else 'FAIL'))


if __name__ == '__main__':
    mode, which = sys.argv[1], sys.argv[2]
    {'count': count, 'read': read}[mode](which)

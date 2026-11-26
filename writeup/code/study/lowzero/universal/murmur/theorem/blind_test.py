"""Blind test of the fixed-level murmuration density M^(N) (level_theorem.md, Theorem 5), weight 2.

Usage (from this folder):
  BUDGET=540 python -I -B blind_test.py count blind     # resumable: repeat until it prints 'all levels done'
  python -I -B blind_test.py read blind                 # decisions D1, D2 (refuses unless all levels are done)
  ... count dry / read dry                              # the same on the 20 discovery levels (in-sample dry run)
The predictions file must have the registered sha256, else the script stops. Data per level: for every prime p <= 4N,
p !| N, T(p) = H(4pN)/2 + sum_{r>=1, r^2 N < 4p} H(4pN - r^2 N^2) - (p + 1) = sum_{f new} eps_f a_f(p) (Skoruppa-Zagier),
E_N(p) = T(p)/dim S_2^new(N); window means over 16 windows of y = p/N in (0, 4] with weights log p.
"""
import sys, os, json, time, math, hashlib
from fractions import Fraction as Fr
import cypari2

HERE = os.path.dirname(os.path.abspath(__file__))
SHA = {'blind': '2f67dc3a2f7fcf86327029c9a0db964cc0bbd21fb11a2fa0d6c6cd27f568a929', 'dry': '99ae368442fa4875b43509487a95bc4993b3045aa1024ea2e9f3dce518ab97db'}
FILES = {'blind': 'predictions_blind.json', 'dry': 'predictions_dry.json'}


def load(which):
    txt = open(os.path.join(HERE, FILES[which]), 'rb').read()
    h = hashlib.sha256(txt).hexdigest()
    if h != SHA[which]:
        sys.exit('REFUSED: %s has sha256 %s, registered %s' % (FILES[which], h, SHA[which]))
    return json.loads(txt)


def count(which):
    doc = load(which)
    pari = cypari2.Pari(); pari.allocatemem(200 * 10 ** 6)
    hcl = pari('D -> qfbhclassno(D)')
    state_f = os.path.join(HERE, 'state_%s.json' % which)
    state = json.load(open(state_f)) if os.path.exists(state_f) else {}
    budget = float(os.environ.get('BUDGET', '540')); t0 = time.time()
    edges = doc['meta']['windows']
    todo = [str(N) for t in sorted(doc['levels']) for N in doc['levels'][t] if str(N) not in state]
    for key in todo:
        if time.time() - t0 > budget:
            break
        N = int(key); tl = time.time()
        dim = int(pari('mfdim([%d,2],0)' % N))
        if dim != doc['pred'][key]['dim']:
            sys.exit('dimension mismatch at N = %d' % N)
        num = [0.0] * len(edges); den = [0.0] * len(edges); cnt = [0] * len(edges)
        for p in pari('primes([2, %d])' % (4 * N)):
            p = int(p)
            if N % p == 0:
                continue
            T = Fr(str(hcl(4 * p * N))) / 2 - (p + 1)
            r = 1
            while r * r * N < 4 * p:
                T += Fr(str(hcl(4 * p * N - r * r * N * N)))
                r += 1
            y = p / N; w = math.log(p)
            for i, (a, b) in enumerate(edges):
                if a < y <= b:
                    num[i] += w * float(T) / dim; den[i] += w; cnt[i] += 1
                    break
        state[key] = {'mean': [n / d for n, d in zip(num, den)], 'nprimes': cnt, 'dim': dim, 'sec': time.time() - tl}
        with open(state_f + '.tmp', 'w') as fh:
            json.dump(state, fh)
        os.replace(state_f + '.tmp', state_f)
        print('N = %d done (%.1f s)' % (N, time.time() - tl), flush=True)
    left = sum(1 for t in doc['levels'] for N in doc['levels'][t] if str(N) not in state)
    print('all levels done' if left == 0 else '%d levels left: run again' % left)


def read(which):
    doc = load(which)
    state_f = os.path.join(HERE, 'state_%s.json' % which)
    state = json.load(open(state_f)) if os.path.exists(state_f) else {}
    dec = doc['meta']['decision']
    missing = [N for t in doc['levels'] for N in doc['levels'][t] if str(N) not in state]
    if missing:
        sys.exit('REFUSED: %d levels not counted yet' % len(missing))
    allpass = True
    for t in sorted(doc['levels']):
        sig = doc['meta']['sigma'][t]
        models = ['own', 'zub'] + ['type_' + x for x in ('q', '2q', '3q', '6q') if x != t]
        c2 = {m: 0.0 for m in models}; c2h = {m: [0.0, 0.0] for m in models}; rms = 0.0; nt = 0
        for N in doc['levels'][t]:
            d = state[str(N)]['mean']; P = doc['pred'][str(N)]
            for m in models:
                for w in range(len(d)):
                    z2 = ((d[w] - P[m][w]) / sig[w]) ** 2
                    c2[m] += z2; c2h[m][w >= len(d) // 2] += z2
            rms += sum((d[w] - P['own'][w]) ** 2 for w in range(len(d))); nt += len(d)
        d1 = c2['own'] <= dec['D1_own_chi2_max']
        d2 = all(c2[m] - c2['own'] >= dec['D2_alt_delta_chi2_min'] and c2[m] > dec['D2_alt_chi2_min'] for m in models if m != 'own')
        allpass &= d1 and d2
        print('type %-3s (%d levels, %d terms): chi2 own %.1f [D1 %s]; ' % (t, len(doc['levels'][t]), nt, c2['own'], 'PASS' if d1 else 'FAIL')
              + ', '.join('%s %.1f' % (m, c2[m]) for m in models if m != 'own') + ' [D2 %s]' % ('PASS' if d2 else 'FAIL'))
        print('          first/last 8 windows: ' + ', '.join('%s %.1f/%.1f' % (m, c2h[m][0], c2h[m][1]) for m in models)
              + '; rms(data - own) %.4f' % math.sqrt(rms / nt))
    print('VERDICT: %s' % ('PASS (D1 and D2 for all four types)' if allpass else 'FAIL'))


if __name__ == '__main__':
    mode, which = sys.argv[1], sys.argv[2]
    {'count': count, 'read': read}[mode](which)

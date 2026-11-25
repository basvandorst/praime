#!/usr/bin/env python3
"""validate_fft.py -- the FFT statistic of mock_sum.py against direct character sums, from scratch.

For q = 10007 (pilot file pilot/q10007.npz, made with --save-p): 200 random characters chi_k (numpy rng 20261125),
supports s = 0.5, 1.0, 1.5 (all 200) and s = 2.0 (the first 25), those present in the file (plus its largest s). Direct sum over primes from a plain full sieve (count checked by PARI primepi),
PARI's discrete log (znlog to the generator of znstar(q,1), checked = least primitive root used by mock_sum) and plain
cos: P(chi_k) = -(2/L) sum_{p^j <= q^s, p != q} log p p^(-j/2) (1 - j log p/(s L))_+ cos(2 pi k j log_g(p)/(q-1)).
Also checks the family size, the parity split (k even <-> chi(-1) = +1 via chareval) and P(chi) = P(conj chi).
"""
import math, sys
import numpy as np
import cypari2

pari = cypari2.Pari()
pari.allocatemem(800 * 10 ** 6)


def main(path='pilot/q10007.npz'):
    z = np.load(path)
    q = int(z['q']); L = math.log(q); sig = z['sig']; Pfft = z['P']
    g = int(pari(f'znprimroot({q})'))
    assert g == int(z['g']), (g, int(z['g']))
    G = pari(f'znstar({q},1)')
    assert int(pari(f'lift(znstar({q},1).gen[1])')) == g
    rng = np.random.default_rng(20261125)
    ks = rng.choice(np.arange(1, q - 1), size=200, replace=False)
    # parity check via chareval at -1
    for k in ks[:20]:
        e = float(pari.chareval(G, pari(f'[{int(k)}]'), q - 1))
        assert (abs(e) < 1e-12) == (k % 2 == 0)
    # dlog of all residues, from PARI
    dl = np.zeros(q, dtype=np.int64)
    for r in range(1, q):
        dl[r] = int(pari.znlog(r, pari(f'Mod({g},{q})')))
    worst = 0.0
    for s in [x for x in (0.5, 1.0, 1.5, 2.0) if np.min(np.abs(sig - x)) < 1e-9] + ([float(sig.max())] if sig.max() < 1.5 - 1e-9 else []):
        j = int(np.argmin(np.abs(sig - s)))
        assert abs(sig[j] - s) < 1e-9
        X = int(math.floor(q ** s))
        sv = np.ones(X + 1, dtype=bool); sv[:2] = False
        for i in range(2, math.isqrt(X) + 1):
            if sv[i]:
                sv[i * i::i] = False
        pr = np.nonzero(sv)[0].astype(np.int64); del sv
        assert len(pr) == int(pari.primepi(X))          # plain full sieve, count checked by PARI
        pr = pr[pr != q]
        terms_l, terms_w = [], []
        lp = np.log(pr.astype(float))
        w1 = lp / np.sqrt(pr) * np.clip(1 - lp / (s * L), 0, None)
        terms_l.append(dl[pr % q]); terms_w.append(w1)
        for p in pr[pr <= math.isqrt(X)]:
            p = int(p); e = 2
            while p ** e <= X:
                fh = max(0.0, 1 - e * math.log(p) / (s * L))
                if fh > 0:
                    terms_l.append(np.array([dl[pow(p, e, q)]])); terms_w.append(np.array([math.log(p) / p ** (e / 2) * fh]))
                e += 1
        ll = np.concatenate(terms_l); ww = np.concatenate(terms_w)
        kk = ks if s < 2 else ks[:25]
        d = []
        for k in kk:
            Pd = -(2 / L) * np.sum(ww * np.cos(2 * np.pi * ((int(k) * ll) % (q - 1)) / (q - 1)))
            d.append(Pd - Pfft[j, k - 1])
            # conjugate
            assert abs(Pfft[j, k - 1] - Pfft[j, (q - 1 - k) - 1]) < 1e-9
        d = np.abs(np.array(d))
        worst = max(worst, d.max())
        print(f's={s}: {len(kk)} characters, max |direct - FFT| = {d.max():.2e} (|P| ~ {np.std(Pfft[j]):.3f})', flush=True)
    print(f'family size {Pfft.shape[1]} = q-2: {Pfft.shape[1] == q - 2}')
    print('PASS' if worst < 1e-8 else 'FAIL')


if __name__ == '__main__':
    main(*sys.argv[1:])

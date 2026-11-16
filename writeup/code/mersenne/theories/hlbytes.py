"""Theories 3-5: the wheel-24 bytes from the Hardy-Littlewood constants, no fitting.
Byte j = primality of 24j + o, o in (-1, 1, 5, 7, 11, 13, 17, 19), MSB first. For a set U of positions,
P(all of U prime) ~ S(U) * prod 1/ln(24j + o), S(U) = prod_r (1 - nu_r(U)/r) / (1 - 1/r)^|U|, nu_r = number of residues
hit by U mod r. P(exactly T prime) = sum over U >= T of (-1)^(|U|-|T|) P(all of U prime) (inclusion-exclusion).
T3: expected count of every byte value up to 10^7 and 10^9 vs the sieve. T4: bytes 223 and 251 (the two prime
septuplet patterns, starting at 24j - 1). T5: entropy and the information between neighbouring bytes, and the limit
for any predictor of the lane words that knows the residues mod 5..47."""
import numpy as np, math, sys
R = 2_000_000
s = np.ones(R + 1, bool); s[:2] = False
for i in range(2, int(R ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
PR = np.nonzero(s)[0]
O8 = [-1, 1, 5, 7, 11, 13, 17, 19]
def singular(offs):
    n = len(offs); NM = 1 << n
    masks = np.arange(NM)
    m = np.array([bin(x).count('1') for x in range(NM)])
    logS = m * math.log(3.0)
    span = max(offs) - min(offs)
    for r in PR[(PR >= 5) & (PR <= span)]:
        rm = np.zeros(NM, np.int64)
        for i, o in enumerate(offs):
            bit = 1 << (n - 1 - i)
            sel = (masks & bit) != 0
            rm[sel] |= 1 << (o % r)
        nu = np.array([bin(int(x)).count('1') for x in rm]) if NM <= 256 else np.vectorize(lambda x: bin(int(x)).count('1'))(rm)
        with np.errstate(divide='ignore'):
            logS = logS + np.log(np.maximum(1 - nu / r, 0)) - m * math.log(1 - 1 / r)
    big = PR[PR > max(span, 4)]
    tail = np.array([np.sum(np.log1p(-k / big) - k * np.log1p(-1 / big)) for k in range(n + 1)])
    logS = logS + tail[m]
    return np.exp(logS), m
def expected(offs, J, center):
    S, m = singular(offs)
    n = len(offs)
    I = np.zeros(n + 1)
    for a in range(1, J + 1, 5_000_000):
        j = np.arange(a, min(a + 5_000_000, J + 1), dtype=np.float64)
        L = 1 / np.log(24 * j + center)
        for k in range(n + 1): I[k] += np.sum(L ** k)
    h = S * I[m]
    for i in range(n):
        bit = 1 << i
        idx = np.arange(1 << n); lo = idx[(idx & bit) == 0]
        h[lo] -= h[lo | bit]
    return h, S
if __name__ == '__main__':
    a7 = np.load('s1e7.npz'); a9 = np.load('s1e9.npz')
    out = []
    def pr(*x):
        print(*x, flush=True)
    S8, m8 = singular(O8)
    adm = S8 > 0
    pr(f'admissible bytes: {adm.sum()} of 256 (the inadmissible ones: {[int(x) for x in np.nonzero(~adm)[0]]})')
    for tag, a, J in (('10^7', a7, int(a7['J'])), ('10^9', a9, int(a9['J']))):
        E, _ = expected(O8, J, 9)
        Ec, _ = None, None
        S0 = 3.0 ** m8
        obs = a['cnt'].astype(float)
        pr(f'\n== T3 at {tag}: {J:,} bytes ==')
        pr(f'sum of expected {E.sum():,.1f} (should be {J:,}); observed bytes on inadmissible values: {obs[~adm].sum():.0f}; expected there: {E[~adm].sum():.3g}')
        big = E >= 5
        chi = np.sum((obs[big] - E[big]) ** 2 / E[big]); dof = big.sum() - 1
        pr(f'chi-square over {big.sum()} values with E >= 5: {chi:.1f} on {dof} dof (z = {(chi - dof) / math.sqrt(2 * dof):+.2f})')
        rel = (obs - E) / np.sqrt(np.maximum(E, 1e-9))
        o = np.argsort(-np.abs(rel * big))[:6]
        pr('largest standardized deviations: ' + ', '.join(f'{int(v)} ({v:08b}): obs {obs[v]:.0f} exp {E[v]:.1f} z {rel[v]:+.1f}' for v in o))
        pr('top 6 bytes: ' + ', '.join(f'{int(v)}: obs {obs[v]:.0f} exp {E[v]:.0f}' for v in np.argsort(-obs)[:6]))
        for v in (223, 251):
            pr(f'T4 byte {v} ({v:08b}): expected {E[v]:.2f}, observed {obs[v]:.0f}')
        p = obs / obs.sum(); q = E / E.sum()
        Hs = -np.sum(p[p > 0] * np.log2(p[p > 0])); Hq = -np.sum(q[q > 0] * np.log2(q[q > 0]))
        pr(f'T5 entropy per byte: observed {Hs:.4f} bits, predicted {Hq:.4f} bits ({Hs / 4:.4f} vs {Hq / 4:.4f} per word)')
        O16 = O8 + [o + 24 for o in O8]
        E16, _ = expected(O16, J - 1, 21)
        P16 = E16 / E16.sum()
        pj = P16.reshape(256, 256); p1 = pj.sum(1); p2 = pj.sum(0)
        with np.errstate(divide='ignore', invalid='ignore'):
            MIp = np.nansum(pj * np.log2(pj / np.outer(p1, p2)))
        po = a['pair'].astype(float).reshape(256, 256); po /= po.sum(); o1 = po.sum(1); o2 = po.sum(0)
        with np.errstate(divide='ignore', invalid='ignore'):
            MIo = np.nansum(np.where(po > 0, po * np.log2(po / np.outer(o1, o2)), 0))
        nn = a['pair'].sum(); cells = np.sum(po > 0)
        bias = (cells - np.sum(o1 > 0) - np.sum(o2 > 0) + 1) / (2 * nn * math.log(2))
        pr(f'T5 information between neighbouring bytes: observed {MIo:.5f} bits (sampling bias about {bias:.5f}), predicted {MIp:.5f} bits')
        pr(f'   check: the twin pair 24j + 11, 24j + 13 has S = {S8[0b00001100]:.5f} = 3^2 * C2 / (3/4) = {9 * 0.6601618158468696 / 0.75:.5f} (C2 the twin prime constant)')

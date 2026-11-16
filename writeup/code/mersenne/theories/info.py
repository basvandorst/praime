"""Theory 5, second part: what a predictor of the lane words can know, compared with the ML run (train 1..10^6, test near
10^6..2*10^6 and far 9*10^6..10^7). A word is the pair (c - 1, c + 1), c = 0 mod 6 (up when c = 0 mod 12, down when 6).
(a) Residues mod 5..47 (what the ML saw): theory-only predictor: a slot with a factor <= 47 is composite; a survivor is prime
with q(n) = 1/(ln n * prod_{r <= 47}(1 - 1/r)) (Mertens), both survivors of a word prime with q^2 * prod_{r > 47} (1 - 1/(r-1)^2).
Scored on the real words; no training data.
(b) History (the previous 8 words): the Hardy-Littlewood joint law of 9 consecutive words (18 positions), mutual information
between the last word and the 8 before it."""
import numpy as np, math
from hlbytes import singular
N = 10_000_001
s = np.ones(N + 2, bool); s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
P47 = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
mert = 0.5 * (2 / 3)
for r in P47: mert *= 1 - 1 / r
pr = np.nonzero(s)[0]; big = pr[(pr > 47) & (pr < 2_000_000)]
c2 = float(np.prod(1 - 1 / (big - 1.0) ** 2))
def score(lo, hi):
    c = np.arange(((lo // 6) + 1) * 6, hi - 1, 6, dtype=np.int64)
    L = c - 1; R = c + 1
    y = 2 * s[L] + s[R]
    def surv(x):
        ok = np.ones(len(x), bool)
        for r in P47: ok &= x % r != 0
        return ok
    sL, sR = surv(L), surv(R)
    q = 1 / (np.log(c) * mert)
    pL = np.where(sL, q, 0.0); pR = np.where(sR, q, 0.0)
    pB = np.where(sL & sR, q * q * c2, 0.0)
    probs = np.stack([1 - pL - pR + pB, pR - pB, pL - pB, pB], 1)
    probs = np.clip(probs, 1e-12, 1)
    bits = -np.log2(probs[np.arange(len(y)), y]).mean()
    H = -(probs * np.log2(probs)).sum(1).mean()
    freq = np.bincount(y, minlength=4) / len(y)
    return bits, H, freq, q.mean()
for tag, lo, hi, ml in (('near', 1_000_000, 2_000_000, 0.836), ('far', 9_000_000, 10_000_000, 0.839)):
    b, H, f, qm = score(lo, hi)
    print(f'(a) {tag}: theory-only predictor (sieve to 47 + Mertens, no training) scores {b:.4f} bits per word (its own expected value {H:.4f}); best ML {ml}; survivor prime share q = {qm:.3f}')
offs = []
for m in range(8, -1, -1):
    offs += [-6 * m - 1, -6 * m + 1]
S, mm = singular(offs)
n = len(offs)
for tag, lo, hi in (('near', 1_000_000, 2_000_000), ('far', 9_000_000, 10_000_000)):
    cc = np.arange(lo, hi, 6, dtype=np.float64)
    I = np.array([np.sum(np.log(cc) ** -k) for k in range(n + 1)]) / len(cc)
    h = S * I[mm]
    idx = np.arange(1 << n)
    for i in range(n):
        bit = 1 << i
        lo_ = idx[(idx & bit) == 0]
        h[lo_] -= h[lo_ | bit]
    h = np.clip(h, 0, None); h /= h.sum()
    last = idx & 3; past = idx >> 2
    pw = np.bincount(last, weights=h, minlength=4)
    Hw = -np.sum(pw * np.log2(pw))
    pp = np.bincount(past, weights=h, minlength=1 << (n - 2))
    with np.errstate(divide='ignore', invalid='ignore'):
        Hj = -np.sum(np.where(h > 0, h * np.log2(h), 0)); Hp = -np.sum(np.where(pp > 0, pp * np.log2(pp), 0))
    MI = Hw + Hp - Hj
    print(f'(b) {tag}: Hardy-Littlewood word entropy {Hw:.4f} bits; information in the previous 8 words about the next: {MI:.4f} bits per word')
print('ML (README of the run): class frequencies 1.492 / 1.408 bits per word (near / far); previous 8 words: logistic 1.472 / 1.397, boosted trees 1.458 / 1.387')

"""Own Hurwitz class number: H(D) = (h(D0)/u(D0)) * sum_{g | f0} g prod_{l | g} (1 - chi(l)/l), -D = D0 f0^2,
h(D0) from PARI quadclassunit (fundamental discriminants only; Buchmann, not qfbhclassno), u = 3, 2, 1."""
from fractions import Fraction as Fr
import cypari2
pari = cypari2.Pari()
pari.allocatemem(200 * 10 ** 6)
_cl = {}

def hfund(d0):
    v = _cl.get(d0)
    if v is None:
        if d0 == -3 or d0 == -4:
            v = 1
        else:
            v = int(pari.quadclassunit(d0)[0])
        _cl[d0] = v
    return v

def kron(a, l):
    return int(pari.kronecker(a, l))

def H_own(D, fac=None):
    """D > 0, D = 0,3 mod 4. fac: factorisation dict of D (optional)."""
    if fac is None:
        f = pari.factor(D)
        fac = {int(f[0][i]): int(f[1][i]) for i in range(len(f[0]))}
    # -D = D0 f0^2: odd primes: f0 exponent floor(e/2); D0 squarefree part times 4/8 adjustments
    core = -1
    f0 = {}
    for l, e in fac.items():
        if l == 2:
            continue
        if e % 2:
            core *= l
        if e // 2:
            f0[l] = e // 2
    e2 = fac.get(2, 0)
    # -D = -(2^e2) * odd; odd part sign included in core (core = -prod l^(e mod 2))
    # write -D = 2^e2 * core * (odd square)
    if e2 % 2 == 0:
        b = e2 // 2
        if core % 4 == 1:
            d0 = core; a2 = b
        else:
            d0 = 4 * core; a2 = b - 1
    else:
        b = (e2 - 1) // 2
        d0 = 8 * core; a2 = b - 1
    assert a2 >= 0
    if a2:
        f0[2] = a2
    u = 3 if d0 == -3 else (2 if d0 == -4 else 1)
    s = Fr(1)
    for l, a in f0.items():
        chi = kron(d0, l)
        loc = Fr(1)
        for j in range(1, a + 1):
            loc += l ** j * (1 - Fr(chi, l))
        s *= loc
    return Fr(hfund(d0), u) * s

if __name__ == '__main__':
    import random, time
    random.seed(1)
    bad = 0
    Ds = [d for d in range(3, 4000) if d % 4 in (0, 3)] + [4 * random.randrange(10 ** 8, 10 ** 9) + random.choice((0, 3)) for _ in range(300)]
    t = time.time()
    for D in Ds:
        a = H_own(D); b = Fr(str(pari.qfbhclassno(D)))
        if a != b:
            bad += 1; print('diff', D, a, b)
    print(f'H_own vs qfbhclassno: {len(Ds)} values, {bad} differences, {time.time()-t:.1f}s')

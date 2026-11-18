import json, math
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nint

mp.dps = 160
data = json.load(open('classpoly.json'))


def forms(D):
    out = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if a == c and b < 0:
                continue
            if math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            out.append((a, b, c))
        a += 1
    return out


def jinv(tau):
    q = exp(2j * pi * tau)
    e4 = mpc(1)
    e6 = mpc(1)
    qn = q
    n = 1
    while abs(qn) > mpf(10) ** (-mp.dps - 5):
        t = qn / (1 - qn)
        e4 += 240 * n**3 * t
        e6 -= 504 * n**5 * t
        n += 1
        qn *= q
    return 1728 * e4**3 / (e4**3 - e6**2)


bad = 0
for D, h, co in data:
    fs = forms(D)
    assert len(fs) == h, (D, h, len(fs))
    poly = [mpc(1)]
    for a, b, c in fs:
        j = jinv((-b + sqrt(mpf(D))) / (2 * a))
        new = [mpc(0)] * (len(poly) + 1)
        for i, v in enumerate(poly):
            new[i] += v
            new[i + 1] -= v * j
        poly = new
    got = [int(nint(v.real)) for v in poly]
    err = max(abs(v.real - g) + abs(v.imag) for v, g in zip(poly, got))
    if got != [int(c) for c in co] or err > mpf(10) ** -20:
        bad += 1
        print('MISMATCH', D)
print('class polynomials recomputed from the j-function with mpmath (', mp.dps, 'digits ):', len(data), 'checked,', bad, 'differ')

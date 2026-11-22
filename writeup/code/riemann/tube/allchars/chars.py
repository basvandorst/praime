"""Primitive Dirichlet characters of conductor q <= 13, one per conjugate pair."""
import cmath, math


def prim_root(q):
    for g in range(2, q):
        if len({pow(g, j, q) for j in range(q - 1)}) == q - 1:
            return g


def chars():
    out = []   # (name, q, values list length q (complex), parity a, real?)
    out.append(('chi3', 3, [0, 1, -1], 1, True))
    out.append(('chi4', 4, [0, 1, 0, -1], 1, True))
    out.append(('chi8', 8, [0, 1, 0, -1, 0, -1, 0, 1], 0, True))
    out.append(('chi-8', 8, [0, 1, 0, 1, 0, -1, 0, -1], 1, True))
    out.append(('chi12', 12, [0, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 1], 0, True))
    # q = 9: generator 2, order 6; primitive k = 1, 2 (5, 4 conjugate)
    for k in (1, 2):
        v = [0] * 9
        for j in range(6):
            v[pow(2, j, 9)] = cmath.exp(2j * math.pi * k * j / 6)
        out.append(('q9 ord%d' % (6 // math.gcd(6, k)), 9, v, (1 if k % 2 else 0), False))
    for q in (5, 7, 11, 13):
        g = prim_root(q); m = q - 1
        for k in range(1, m // 2 + 1):          # k and m-k conjugate
            v = [0] * q
            for j in range(m):
                v[pow(g, j, q)] = cmath.exp(2j * math.pi * k * j / m)
            order = m // math.gcd(m, k)
            real = (order <= 2)
            if real:
                v = [complex(round(z.real)) for z in v]
            par = 0 if (k % 2 == 0) else 1      # chi(-1) = chi(g^(m/2)) = (-1)^k
            out.append(('q%d ord%d k%d' % (q, order, k), q, v, par, real))
    # clean tiny imaginary parts
    res = []
    for name, q, v, a, real in out:
        v = [complex(z) for z in v]
        v = [complex(round(z.real, 15), round(z.imag, 15)) for z in v]
        if real:
            v = [z.real for z in v]
        res.append((name, q, v, a, real))
    return res


if __name__ == '__main__':
    for c in chars():
        print(c[0], c[1], 'odd' if c[3] else 'even', 'real' if c[4] else 'complex')
    print(len(chars()))

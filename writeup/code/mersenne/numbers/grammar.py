"""Grammar of integer maps x -> f(x) (optionally using the index n), vectorised over int64 arrays.

Convention: -1 marks an invalid value; every primitive maps invalid -> invalid and any result
outside [1, VMAX] -> invalid. VMAX = 2^53 keeps float arithmetic exact enough.
Each primitive is (name, func, tags); tags: 'core2' = used in depth-2 compositions,
'core3' = used in depth-3 compositions, 'n' = uses the index n.
"""
import numpy as np, math
from fractions import Fraction

VMAX = 2 ** 53
I64 = np.int64


def clean(y, x):
    y = np.asarray(y)
    if y.dtype != I64:
        yf = np.asarray(y, dtype=np.float64)
        bad = ~np.isfinite(yf) | (yf < 1) | (yf > VMAX)
        y = np.where(bad, -1, yf).astype(I64)
        return np.where(x < 1, -1, y)
    return np.where((x < 1) | (y < 1) | (y > VMAX), -1, y)


def bitlen(x):
    return np.frexp(np.maximum(x, 1).astype(np.float64))[1].astype(I64)


def popcount(x):
    return np.bitwise_count(np.maximum(x, 0).astype(np.uint64)).astype(I64)


_REVTAB = np.array([int('{:08b}'.format(i)[::-1], 2) for i in range(256)], np.uint8)


def rev64(x):
    u = np.ascontiguousarray(np.maximum(x, 0).astype(np.uint64))
    b = u.view(np.uint8).reshape(u.shape + (8,))
    r = _REVTAB[b][..., ::-1].copy()
    return r.view(np.uint64).reshape(u.shape)


def revbin(x):
    L = bitlen(x).astype(np.uint64)
    return (rev64(x) >> (np.uint64(64) - L)).astype(I64)


def ndig(x):
    d = np.floor(np.log10(np.maximum(x, 1).astype(np.float64))).astype(I64) + 1
    p = 10 ** np.minimum(d, 18)
    d = np.where(np.maximum(x, 1) >= p, d + 1, d)
    p = 10 ** np.minimum(d - 1, 18)
    d = np.where(np.maximum(x, 1) < p, d - 1, d)
    return d


def digit_op(x, kind):
    t = np.maximum(x, 0).copy(); r = np.zeros_like(t) if kind != 'prod' else np.ones_like(t)
    for _ in range(17):
        m = t > 0
        d = t % 10
        if kind == 'rev':
            r = np.where(m, r * 10 + d, r)
        elif kind == 'sum':
            r = r + d
        elif kind == 'sq':
            r = r + d * d
        elif kind == 'prod':
            r = np.where(m, r * d, r)
        t = t // 10
        if not m.any():
            break
    return r


def isqrt_v(x):
    s = np.floor(np.sqrt(np.maximum(x, 0).astype(np.float64))).astype(I64)
    s = np.where(s * s > x, s - 1, s)
    s = np.where((s + 1) * (s + 1) <= x, s + 1, s)
    return s


def next_in(S):
    S = np.asarray(S, I64)
    def f(x, n):
        i = np.searchsorted(S, x, side='right')
        ok = i < len(S)
        return np.where(ok, S[np.minimum(i, len(S) - 1)], -1)
    return f


def gosper(x):
    c = x & -x
    r = x + c
    return (((r ^ x) >> 2) // np.maximum(c, 1)) | r


def next_popcount_k(x, k):
    y = x + 1
    for _ in range(64):
        pc = popcount(y)
        hi = pc > k
        if not hi.any():
            break
        y = np.where(hi, (y | (y - 1)) + 1, y)   # clear lowest run of ones, carry
    for _ in range(64):
        pc = popcount(y)
        lo = pc < k
        if not lo.any():
            break
        y = np.where(lo, y | (y + 1), y)          # set lowest zero bit
    return y


def bin_palindromes(maxbits=40):
    out = []
    for L in range(1, maxbits + 1):
        h = (L + 1) // 2
        a = np.arange(1 << (h - 1), 1 << h, dtype=I64)
        low = revbin(a)
        if L % 2 == 0:
            v = (a << h) | low
        else:
            v = (a << (h - 1)) | (low & ((1 << (h - 1)) - 1)) if h > 1 else a
        out.append(v)
    return np.unique(np.concatenate(out))


def dec_palindromes(maxd=12):
    out = [np.arange(1, 10, dtype=I64)]
    for D in range(2, maxd + 1):
        h = (D + 1) // 2
        a = np.arange(10 ** (h - 1), 10 ** h, dtype=I64)
        r = digit_op(a, 'rev')
        if D % 2 == 0:
            v = a * 10 ** h + r
        else:
            v = a * 10 ** (h - 1) + r % 10 ** (h - 1)
        out.append(v)
    return np.unique(np.concatenate(out))


def build():
    P = []
    def add(name, fn, *tags):
        P.append((name, fn, set(tags)))

    # --- additive / linear -----------------------------------------------------------
    for c in range(-30, 31):
        if c:
            add(f'x{c:+d}', (lambda c: lambda x, n: x + c)(c), *(('core2',) if abs(c) <= 12 else ()),
                *(('core3',) if c in (-2, -1, 1, 2) else ()))
    for a in (2, 3, 4, 5, 6, 7, 8, 10, 16):
        for d in range(-15, 16):
            tg = ('core2',) if (abs(d) <= 3 and a in (2, 3, 4, 6, 8, 10)) else ()
            tg3 = ('core3',) if (a in (2, 4) and d in (-1, 0, 1)) else ()
            add(f'{a}x{d:+d}', (lambda a, d: lambda x, n: a * x + d)(a, d), *tg, *tg3)
    for c in range(2, 11):
        add(f'floor(x/{c})', (lambda c: lambda x, n: x // c)(c), 'core2', *(('core3',) if c == 2 else ()))
        add(f'round(x/{c})', (lambda c: lambda x, n: (2 * x + c) // (2 * c))(c), 'core2')
    # --- multiplicative by constants -------------------------------------------------------
    R = {}
    for q in range(1, 7):
        for p in range(q + 1, 4 * q + 1):
            fr = Fraction(p, q)
            if fr.denominator == q and fr not in R:
                R[fr] = f'{p}/{q}' if q > 1 else f'{p}'
    consts = {'e': math.e, 'pi': math.pi, 'phi': (1 + 5 ** .5) / 2, 'phi^2': ((1 + 5 ** .5) / 2) ** 2,
              'sqrt2': 2 ** .5, 'sqrt3': 3 ** .5, 'sqrt5': 5 ** .5, '2^(1/3)': 2 ** (1 / 3),
              '2^(1/4)': 2 ** .25, '2^(1/5)': 2 ** .2, '2^(3/2)': 2 ** 1.5, 'e^gamma': math.exp(0.5772156649015329),
              'ln10': math.log(10), 'e/2': math.e / 2, 'pi/2': math.pi / 2, '3/e': 3 / math.e,
              '2^(1/e^gamma)': 2 ** math.exp(-0.5772156649015329), 'e^(1/e)': math.e ** (1 / math.e),
              'ln2*2': 2 * math.log(2)}
    mults = [(nm, fr) for fr, nm in R.items() if fr != int(fr)] + list(consts.items())
    mults += [('1/' + nm, 1 / v) for nm, v in mults]
    core_m = {'3/2', '4/3', '5/3', '5/2', 'e', 'pi', 'phi', 'sqrt2', '2^(1/e^gamma)', '1/3/2', '1/phi', '1/sqrt2', '1/e'}
    for nm, v in mults:
        for mode, op in (('floor', np.floor), ('round', np.rint), ('ceil', np.ceil)):
            tg = ('core2',) if (nm in core_m and mode != 'ceil') else ()
            tg3 = ('core3',) if (nm in ('3/2', 'phi', '2^(1/e^gamma)') and mode == 'round') else ()
            if isinstance(v, Fraction):   # exact integer arithmetic for rational multipliers
                add(f'{mode}(x*{nm})', (lambda p, q, mode: lambda x, n: _ratmul(x, p, q, mode))(v.numerator, v.denominator, mode), *tg, *tg3)
            else:
                add(f'{mode}(x*{nm})', (lambda v, op: lambda x, n: op(x.astype(np.float64) * v))(v, op), *tg, *tg3)
    # --- powers ----------------------------------------------------------------------------
    for nm, e_ in (('1.1', 1.1), ('1.2', 1.2), ('5/4', 1.25), ('4/3', 4 / 3), ('3/2', 1.5), ('5/3', 5 / 3),
                   ('phi', (1 + 5 ** .5) / 2), ('sqrt2', 2 ** .5), ('e', math.e), ('2.5', 2.5), ('3', 3.0),
                   ('1/2', .5), ('2/3', 2 / 3), ('3/4', .75)):
        for mode, op in (('floor', np.floor), ('round', np.rint)):
            add(f'{mode}(x^{nm})', (lambda e_, op: lambda x, n: op(np.power(x.astype(np.float64), e_)))(e_, op),
                *(('core2',) if nm in ('3/2', '4/3', '1/2') else ()))
    for c in range(-10, 11):
        add(f'x^2{c:+d}', (lambda c: lambda x, n: np.where(x < 2 ** 26, x * x + c, -1))(c), *(('core2',) if abs(c) <= 2 else ()))
    for c in (-2, -1, 0, 1, 2):
        add(f'x^2+x{c:+d}', (lambda c: lambda x, n: np.where(x < 2 ** 26, x * x + x + c, -1))(c))
        add(f'x(x+1)/2{c:+d}', (lambda c: lambda x, n: np.where(x < 2 ** 26, x * (x + 1) // 2 + c, -1))(c))
        add(f'x^3{c:+d}', (lambda c: lambda x, n: np.where(x < 2 ** 17, x ** 3 + c, -1))(c))
    for c in (-3, -2, -1, 0, 1, 2, 3):
        add(f'2^x{c:+d}', (lambda c: lambda x, n: np.where(x < 53, np.left_shift(1, np.minimum(x, 52)) + c, -1))(c),
            *(('core2',) if c in (-1, 1) else ()))
    add('3^x', lambda x, n: np.where(x < 33, 3 ** np.minimum(x, 32), -1))
    add('x*2^x', lambda x, n: np.where(x < 45, x << np.minimum(x, 44), -1))
    add('2^x-x', lambda x, n: np.where(x < 53, (1 << np.minimum(x, 52)) - x, -1))
    add('fib(x)', lambda x, n: np.where(x < 75, _FIB[np.minimum(x, 74)], -1))
    add('x!', lambda x, n: np.where(x < 19, _FACT[np.minimum(x, 18)], -1))
    # --- self-referential increments -----------------------------------------------------
    add('x+popcount(x)', lambda x, n: x + popcount(x), 'core2')
    add('x+bitlen(x)', lambda x, n: x + bitlen(x), 'core2')
    add('x+digitsum(x)', lambda x, n: x + digit_op(x, 'sum'), 'core2')
    add('x+digitprod(x)', lambda x, n: x + digit_op(x, 'prod'), 'core2')
    add('x+revdec(x)', lambda x, n: x + digit_op(x, 'rev'), 'core2')
    add('x+revbin(x)', lambda x, n: x + revbin(x), 'core2')
    add('x+isqrt(x)', lambda x, n: x + isqrt_v(x), 'core2')
    add('x+floor(ln x)', lambda x, n: x + np.floor(np.log(np.maximum(x, 1))).astype(I64), 'core2')
    add('x+bitlen(x)^2', lambda x, n: x + bitlen(x) ** 2, 'core2')
    add('x+ndig(x)', lambda x, n: x + ndig(x), 'core2')
    add('x*bitlen(x)', lambda x, n: x * bitlen(x), 'core2')
    add('x*ndig(x)', lambda x, n: x * ndig(x), 'core2')
    add('x+x>>1|1', lambda x, n: (x + (x >> 1)) | 1)
    # --- index-based (use n) ----------------------------------------------------------------
    for c in (1, 2, 3, 4, 6, 10):
        add(f'x+{c}n', (lambda c: lambda x, n: x + c * n)(c), 'n', 'core2')
        add(f'x+{c}n^2', (lambda c: lambda x, n: x + c * n * n)(c), 'n')
    for c in (2, 3, 4):
        add(f'{c}x+n', (lambda c: lambda x, n: c * x + n)(c), 'n', 'core2')
        add(f'{c}x-n', (lambda c: lambda x, n: c * x - n)(c), 'n', 'core2')
    add('x*n', lambda x, n: x * n, 'n', 'core2')
    add('x*n+1', lambda x, n: x * n + 1, 'n')
    add('x*n-1', lambda x, n: x * n - 1, 'n')
    add('x+2^n', lambda x, n: x + (1 << np.minimum(n, 52)), 'n', 'core2')
    add('x+fib(n)', lambda x, n: x + _FIB[np.minimum(n, 74)], 'n')
    add('x+n!', lambda x, n: np.where(n < 19, x + _FACT[np.minimum(n, 18)], -1), 'n')
    add('x+prime(n)', lambda x, n: x + _PR[n], 'n')
    # --- bit operations ----------------------------------------------------------------------
    add('revbin', lambda x, n: revbin(x), 'core2', 'core3')
    for k in (1, 2, 3):
        add(f'rotl{k}', (lambda k: lambda x, n: _rotl(x, k))(k), *(('core2',) if k == 1 else ()), *(('core3',) if k == 1 else ()))
        add(f'rotr{k}', (lambda k: lambda x, n: _rotl(x, -k))(k), *(('core2',) if k == 1 else ()), *(('core3',) if k == 1 else ()))
    add('compl', lambda x, n: x ^ ((1 << bitlen(x)) - 1), 'core2', 'core3')
    add('gray', lambda x, n: x ^ (x >> 1), 'core2')
    add('invgray', lambda x, n: _invgray(x), 'core2')
    add('x|x>>1', lambda x, n: x | (x >> 1), 'core2')
    add('x|x<<1', lambda x, n: x | (x << 1), 'core2')
    add('x^x<<1', lambda x, n: x ^ (x << 1), 'core2')
    add('nextsamepop', lambda x, n: gosper(x), 'core2', 'core3')
    add('prepend1', lambda x, n: x + (1 << bitlen(x)), 'core2', 'core3')
    add('prepend10', lambda x, n: x + (2 << bitlen(x)), 'core2')
    add('prepend11', lambda x, n: x + (3 << bitlen(x)), 'core2')
    add('drophigh', lambda x, n: x - (1 << (bitlen(x) - 1)), 'core2')
    add('bincat(x,x)', lambda x, n: np.where(bitlen(x) < 26, (x << bitlen(x)) | x, -1), 'core2')
    add('binpal(x)', lambda x, n: np.where(bitlen(x) < 26, (x << bitlen(x)) | revbin(x), -1), 'core2')
    add('binpal_odd(x)', lambda x, n: np.where(bitlen(x) < 26, (x << (bitlen(x) - 1)) | (revbin(x) & ((1 << (bitlen(x) - 1)) - 1)), -1), 'core2')
    for k in range(0, 28):
        add(f'setbit{k}', (lambda k: lambda x, n: x | (1 << k))(k), *(('core2',) if k <= 4 else ()))
        add(f'clrbit{k}', (lambda k: lambda x, n: x & ~(1 << k))(k), *(('core2',) if k <= 4 else ()))
        add(f'flipbit{k}', (lambda k: lambda x, n: x ^ (1 << k))(k))
    add('popcount', lambda x, n: popcount(x), 'core2')
    add('bitlen', lambda x, n: bitlen(x), 'core2')
    add('2^bitlen-1', lambda x, n: (1 << bitlen(x)) - 1, 'core2')
    add('2^bitlen+1', lambda x, n: (1 << bitlen(x)) + 1, 'core2')
    # --- decimal ---------------------------------------------------------------------------
    add('revdec', lambda x, n: digit_op(x, 'rev'), 'core2', 'core3')
    add('digitsum', lambda x, n: digit_op(x, 'sum'), 'core2')
    add('digitprod', lambda x, n: digit_op(x, 'prod'))
    add('digitsqsum', lambda x, n: digit_op(x, 'sq'))
    add('deccat(x,x)', lambda x, n: np.where(x < 10 ** 8, x * 10 ** ndig(x) + x, -1), 'core2')
    add('deccat(x,revdec x)', lambda x, n: np.where(x < 10 ** 8, x * 10 ** ndig(x) + digit_op(x, 'rev'), -1), 'core2')
    add('decrotl', lambda x, n: (x % 10 ** (ndig(x) - 1)) * 10 + x // 10 ** (ndig(x) - 1), 'core2')
    add('decrotr', lambda x, n: (x % 10) * 10 ** (ndig(x) - 1) + x // 10, 'core2')
    for d in range(1, 10):
        add(f'prepend_dec{d}', (lambda d: lambda x, n: d * 10 ** ndig(x) + x)(d), *(('core2',) if d <= 3 else ()))
    # --- next number > x with property P -------------------------------------------------------
    for m in list(range(2, 13)) + [16, 18, 20, 24, 30, 60]:
        for r in range(m):
            tg = ('core2',) if m <= 6 else ()
            tg3 = ('core3',) if (m, r) in ((4, 1), (4, 3), (6, 1), (6, 5)) else ()
            add(f'next≡{r}(mod {m})', (lambda m, r: lambda x, n: x + 1 + ((r - x - 1) % m))(m, r), *tg, *tg3)
    add('next binpal', next_in(_BPAL), 'core2', 'core3')
    add('next decpal', next_in(_DPAL), 'core2', 'core3')
    for k in range(1, 9):
        add(f'next pop={k}', (lambda k: lambda x, n: next_popcount_k(x, k))(k), *(('core2',) if k <= 4 else ()))
    add('next pop=pop(x)+1', lambda x, n: _next_pop_rel(x, 1), 'core2')
    add('next pop=pop(x)-1', lambda x, n: _next_pop_rel(x, -1), 'core2')
    add('next pop<=2', next_in(_POP2), 'core2', 'core3')
    for c in range(-3, 4):
        add(f'next 2^k{c:+d}', next_in(np.unique([(1 << k) + c for k in range(2, 54) if (1 << k) + c > 0])),
            'core2', *(('core3',) if c in (-1, 1) else ()))
    add('next 2^k±1', next_in(_PM1), 'core2', 'core3')
    add('next 2^k±small(<=3)', next_in(_PM3), 'core2')
    add('next 4^k±3', next_in(np.unique([v for k in range(1, 27) for v in ((4 ** k) - 3, (4 ** k) + 3)])), 'core2')
    add('next square', lambda x, n: (isqrt_v(x) + 1) ** 2, 'core2')
    add('next square+1', lambda x, n: (isqrt_v(x - 1) + 1) ** 2 + 1, 'core2')
    add('next square-1', lambda x, n: (isqrt_v(x + 1) + 1) ** 2 - 1, 'core2')
    add('next triangular', next_in(_TRI), 'core2')
    add('next fibonacci', next_in(_FIB[2:]), 'core2')
    add('next 3^k', next_in([3 ** k for k in range(0, 34)]), 'core2')
    add('next 3^k±2', next_in(np.unique([3 ** k + s for k in range(1, 34) for s in (-2, 2)])), 'core2')
    add('next 10^k±small', next_in(np.unique([10 ** k + s for k in range(1, 16) for s in range(-3, 4)])))
    add('next k*2^m+1 (k<=7 odd)', next_in(_PROTH), 'core2')
    add('next k*2^m-1 (k<=7 odd)', next_in(_RIESEL), 'core2')
    return P


def _ratmul(x, p, q, mode):
    x = np.where(x < (1 << 57) // p, x, -1)
    t = x * p
    fl = t // q
    if mode == 'floor':
        y = fl
    elif mode == 'ceil':
        y = -((-t) // q)
    else:  # round half to even, like numpy rint and Python round(Fraction)
        rem = t - fl * q
        y = np.where(2 * rem > q, fl + 1, np.where(2 * rem == q, fl + (fl & 1), fl))
    return np.where(x < 0, -1, y)


def _rotl(x, k):
    L = bitlen(x)
    k = np.full_like(x, k) % L
    mask = (1 << L) - 1
    return ((x << k) & mask) | (x >> (L - k))


def _invgray(x):
    y = x.copy()
    for s in (1, 2, 4, 8, 16, 32):
        y = y ^ (y >> s)
    return y


def _next_pop_rel(x, d):
    pc = popcount(x) + d
    out = np.full_like(x, -1)
    for k in range(1, 54):
        m = pc == k
        if m.any():
            out[m] = next_popcount_k(x[m], k)
    return out


_FIB = np.array([0, 1] + [0] * 73, dtype=object)
for _i in range(2, 75):
    _FIB[_i] = _FIB[_i - 1] + _FIB[_i - 2]
_FIB = np.array([int(v) for v in _FIB], dtype=I64)
_FACT = np.array([math.factorial(i) for i in range(19)], dtype=I64)
_TRI = np.array([k * (k + 1) // 2 for k in range(1, 1 << 22)], dtype=I64)
_PM1 = np.unique([(1 << k) + s for k in range(1, 54) for s in (-1, 1)])
_PM3 = np.unique([(1 << k) + s for k in range(2, 54) for s in range(-3, 4)])
_POP2 = np.unique([(1 << i) | (1 << j) for i in range(54) for j in range(i, 54)])
_PROTH = np.unique([k * (1 << m) + 1 for k in (1, 3, 5, 7) for m in range(1, 51)])
_RIESEL = np.unique([k * (1 << m) - 1 for k in (1, 3, 5, 7) for m in range(1, 51)])
_PR = np.array([0, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
                101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
                199, 211, 223, 227, 229, 233, 239, 241, 251], dtype=I64)
_BPAL = bin_palindromes(40)
_DPAL = dec_palindromes(12)


def apply(fn, x, n):
    with np.errstate(all='ignore'):
        return clean(fn(np.where(x < 1, 1, x), n), x)


if __name__ == '__main__':
    P = build()
    print(len(P), 'primitives;', sum('core2' in t for _, _, t in P), 'core2;', sum('core3' in t for _, _, t in P), 'core3')
    x = np.array([1, 2, 3, 5, 7, 13, 61, 127, 521, 9689, 136279841], I64); n = np.arange(1, len(x) + 1)
    for nm, fn, _ in P:
        if nm in ('revbin', 'rotl1', 'rotr1', 'compl', 'nextsamepop', 'next binpal', 'next decpal', 'next pop=3',
                  'revdec', 'decrotl', 'next≡1(mod 4)', 'binpal_odd(x)', 'invgray', 'next pop=pop(x)+1', 'next square+1',
                  'round(x*3/2)', 'x+prime(n)', 'digitsum', 'next k*2^m+1 (k<=7 odd)'):
            print(nm, list(apply(fn, x, n)))
    print('bin pal count', len(_BPAL), _BPAL[:12], 'dec pal', len(_DPAL), _DPAL[:12])

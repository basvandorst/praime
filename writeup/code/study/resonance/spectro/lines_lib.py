"""lines_lib.py: predicted amplitudes of the chirp lines r = a/b on PRIMES (spectro.md, section 2).

  sum_p lambda(p) ln p g(p) e(-2 sqrt(r p/N))  ~  c(r) * MT(r),   MT(r) = (3/2) K1 r^{-1/4} int ... (spectro_lib.main_term_unit)
  c(r) = prod_l L_l(r),  good l (chi(l) != 0), e_a = v_l(a), e_b = v_l(b):
      e_b = 0:  lambda(l^e_a) - chi(l) lambda(l^(e_a+2)) / l
      e_b = 1:  -chi(l) lambda(l) / l
      e_b = 2:  chi(l) / (l (l - 1))
      e_b >= 3: 0
  bad l: L_l = A_l(k), k = v_l(r), measured EXACTLY on integers sieved at l alone (bad_factors()); for k = 0 it is c_l.
  The product over l !| ab is the single-prime constant c (round 153); so c(r) = c * prod_{l | ab} L_l(r) / L_l(1)."""
import math, os, json
import numpy as np
import spectro_lib as S

def vval(x, l):
    k = 0
    while x % l == 0: x //= l; k += 1
    return k

def factor(n):
    out = {}; q = 2
    while q * q <= n:
        while n % q == 0: out[q] = out.get(q, 0) + 1; n //= q
        q += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out

_LOC = {}
def loc(label):
    if label not in _LOC: _LOC[label] = S.local_lams(label, 2000)
    return _LOC[label]

def lam_pk(label, l, k):
    lp, chi = loc(label)[l]
    return S.lam_pow_table(lp, chi, max(k, 2) + 2)[k]

def L_good(label, l, ea, eb):
    lp, chi = loc(label)[l]
    if eb == 0: return lam_pk(label, l, ea) - chi * lam_pk(label, l, ea + 2) / l
    if eb == 1: return -chi * lp / l
    if eb == 2: return chi / (l * (l - 1.0))
    return 0.0

def bad_factor_file(label):
    return os.path.join(S.HERE, 'badfac_%s.json' % label.replace(':', '_'))

def c_line(label, a, b, bad=None):
    """c(r) for r = a/b in lowest terms (real), using bad-prime factors from bad (dict 'l|k' -> value) or the file."""
    I = S.label_info(label); c = I['c']
    if bad is None:
        fn = bad_factor_file(label); bad = json.load(open(fn)) if os.path.exists(fn) else {}
    out = c
    for l in sorted(set(factor(a)) | set(factor(b))):
        ea, eb = vval(a, l), vval(b, l)
        lp, chi = loc(label)[l]
        if chi != 0.0:
            out *= L_good(label, l, ea, eb) / L_good(label, l, 0, 0)
        else:
            key = '%d|%d' % (l, ea - eb)
            if key not in bad: return float('nan')
            out *= bad[key] / bad['%d|0' % l]
    return out

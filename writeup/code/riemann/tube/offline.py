"""What an off-line zero would do to the primes, read through the tube of fact 2048.
G(y) - ln 2pi + T(y) = sum_rho y^rho/(rho(1-rho)); on the line every term has size sqrt(y)/|rho|^2 and the sum is at most
C sqrt(y), C = sum 1/|rho|^2 = 2 + gamma - ln 4 pi. A zero rho0 = beta + i t (beta > 1/2) adds, with its conjugate, a term
of amplitude A(y) = 2 y^beta/|rho0 (1 - rho0)| (its partners 1 - rho0 give y^(1-beta), smaller). Once A(y) > 2 C sqrt(y)
the tube must break at the term's peaks (the rest is at most C sqrt(y) in size, even if all other zeros are on the
line; we take C for the rest, slightly generous). So the tube is certainly broken beyond
  y*(beta, t) = (C |rho0 (1 - rho0)|)^(1/(beta - 1/2)),
and conversely the primes up to Y exclude, by this argument, only zeros with y*(beta, t) <= Y.
usage: offline.py"""
import json, math
C = 2 + 0.5772156649015329 - math.log(4 * math.pi)


def ystar_log10(beta, t):
    rho = complex(beta, t); m = abs(rho * (1 - rho))
    return math.log10(C * m) / (beta - 0.5)


rows = []
for t in (14.13, 1e3, 1e6, 1e9, 3e12):
    rows.append({'t': t, 'log10_ystar': {str(b): ystar_log10(b, t) for b in (0.51, 0.55, 0.6, 0.75, 0.9, 1.0)}})
# the largest height t that the primes up to Y exclude for a given beta: C |rho(1-rho)| <= Y^(beta - 1/2)
excl = {}
for Y in (1e8, 1e18, 1e29):
    excl[str(Y)] = {str(b): math.sqrt(max(Y ** (b - 0.5) / C, 0.0)) for b in (0.6, 0.75, 0.9, 1.0)}
out = {'C': C, 'rows': rows, 'max_height_excluded_by_primes_up_to_Y': excl,
       'note': 'violation at every peak beyond y*, proved if rho0 quadruple is the only off-line one; not a statement about detection below y*'}
print('C =', C)
for r in rows:
    print('t = %-8g' % r['t'], ' '.join('b=%s: 10^%.0f' % (b, v) for b, v in r['log10_ystar'].items()))
for Y, d in excl.items():
    print('primes up to', Y, 'exclude zeros below height', {b: round(v, 1) for b, v in d.items()})
import os
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'offline.json'), 'w'), indent=1)

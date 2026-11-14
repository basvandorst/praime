"""The information law of the explicit formula: digits of the zeros buy height (RIEMANN.md, section 20).

Every run is the certified decoder lpchain.py (section 14) on a stretch [S, S + L) of odd numbers, the primes below S given.
Bits per sample: b = log2(ln x / (delta/K0)), delta the rigorous budget of cert.py at the middle x of the stretch.

  info_law.py budget x d Omega               b at x for zeros rounded to d decimals (d = 28: the full strings)
  info_law.py digits d Omega S [limit [out]] zeros rounded to d decimals; the smallest even beta whose tail is below a
                                             tenth of the precision term, a = 1.2 beta / (Omega pi); stretch of 200
  info_law.py check-ambiguous d out_ambiguous.json   info_law_check.py on an 'ambiguous' stop, with the same zeros
  info_law.py window Omega S L frac [limit]  28-digit zeros, a = frac * x, beta = Omega pi a / 1.2
  info_law.py bisect-digits d lo hi steps S  smallest Omega that proves the stretch (lo assumed to fail)
  info_law.py bisect-window S L frac lo hi steps
  info_law.py stages X c frac Tmin           stages of a chain from 3 with the height max(Tmin, c sqrt(x)) at each stage top

The height at x is Omega * pi * x. A failed run ends either 'ambiguous' (two 0/1 patterns fit every sample: the information
limit, checked by info_law_check.py) or 'search limit' (the branch and bound ran out of nodes).
"""
import sys, os, math, json, subprocess, tempfile, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cert import sample
from kbm import K0
ZFULL = os.path.join(HERE, '../../../data/riemann/zeros_hp2419.txt')
PY = sys.executable

def zeros(d):
    """(file, eps): the 28-digit zeros rounded to d decimals, written to a temporary file."""
    if d >= 28:
        return ZFULL, 2e-16
    S = [l.strip() for l in open(ZFULL) if l.strip()]
    f = os.path.join(tempfile.gettempdir(), f'info_law_zeros_{d}.txt')
    if not os.path.exists(f):
        open(f, 'w').write('\n'.join(format(float(s), f'.{d}f') for s in S) + '\n')
    return f, 0.5 * 10.0 ** (-d) + 2e-16

def bits(x, G, T, a, beta, eps):
    z, dl, info = sample(x, G, T, a, beta, 3, eps)
    rel = dl / K0(beta, 3)
    return math.log2(math.log(x) / rel), rel, info

def lp(zf, eps, end, a, beta, k, S, limit, out=None):
    env = dict(os.environ, ZFILE=os.path.abspath(zf), EPS=repr(eps), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    opts = dict(truth_below=S, start=S, maxbits=40, limit=limit)
    if out: opts['out'] = out
    p = subprocess.run([PY, os.path.join(HERE, 'lpchain.py'), str(end), json.dumps([[end, round(a, 3), round(beta, 4), 3, 1.2, k]]), json.dumps(opts)],
                       capture_output=True, text=True, env=env)
    lines = p.stdout.strip().splitlines()
    return (lines[-2], lines[-1]) if len(lines) >= 2 else (p.stdout[-300:], p.stderr[-300:])

def verdict(res, summ):
    if "'stop': None" in res and 'wrong: []' in summ: return 'ok'
    return 'ambiguous' if 'ambiguous' in res else 'search' if 'search' in res else 'other'

def digits(d, Om, S, limit=1500, out=None):
    zf, eps = zeros(d)
    G = np.array([l.strip() for l in open(zf) if l.strip()], dtype=np.longdouble)
    x = S + 100.0; T = Om * math.pi * x
    for beta in range(8, 46, 2):
        a = 1.2 * beta / (Om * math.pi)
        b, rel, info = bits(x, G, T, a, beta, eps)
        if info['tail'] <= 0.1 * info['prec'] or beta == 44:
            break
    t0 = time.time(); res, summ = lp(zf, eps, S + 200, a, beta, int(round(a / 2)), S, limit, out)
    return dict(d=d, Omega=Om, S=S, beta=beta, a=round(a, 1), bits=round(b, 2), Om_b=round(Om * b, 3), verdict=verdict(res, summ),
                result=res, summary=summ, secs=round(time.time() - t0))

def window(Om, S, L, frac, limit=1500, out=None):
    G = np.array([l.strip() for l in open(ZFULL) if l.strip()], dtype=np.longdouble)
    x = S + L / 2; a = frac * x; beta = Om * math.pi * a / 1.2; T = Om * math.pi * x
    b, rel, info = bits(x, G, T, a, beta, 2e-16)
    t0 = time.time(); res, summ = lp(ZFULL, 2e-16, S + L, a, beta, int(round(a / 2)), S, limit, out)
    return dict(Omega=Om, S=S, L=L, a=round(a, 1), beta=round(beta, 2), bits=round(b, 2), Om_b=round(Om * b, 3),
                T_over_sqrt_x=round(T / math.sqrt(x), 3), zeros=int(np.sum(G <= Om * math.pi * (S + L))), verdict=verdict(res, summ),
                result=res, summary=summ, secs=round(time.time() - t0))

def bisect(run, lo, hi, steps):
    r = run(hi); print(json.dumps(r), flush=True)
    while r['verdict'] != 'ok':
        lo, hi = hi, hi * 1.3; r = run(hi); print(json.dumps(r), flush=True)
    for _ in range(steps):
        mid = math.sqrt(lo * hi); r = run(mid); print(json.dumps(r), flush=True)
        if r['verdict'] == 'ok': hi = mid
        else: lo = mid
    print(f'RESULT fails at {lo:.5f}, proves at {hi:.5f}')

if __name__ == '__main__':
    cmd, A = sys.argv[1], sys.argv[2:]
    if cmd == 'budget':
        x, d, Om = float(A[0]), int(A[1]), float(A[2]); zf, eps = zeros(d)
        G = np.array([l.strip() for l in open(zf) if l.strip()], dtype=np.longdouble)
        for beta in range(8, 46, 2):
            a = 1.2 * beta / (Om * math.pi); b, rel, info = bits(x, G, Om * math.pi * x, a, beta, eps)
            print(f'beta {beta}: a = {a:.1f}, delta/K0 = {rel:.3e} (tail {info["tail"]:.1e}, precision {info["prec"]:.1e}), b = {b:.2f}')
    elif cmd == 'digits':      # digits d Omega S [limit [out.json]]: with out, an 'ambiguous' stop writes out_ambiguous.json
        print(json.dumps(digits(int(A[0]), float(A[1]), int(A[2]), int(A[3]) if len(A) > 3 else 1500, A[4] if len(A) > 4 else None)))
    elif cmd == 'check-ambiguous':   # check-ambiguous d run_ambiguous.json: info_law_check.py with the same zeros and eps
        zf, eps = zeros(int(A[0]))
        env = dict(os.environ, ZFILE=os.path.abspath(zf), EPS=repr(eps))
        subprocess.run([PY, os.path.join(HERE, 'info_law_check.py'), A[1]], env=env)
    elif cmd == 'window':
        print(json.dumps(window(float(A[0]), int(A[1]), int(A[2]), float(A[3]), *(int(v) for v in A[4:5]))))
    elif cmd == 'bisect-digits':
        d, lo, hi, steps, S = int(A[0]), float(A[1]), float(A[2]), int(A[3]), int(A[4])
        bisect(lambda Om: digits(d, Om, S), lo, hi, steps)
    elif cmd == 'bisect-window':
        S, L, frac, lo, hi, steps = int(A[0]), int(A[1]), float(A[2]), float(A[3]), float(A[4]), int(A[5])
        bisect(lambda Om: window(Om, S, L, frac), lo, hi, steps)
    elif cmd == 'stages':
        X, c, frac, Tmin = int(A[0]), float(A[1]), float(A[2]), float(A[3])
        b = [3]
        while b[-1] < X:
            b.append(min(X, max(b[-1] + 2, int(round(b[-1] * 1.25)))))
        st = []
        for s, e in zip(b[:-1], b[1:]):
            a = frac * s; T = max(Tmin, c * math.sqrt(e))
            st.append([e, round(a, 4), round(T * a / (1.2 * e), 4), 3, 1.2, max(1, int(round(a / 2)))])
        print(json.dumps(st))

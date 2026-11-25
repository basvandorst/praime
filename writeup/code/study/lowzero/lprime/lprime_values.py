"""lprime_values.py: L'(rho, chi_D) at the zeros of real primitive Dirichlet L-functions, from a zero file.

For every fundamental discriminant D in the zero file(s) (keys 'D<int>', arrays of zero heights gamma > 0) and every
zero with LO <= gamma <= HI, PARI computes
    d   = L'(1/2 + i gamma, chi_D)      lfun(lfuninit(lfuncreate(D), [1/2, 0, HI + 5], 1), 1/2 + I*gamma, 1)
    v   = |L(1/2 + i gamma, chi_D)|     the same init, derivative order 0 (a check that gamma is a zero)
and stores gamma, log|d|, arg d, |v| in OUTDIR/D<D>.npz (one file per discriminant, written atomically), so a run can
stop at any time and continue. Every run stops before SECONDS (default 540) have passed.

usage:  lprime_values.py OUTDIR ZEROFILE [ZEROFILE ...] [--seconds S] [--lo 100] [--hi 600]
        lprime_values.py OUTDIR --pack PACKFILE     (one npz with every discriminant + a manifest; checks completeness)
Self-contained: numpy and cypari2 only.
"""
import sys, os, time, json, hashlib
import numpy as np


def parse(argv):
    a = {'seconds': 540.0, 'lo': 100.0, 'hi': 600.0, 'pack': None, 'files': []}
    out = argv[1]; i = 2
    while i < len(argv):
        if argv[i] == '--seconds': a['seconds'] = float(argv[i + 1]); i += 2
        elif argv[i] == '--lo': a['lo'] = float(argv[i + 1]); i += 2
        elif argv[i] == '--hi': a['hi'] = float(argv[i + 1]); i += 2
        elif argv[i] == '--pack': a['pack'] = argv[i + 1]; i += 2
        else: a['files'].append(argv[i]); i += 1
    return out, a


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''): h.update(blk)
    return h.hexdigest()


def compute(out, a):
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6); pari('default(realprecision, 38)')
    Z = {}
    for f in a['files']:
        z = np.load(f)
        for k in z.files: Z[k] = np.sort(np.asarray(z[k], float))
    keys = sorted(Z, key=lambda k: (abs(int(k[1:])), int(k[1:])))
    os.makedirs(out, exist_ok=True)
    meta = os.path.join(out, 'manifest.json')
    if not os.path.exists(meta):
        json.dump({'zero_files': [os.path.abspath(f) for f in a['files']], 'zero_sha256': [sha(f) for f in a['files']],
                   'lo': a['lo'], 'hi': a['hi'], 'keys': keys, 'pari': str(pari('version()')), 'realprecision': 38},
                  open(meta, 'w'), indent=1)
    t0 = time.time(); last = 0.0; done = 0
    for k in keys:
        fn = os.path.join(out, k + '.npz')
        if os.path.exists(fn): continue
        if time.time() - t0 + 1.5 * last > a['seconds']:
            break
        t1 = time.time(); D = int(k[1:]); g = Z[k]; b = g[(g >= a['lo']) & (g <= a['hi'])]
        Li = pari('lfuninit(lfuncreate(%d), [1/2, 0, %d], 1)' % (D, int(a['hi']) + 5))
        dre = np.empty(len(b)); dim = np.empty(len(b)); av = np.empty(len(b))
        for j, x in enumerate(b):
            s = pari('1/2 + I*%.17g' % x)
            d = complex(pari.lfun(Li, s, 1)); v = complex(pari.lfun(Li, s))
            dre[j] = d.real; dim[j] = d.imag; av[j] = abs(v)
        ad = np.hypot(dre, dim)
        tmp = fn + '.tmp.npz'
        np.savez(tmp, D=D, gamma=b, logabs=np.log(ad), arg=np.arctan2(dim, dre), absL=av, seconds=time.time() - t1)
        os.replace(tmp, fn)
        last = time.time() - t1; done += 1
        print('D=%d  %d zeros  %.1f s  mean log|L\'| %.4f  max |L|/|L\'| %.1e' % (D, len(b), last, np.log(ad).mean(), (av / ad).max()), flush=True)
    left = [k for k in keys if not os.path.exists(os.path.join(out, k + '.npz'))]
    print('this run: %d discriminants in %.0f s; %d of %d left' % (done, time.time() - t0, len(left), len(keys)), flush=True)


def pack(out, a):
    m = json.load(open(os.path.join(out, 'manifest.json')))
    keys = m['keys']; miss = [k for k in keys if not os.path.exists(os.path.join(out, k + '.npz'))]
    if miss: sys.exit('not complete: %d discriminants missing, e.g. %s' % (len(miss), miss[:5]))
    arrs = {}
    for k in keys:
        z = np.load(os.path.join(out, k + '.npz'))
        for f in ('gamma', 'logabs', 'arg', 'absL'): arrs['%s__%s' % (k, f)] = z[f]
    arrs['manifest'] = np.array(json.dumps(m))
    np.savez_compressed(a['pack'], **arrs)
    print('packed %d discriminants into %s (sha256 %s)' % (len(keys), a['pack'], sha(a['pack'])))


if __name__ == '__main__':
    if len(sys.argv) < 3: sys.exit(__doc__)
    out, a = parse(sys.argv)
    if a['pack']: pack(out, a)
    else: compute(out, a)

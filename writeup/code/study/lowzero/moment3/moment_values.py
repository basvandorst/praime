"""moment_values.py: |L(1/2 + it, chi_D)| on the fixed grid t_k = 100 + k/20, k = 0..18000 (step h = 0.05, [100, 1000]),
one file OUTDIR/absL_<D>.npy (float64, 18001 values) per fundamental discriminant D (D = 1: zeta).  Resumable: a D whose
file exists is skipped; the run stops before the time budget is exceeded (default 540 s, i.e. under 10 minutes).

usage:  moment_values.py OUTDIR D [D ...] [--budget 540]
        moment_values.py OUTDIR --set SETS.json NAME [--budget 540]      (NAME a key of SETS.json holding a list of D)
Method: PARI lfuninit(lfuncreate(D), [1/2, 0, 1001]) and lfun at each grid point, realprecision 19 (checked against
38 digits: max relative difference 2e-16; and against an independent Hurwitz-zeta sum, check_values.py).
Writes OUTDIR/manifest.json (sha256 of every file, grid, PARI version).
"""
import sys, os, json, time, hashlib
import numpy as np
import cypari2

T0, H, NPT = 100.0, 0.05, 18001          # grid: t_k = T0 + k*H, k = 0..NPT-1  (last point 1000.0)


def grid():
    return T0 + H * np.arange(NPT)


def compute(pari, D):
    pari('Li = lfuninit(lfuncreate(%d), [1/2, 0, 1001])' % D)
    v = pari('vector(%d, k, abs(lfun(Li, 1/2 + I*(%d + (k-1)/%d))))' % (NPT, int(T0), int(round(1 / H))))
    return np.array([float(x) for x in v])


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def main(argv):
    argv = list(argv)
    budget = 540.0
    if '--budget' in argv:
        i = argv.index('--budget'); budget = float(argv[i + 1]); del argv[i:i + 2]
    out = argv[1]; os.makedirs(out, exist_ok=True)
    if '--set' in argv:
        i = argv.index('--set'); Ds = json.load(open(argv[i + 1]))[argv[i + 2]]
    else:
        Ds = [int(a) for a in argv[2:]]
    pari = cypari2.Pari(); pari.allocatemem(400 * 10**6); pari.set_real_precision(19)
    start = time.time(); per = 60.0; done = 0
    for D in Ds:
        f = os.path.join(out, 'absL_%d.npy' % D)
        if os.path.exists(f):
            continue
        if time.time() - start + 1.3 * per > budget:
            print('budget reached; rerun to continue'); break
        t1 = time.time()
        v = compute(pari, D)
        np.save(f + '.tmp.npy', v); os.replace(f + '.tmp.npy', f)
        per = max(per * 0.5, time.time() - t1); done += 1
        print('D = %6d  %.1f s  mean|L|^2 = %.4f' % (D, time.time() - t1, float(np.mean(v ** 2))), flush=True)
    man = {'grid': {'t0': T0, 'h': H, 'n': NPT}, 'pari': str(pari.version()), 'realprecision': 19, 'files': {}}
    for fn in sorted(os.listdir(out)):
        if fn.startswith('absL_') and fn.endswith('.npy'):
            man['files'][fn] = sha(os.path.join(out, fn))
    json.dump(man, open(os.path.join(out, 'manifest.json'), 'w'), indent=1)
    left = [D for D in Ds if not os.path.exists(os.path.join(out, 'absL_%d.npy' % D))]
    print('computed %d this run; %d left' % (done, len(left)))


if __name__ == '__main__':
    main(sys.argv)

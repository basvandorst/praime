"""blind zero runs (registered), resumable, in bounded chunks.
usage: bz4_zeros.py --budget 540 [--pred predictions_blind4.json] [--out zeros_blind4] [--only KEY ...]
Runs the manifest of the frozen prediction file in its order (Rankin-Selberg, genus 2, Sym^3; cheapest first). For each key:
skip if <out>/<key>.npz exists; before starting, append to <out>/<key>.attempts; a key with 2 earlier unfinished attempts
is written as REJECTED (timeout) instead of being run again; it is started only if elapsed + 1.5 * est_secs <= budget
(the first key of an invocation always starts). The zero routine and the corrected completeness rule are bz4_lib.zeros (registered).
Each result is written at once (g, info). Rerun the same command until it prints 'ALL DONE'.
Command: (ulimit -v 1500000; timeout 600 nice .venv/bin/python -I bz4_zeros.py --budget 540)"""
import sys, os, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from bz4_lib import zeros
ap = argparse.ArgumentParser(); ap.add_argument('--budget', type=float, default=540); ap.add_argument('--pred', default='predictions_blind4.json')
ap.add_argument('--out', default='zeros_blind4'); ap.add_argument('--only', nargs='*', default=None)
a = ap.parse_args()
t0 = time.time()
pari = cypari2.Pari(); pari.allocatemem(1100 * 10 ** 6)
man = json.load(open(a.pred))['manifest']
os.makedirs(a.out, exist_ok=True)
todo = [m for m in man if (a.only is None or m['key'] in a.only) and not os.path.exists(os.path.join(a.out, m['key'] + '.npz'))]
started = 0
for m in todo:
    key = m['key']; att = os.path.join(a.out, key + '.attempts')
    n_att = len(open(att).read().split()) if os.path.exists(att) else 0
    if n_att >= 2:
        info = dict(spec=m['spec'], T2=m['T2'], status='REJECTED', reason='two unfinished attempts (timeout)')
        np.savez_compressed(os.path.join(a.out, key + '.npz'), g=np.zeros(0), info=json.dumps(info)); print(key, 'REJECTED (timeout)', flush=True); continue
    el = time.time() - t0
    if started and el + 1.5 * m['est_secs'] > a.budget:
        print('budget reached after %.0f s: %d keys left; rerun the same command to continue' % (el, len([x for x in todo if not os.path.exists(os.path.join(a.out, x['key'] + '.npz'))])), flush=True)
        sys.exit(0)
    open(att, 'a').write('%d\n' % int(time.time()))
    started += 1
    g, info = zeros(pari, m['spec'], m['T2'], poles=m['poles'])
    np.savez_compressed(os.path.join(a.out, key + '.npz'), g=g, info=json.dumps(info))
    os.remove(att)
    print('%-34s %-9s n=%5d drift %+.3f r_mean %+.3f (c_exp %+.3f, m0 %d) added %d/%d dup %d  %5.1fs  [elapsed %.0f s]' % (
        key, info['status'], info['n'], info['drift'], info['r_mean'], info['c_exp'], info['m0'], info['added'], info['added_extra'],
        info['duplicates'], info['secs'], time.time() - t0), flush=True)
left = [x for x in man if (a.only is None or x['key'] in a.only) and not os.path.exists(os.path.join(a.out, x['key'] + '.npz'))]
print('ALL DONE' if not left else '%d keys left' % len(left), flush=True)

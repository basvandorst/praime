"""verify the corrected completeness rule on round-151 zero files (read-only): the two characters mod 17 (rejected in round 151)
and all 30 degree-4 sets (Sym^3, genus 2), and on the X7/X5 characters of round 1."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from bz4_lib import build, params, arg_centre, check
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
Z = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'blind', 'zeros_blind')
keys = ['T1_17_5', 'T1_17_7'] + sorted(os.path.basename(f)[:-4] for f in os.listdir(Z) if f.startswith('T4'))
for key in keys:
    z = np.load(os.path.join(Z, key + '.npz')); info = json.loads(str(z['info'])); g = z['g']
    spec = info['spec']
    L = build(pari, spec); Q, k, mu = params(pari, L)
    m0 = int(pari.lfunorderzero(L)); A0, absL = arg_centre(pari, L, k, m0, real=(spec[0] != "char"))
    cexp = 0.5 - m0 / 2 - A0 / math.pi
    g = g[g <= info['T2'] + 1.0]
    drift, rmean, ok = check(g, Q, mu, cexp)
    print('%-28s round-151 status %-9s old c_exp %+.4f | A0 %+.4f m0 %d new c_exp %+.4f r_mean %+.4f drift %+.3f -> %s' % (
        key, info['status'], info['c_exp'], A0, m0, cexp, rmean, drift, 'COMPLETE' if ok else 'REJECTED'), flush=True)

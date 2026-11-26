import json, math
rows = json.load(open('/tmp/claude-1000/referee_martin/rows.json'))
R = {(r[1], r[2], r[3]): r for r in rows}
def comp(l):
    pp = R[(l,1,1)][7]; mp = R[(l,-1,1)][7]
    return (pp+mp)/2, (pp-mp)/2   # Dq-part, DN-part of the limit (eq = +1)
models = {'own M_lim': lambda l,e3,eq: R[(l,e3,eq)][7],
          'own prime-sum M_X': lambda l,e3,eq: R[(l,e3,eq)][8],
          'zero': lambda l,e3,eq: 0.0,
          'W_N only (no W_q term)': lambda l,e3,eq: eq*e3*comp(l)[1],
          'W_q only (no W_N term)': lambda l,e3,eq: eq*comp(l)[0],
          'sign at 3 swapped': lambda l,e3,eq: R[(l,-e3,eq)][7],
          'sign at q swapped': lambda l,e3,eq: R[(l,e3,-eq)][7],
          'both swapped': lambda l,e3,eq: R[(l,-e3,-eq)][7]}
for name, f in models.items():
    c = sum(((r[4]-f(r[1],r[2],r[3]))/r[6])**2 for r in rows)
    print('chi2 %-26s %8.1f  on %d terms' % (name, c, len(rows)))
print('rms(A - M_lim) = %.3f ; rms se = %.3f ; rms(A_fixed) = %.3f ; rms(M_lim) = %.3f' % (
    math.sqrt(sum((r[4]-r[7])**2 for r in rows)/len(rows)), math.sqrt(sum(r[6]**2 for r in rows)/len(rows)),
    math.sqrt(sum(r[5]**2 for r in rows)/len(rows)), math.sqrt(sum(r[7]**2 for r in rows)/len(rows))))
print('window means over the 6 primes l per target (A, A_fixed, M_lim, se_mean):')
for xt in sorted(set(r[0] for r in rows)):
    for e3 in (1,-1):
        for eq in (1,-1):
            sel = [r for r in rows if r[0]==xt and r[2]==e3 and r[3]==eq]
            A = sum(r[4] for r in sel)/len(sel); Af = sum(r[5] for r in sel)/len(sel); M = sum(r[7] for r in sel)/len(sel)
            se = math.sqrt(sum(r[6]**2 for r in sel))/len(sel)
            print('  x=%.1f (%+d,%+d): A %+.3f  fixed %+.3f  M_lim %+.3f  z %+.2f' % (xt, e3, eq, A, Af, M, (A-M)/se))

"""Collect the runs of famcrit3.py over several moduli q and test the two laws:
 uncertainty:  the height needed times the kernel width, T c, is the same for every q (c = ln(1 + q/n) at the largest n);
 conservation: the zeros used, about (q-1) (T/pi) ln(qT/(2 pi e)), equal (1.2 beta/pi) X ln(1.2 beta X/(2 pi e)) for q << X.
Usage: uncertainty_table.py fc3_*.json"""
import sys, json, math
rows = []
for f in sys.argv[1:]:
    rows += json.load(open(f))
X = 10000
by_q = {}
for r in rows: by_q.setdefault(r['q'], []).append(r)
print('q        c         beta  T        T*c    zeros used  certified  correct   max err/(ln n/4)  max delta/(ln n/4)')
for q in sorted(by_q):
    for r in sorted(by_q[q], key=lambda r: r['beta']):
        print(f"{q:<8d} {r['c']:.6f}  {r['beta']:4.0f}  {r['T']:8.2f} {r['Tc']:6.2f} {r['zeros']:10,d}  {r['certified']:5d}/{r['n']:<5d} {r['correct']:5d}/{r['n']:<5d} {r['max_err_over_quarter']:.2e}          {r['max_delta_over_quarter']:.3f}")
print()
print('smallest beta (height) at which every n < 10,000 is certified, and at which every decision is right:')
for q in sorted(by_q):
    rs = sorted(by_q[q], key=lambda r: r['beta'])
    cert = next((r for r in rs if r['certified'] == r['n']), None)
    corr = next((r for r in rs if r['correct'] == r['n']), None)
    inv = lambda b: (1.2 * b / math.pi) * X * math.log(1.2 * b * X / (2 * math.pi * math.e))
    s = f'q = {q}: '
    if cert: s += f"certified from beta = {cert['beta']:.0f}: T = {cert['T']:.2f}, T c = {cert['Tc']:.2f}, zeros {cert['zeros']:,} (conservation formula {inv(cert['beta']):,.0f}); "
    if corr: s += f"all correct from beta = {corr['beta']:.0f}: T = {corr['T']:.2f}, T c = {corr['Tc']:.2f}, zeros {corr['zeros']:,}"
    print(s)

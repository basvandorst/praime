"""Second half of page_check.js: the Lucas-Lehmer residues the page's FFT produced, recomputed from scratch with gmpy2.
  node page_check.js > page_check.json && python3 page_check.py < page_check.json"""
import sys, json, gmpy2
rows = json.loads(sys.stdin.read()); bad = 0
for r in rows:
    q = r['q']; M = (gmpy2.mpz(1) << q) - 1; s = gmpy2.mpz(4)
    for _ in range(q - 2):
        s = (s * s - 2) % M
    ok = str(s & ((1 << 64) - 1)) == r['lo64'] and str(s % ((1 << 61) - 1)) == r['m61']
    bad += not ok
    print(f"q {q}: FFT length {r['n']}, largest rounding error {r['maxErr']:.4f}, residue {'0 (prime)' if s == 0 else 'nonzero'}, page = gmpy2: {ok}")
print(f'{len(rows)} residues, {bad} differ')

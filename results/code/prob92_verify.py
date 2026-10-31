# Problem 92: a Sierpinski and a Riesel number, each with at least 8 distinct (minimal) covering sets.
# Every class is recomputed from k itself: prime p divides k*2^n + s exactly when 2^n = -s/k (mod p).
import sys, ast, numpy as np
from sympy import n_order, isprime
def check(path, s):
  lines = open(path).read().split("\n")
  k = int(lines[0]); tails = ast.literal_eval(lines[2])
  L = 720720
  core = [3, 5, 7, 13]
  def cls(p):                                    # class r mod e with p | k 2^n + 1 iff n = r (mod e), or None
      e = n_order(2, p); assert L % e == 0
      t = (-s * pow(k, -1, p)) % p; x = 1
      for r in range(e):
          if x == t: return e, r
          x = x * 2 % p
      return None
  def covers(ps):
      m = np.zeros(L, bool)
      for p in ps:
          c = cls(p); assert c is not None
          m[c[1]::c[0]] = True
      return m.all()
  assert k % 2 == 1 and k > 10**8
  sets = [tuple(core + [t[0] for t in tail]) for tail in tails]
  for S in sets:
      assert all(isprime(p) for p in S)
      assert covers(S)                                       # covering set
      assert all(not covers([q for q in S if q != p]) for p in S)   # minimal
  assert len(set(sets)) == len(sets) >= 8
  # direct check for small n: k*2^n + s has a prime factor from the first set and is larger than it
  for n in range(1, 3000):
      N = k * 2**n + s
      assert any(N % p == 0 and N > p for p in sets[0])
  print(len(sets), "distinct minimal covering sets; k has", len(str(k)), "digits")
  print("verified")

check("results/data/prob92_k_S.txt", 1)
check("results/data/prob92_k_R.txt", -1)
print("all claims verified")

# a 27-digit Sierpinski number with five distinct minimal covering sets (Gurisatti's four + one disjoint tail)
k5 = 660664662814032650416080443
L5 = 720720
def cls5(p):
    e = n_order(2, p); t = (-pow(k5, -1, p)) % p; x = 1
    for r in range(e):
        if x == t: return e, r
        x = x * 2 % p
def covers5(ps):
    m = np.zeros(L5, bool)
    for p in ps:
        c = cls5(p)
        if c is None: return False
        m[c[1]::c[0]] = True
    return m.all()
sets5 = [(3, 5, 7, 13, 17, 241), (3, 5, 7, 13, 97, 241, 257), (3, 5, 7, 13, 19, 37, 73), (3, 5, 7, 13, 19, 37, 109), (3, 5, 7, 13, 11, 31, 41, 151, 331)]
assert all(covers5(S) and all(not covers5([q for q in S if q != p]) for p in S) for S in sets5)
assert all(any((k5 * 2**n + 1) % p == 0 for p in sets5[0]) for n in range(1, 3000))
print("k5 =", k5, "has 5 distinct minimal covering sets")
print("all claims verified")

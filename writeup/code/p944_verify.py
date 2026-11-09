# Puzzle 944: P palprime, intertwining the digit D gives A + D B. Checks every k = 5 solution in data/p944_results.txt
# (P palindrome and prime, A + D B prime for D = 0..5 with PARI isprime proofs, A + 6 B composite), Rivera's P,
# Cleary's two 8-of-10 examples, and the identity N_D = A + D B against direct digit intertwining.
# Run from the repo root: .venv/bin/python3 writeup/code/p944_verify.py
import re
from sympy import isprime
import cypari2
pari = cypari2.Pari()
def intertwine(P, D): return int(str(D).join(str(P)))
def spread(P): return int("0".join(str(P)))
def B(n): return 10 * (100**(n - 1) - 1) // 99
def run(P):
    n = len(str(P)); A = spread(P); k = -1
    while k < 9 and pari.isprime(A + (k + 1) * B(n)) == 1: k += 1
    return k
for P in (3273534353723,) + tuple(int(m.group(1)) for m in re.finditer(r"k=5 P=(\d+)", open("writeup/data/p944_results.txt").read())):
    s = str(P); assert s == s[::-1] and pari.isprime(P) == 1
    for D in range(10): assert intertwine(P, D) == spread(P) + D * B(len(s))
    assert run(P) == 5, P
    print(f"P={P}: palprime, A + D B prime for D = 0..5 (proven), composite for D = 6")
for P, pattern in ((7808765306035678087, "1111011101"), (7792235971795322977, "1111011101")):   # Cleary: D = 0..3, 5, 6, 7, 9
    got = "".join("1" if pari.isprime(intertwine(P, D)) == 1 else "0" for D in range(10)); assert got == pattern, got
print("Cleary's 19-digit examples: primes exactly for D = 0, 1, 2, 3, 5, 6, 7, 9")
for n in (13, 19, 25): print(f"n={n}: B divisible by", [q for q in (2, 3, 5, 7, 11, 13, 37, 101, 9901) if B(n) % q == 0])
print("all checks passed")

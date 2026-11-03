"""Puzzle 1145: repeat the bounded covering-set search.

    gcc -O2 -o writeup/code/p1145_dfs writeup/code/p1145_dfs.c -lm
    python3 writeup/code/p1145_search.py pool          # check the prime pool (seconds)
    python3 writeup/code/p1145_search.py run 360 1e29 4  # one period: N, bound on the product, parts
    python3 writeup/code/p1145_search.py all           # every run in the write-up (about 10 minutes)

Pool: for odd n = 2m + 1, a prime p can cover a residue class of m only if (-6/p) = 1,
and the class has period ord_p(4). ord_p(4) | N exactly when p divides 4^N - 1, which is
the product of Phi_e(2) over e | 2N. writeup/data/p1145_cyclo.txt holds the complete
factorisation of every Phi_e(2) needed (PARI factor(); e = 1440 from factordb.com).
The pool step checks every line: the factors multiply to Phi_e(2) and each one is a
proven prime. So the pool for period N is every admissible prime below 2^64 with
ord_p(4) | N, with nothing left out.

Input for p1145_dfs.c (stdin): "N ln(bound) count", then per prime "p ord_p(4)" and a
line with ord_p(4) numbers s_r, s_r^2 = -(6 * 4^r)^(-1) mod p. Arguments: time limit (s),
report threshold, exact-capacity depth, MRV flag, symmetry flag (fixes the classes of 5
and 7; the 6 shifts are restored in the CRT step), parts, part, split depth, knapsack flag.
Output per prime set: number of covers and the least a over all covers and all sign
choices (a = +-1 mod 6 and a square root per prime). The parts split the tree at depth 2;
their union is the whole search.
"""
import math
import subprocess
import sys
from pathlib import Path

import cypari2
from sympy import cyclotomic_poly, divisors, isprime, legendre_symbol, n_order, sqrt_mod

HERE = Path(__file__).resolve().parent
CYCLO = HERE.parent / "data" / "p1145_cyclo.txt"
DFS = HERE / "p1145_dfs"
PMAX = 2**64
BASE = (5, 7, 11, 31, 97, 151, 241, 673, 1321)
EXPECTED_360 = [
    (153479820268467961, 34560, (73, 631, 23311)),
    (2102924150668420411, 34560, (73, 631, 54001)),
    (11332178567623762313, 103680, (73, 433, 631, 23311)),
    (18146510236193475391, 103680, (73, 577, 631, 23311)),
    (34465873137135828277, 69120, (73, 433, 631, 38737)),
    (36893778994594030543, 103680, (73, 433, 631, 54001)),
    (155473688074716734167, 34560, (73, 23311, 54001)),
    (217094373626187598283, 34560, (631, 23311, 54001)),
    (222263303128057762717, 69120, (73, 433, 631, 61681)),
    (345709756933378561391, 34560, (73, 631, 18837001)),
    (706834171741275231299, 103680, (73, 577, 631, 54001)),
]
RUNS = [(24, 1e30, 1), (36, 1e30, 1), (48, 1e30, 1), (72, 1e30, 1), (96, 1e30, 1),
        (120, 1e30, 1), (144, 1e30, 1), (180, 1e30, 1), (240, 1e30, 1), (288, 1e30, 1),
        (480, 1e27, 1), (360, 1e29, 4), (720, 1e26, 16)]

pari = cypari2.Pari()


def read_cyclo():
    table = {}
    for line in CYCLO.read_text().splitlines():
        e, rest = line.split(":")
        table[int(e)] = [tuple(map(int, f.split("^"))) if "^" in f else (int(f), 1) for f in rest.split()]
    return table


def check_cyclo(table):
    for e, factors in table.items():
        assert math.prod(p**k for p, k in factors) == int(cyclotomic_poly(e, 2)), e
        for p, k in factors:
            assert isprime(p) if p < PMAX else bool(pari.isprime(p)), (e, p)
    return table


def pool(N, table):
    primes = set()
    for e in divisors(2 * N):
        assert e in table, f"Phi_{e}(2) missing from {CYCLO.name}"
        primes |= {p for p, k in table[e] if 3 < p < PMAX and legendre_symbol(-6 % p, p) == 1}
    return sorted(primes)


def dfs_input(N, bound, primes):
    lines = [f"{N} {math.log(bound)} {len(primes)}"]
    for p in primes:
        o = n_order(4, p)
        assert N % o == 0
        lines.append(f"{p} {o}")
        lines.append(" ".join(str(sqrt_mod(-pow(6 * pow(4, r, p), -1, p) % p, p)) for r in range(o)))
    return "\n".join(lines) + "\n"


def run(N, bound, parts, table, seconds=3600):
    primes = pool(N, table)
    text = dfs_input(N, bound, primes)
    sets = {}
    for part in range(parts):
        args = ["nice", "-n", "19", str(DFS), str(seconds), "0", "5", "1", "1", str(parts), str(part), "2", "1"]
        out = subprocess.run(args, input=text, capture_output=True, text=True, check=True).stdout.splitlines()
        done = [l for l in out if l.startswith("DONE")]
        assert done and done[0].endswith("complete"), (N, part, done)
        for l in out:
            if not l.startswith("SET"):
                continue
            w = l.split()
            key = tuple(int(x) for x in w[w.index("primes") + 1:])
            covers, a = int(w[1].split("=")[1]), int(w[2].split("=")[1])
            old = sets.get(key, (0, a))
            sets[key] = (old[0] + covers, min(old[1], a))
        print(f"  N = {N} part {part + 1}/{parts}: {done[0].split()[-2]}", flush=True)
    return primes, sorted(sets.items(), key=lambda kv: kv[1][1])


def report(N, bound, primes, sets):
    print(f"N = {N}: pool of {len(primes)} primes, bound {bound:.0e}: {len(sets)} covering sets")
    for key, (covers, a) in sets:
        print(f"  a = {a}  covers = {covers}  product = {math.prod(key):.3e}  primes {' '.join(map(str, key))}")


def main():
    table = check_cyclo(read_cyclo())
    print(f"{CYCLO.name}: {len(table)} factorisations checked (product and primality)")
    mode = sys.argv[1] if len(sys.argv) > 1 else "pool"
    if mode == "pool":
        for N, bound, parts in RUNS:
            print(f"N = {N}: {' '.join(f'{p}[{n_order(4, p)}]' for p in pool(N, table))}")
        return
    if not DFS.exists():
        sys.exit(f"compile first: gcc -O2 -o {DFS} {DFS}.c -lm")
    if mode == "run":
        N, bound, parts = int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
        report(N, bound, *run(N, bound, parts, table))
        return
    results = {}
    for N, bound, parts in RUNS:
        primes, sets = run(N, bound, parts, table)
        report(N, bound, primes, sets)
        results[N] = sets
    for N in (24, 36, 48, 72, 96, 120, 144, 180, 240, 288, 480):
        assert results[N] == [], N
    expected = [(tuple(sorted(BASE + extra)), (covers, a)) for a, covers, extra in EXPECTED_360]
    assert results[360] == expected
    assert [(key, a) for key, (covers, a) in results[720]] == [(key, a) for key, (covers, a) in expected[:2]]
    print("all claims of the write-up reproduced")


if __name__ == "__main__":
    main()

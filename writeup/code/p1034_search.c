/* Puzzles 1034 / 1091.  A chain p, p + 6k, p + 6k + 6(k+1), ... of primes is a run of consecutive n with
   f(n) = 3n^2 + 3n + c prime (p = f(n0), 6k = 6(n0+1)).  f has no root mod q (q >= 5) iff D = 9 - 12c is a
   non-residue mod q; the Hardy-Littlewood constant C(c) = 3 prod_{q>=5} (1 - w_q/q)/(1 - 1/q) measures the
   prime density of f (w_q = 1 + (D|q)); c must be odd and prime to 3 (else f is always even or divisible by 3).
   For every c in [-CMAX, CMAX] with (D|q) = -1 for q = 5..23 and C(c) (over q <= 2000) >= CMIN, scan n = 0..N
   and report every run of >= RMIN consecutive prime values (deterministic Miller-Rabin, 64-bit).
   Usage: p1034_search CMAX N CMIN RMIN        gcc -O2 -fopenmp -o p1034_search p1034_search.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
typedef unsigned __int128 u128;
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) { return (u128)a * b % m; }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) { uint64_t r = 1; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static int isprime64(uint64_t n) {
    if (n < 2) return 0;
    static const uint64_t sp[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) { if (n == sp[i]) return 1; if (n % sp[i] == 0) return 0; }
    uint64_t d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        uint64_t x = powmod(sp[i], d, n); if (x == 1 || x == n - 1) continue;
        int ok = 0; for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}
static int legendre(int64_t a, uint32_t q) {      /* (a|q) for odd prime q */
    a %= (int64_t)q; if (a < 0) a += q; if (a == 0) return 0;
    return powmod(a, (q - 1) / 2, q) == 1 ? 1 : -1;
}
int main(int argc, char **argv) {
    int64_t CMAX = atoll(argv[1]); int N = atoi(argv[2]); double CMIN = atof(argv[3]); int RMIN = atoi(argv[4]);
    static uint32_t qs[400]; int nq = 0;
    for (uint32_t q = 5; q <= 2000; q++) { int pr = 1; for (uint32_t d = 2; d * d <= q; d++) if (q % d == 0) { pr = 0; break; } if (pr) qs[nq++] = q; }
    long long ncand = 0, nscan = 0; double bestC = 0; int bestrun = 0;
    #pragma omp parallel for schedule(dynamic, 1000) reduction(+:ncand,nscan)
    for (int64_t c = -CMAX; c <= CMAX; c++) {
        if (!(c & 1) || c % 3 == 0) continue;
        int64_t D = 9 - 12 * c;
        int ok = 1;
        for (int i = 0; i < 7; i++) if (legendre(D, qs[i]) != -1) { ok = 0; break; }   /* q = 5 .. 23 */
        if (!ok) continue;
        ncand++;
        double C = 3.0;
        for (int i = 0; i < nq; i++) { int w = 1 + legendre(D, qs[i]); C *= (1.0 - (double)w / qs[i]) / (1.0 - 1.0 / qs[i]); }
        if (C < CMIN) continue;
        nscan++;
        int run = 0;
        for (int n = 0; n <= N; n++) {
            int64_t v = 3LL * n * n + 3LL * n + c;
            int pr = v > 1 && isprime64((uint64_t)v);
            if (pr) run++;
            if (!pr || n == N) {
                if (run >= RMIN) {
                    int n0 = (pr ? n : n - 1) - run + 1; int64_t p = 3LL * n0 * n0 + 3LL * n0 + c;
                    #pragma omp critical
                    { printf("run %d: c=%lld C=%.2f n0=%d p=%lld step=%d\n", run, (long long)c, C, n0, (long long)p, 6 * (n0 + 1)); fflush(stdout);
                      if (run > bestrun) bestrun = run; }
                }
                run = 0;
            }
        }
        #pragma omp critical
        if (C > bestC) bestC = C;
    }
    printf("c in [-%lld, %lld]: %lld pass the residue filter, %lld scanned (C >= %.1f), largest C %.2f, longest run %d\n",
           (long long)CMAX, (long long)CMAX, ncand, nscan, CMIN, bestC, bestrun);
    return 0;
}

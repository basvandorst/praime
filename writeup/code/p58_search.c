/* Puzzle 58 (and 54): k consecutive primes p with 4p^2 + 1 prime.  Earliest k = 8 is unknown; it lies in
   (2e11, 1719004000057] (k = 7 minimal below 2e11, a k = 9 run starts at 1719004000057).
   Segmented sieve over [LO, HI); for each prime p test 4p^2 + 1 against the primes q = 1 (mod 4) up to 200
   (q | 4p^2 + 1 iff p^2 = -1/4 mod q, two residues), count consecutive passing primes, and print every run of
   >= K consecutive candidates (the run and the preceding prime) for exact verification in Python.
   Usage: p58_search LO HI K [threads]          gcc -O2 -fopenmp -o p58_search p58_search.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#define NQ 24
static uint32_t Q[NQ]; static uint32_t R1[NQ], R2[NQ]; static int nq = 0;
static uint64_t powmod32(uint64_t a, uint64_t e, uint64_t m) { uint64_t r = 1; a %= m; while (e) { if (e & 1) r = r * a % m; a = a * a % m; e >>= 1; } return r; }
static int passes(uint64_t p) {
    for (int i = 0; i < nq; i++) { uint32_t r = p % Q[i]; if (r == R1[i] || r == R2[i]) return 0; }
    return 1;
}
int main(int argc, char **argv) {
    uint64_t LO = strtoull(argv[1], 0, 10), HI = strtoull(argv[2], 0, 10); int K = atoi(argv[3]); int th = argc > 4 ? atoi(argv[4]) : 2;
    /* residues: p^2 = -1/4 (mod q) for primes q = 1 mod 4 */
    for (uint32_t q = 5; q <= 200 && nq < NQ; q += 4) {
        int pr = 1; for (uint32_t d = 2; d * d <= q; d++) if (q % d == 0) { pr = 0; break; }
        if (!pr) continue;
        uint64_t inv4 = powmod32(4, q - 2, q), t = (q - inv4) % q;       /* -1/4 */
        uint64_t r = 0; for (uint64_t x = 1; x < q; x++) if (x * x % q == t) { r = x; break; }
        if (!r) continue;
        Q[nq] = q; R1[nq] = r; R2[nq] = q - r; nq++;
    }
    fprintf(stderr, "%d filter primes, first %u last %u\n", nq, Q[0], Q[nq - 1]);
    uint64_t Rt = 1; while (Rt * Rt < HI) Rt++;
    char *small = calloc(Rt + 1, 1); uint32_t *bp = malloc(4 * (Rt / 4 + 100)); int nb = 0;
    for (uint64_t i = 3; i <= Rt; i += 2) if (!small[i]) { bp[nb++] = i; for (uint64_t j = i * i; j <= Rt; j += 2 * i) small[j] = 1; }
    const uint64_t SEG = 1ull << 24;                       /* odd numbers per segment */
    omp_set_num_threads(th);
    #pragma omp parallel
    {
        int t = omp_get_thread_num(), T = omp_get_num_threads();
        uint64_t span = (HI - LO) / T, blo = LO + span * t, bhi = (t == T - 1) ? HI : blo + span;
        uint8_t *seg = malloc(SEG);
        int run = 0; uint64_t runp[64]; uint64_t prevp = 0;
        /* start a little before blo so that runs crossing the boundary are seen by this thread when they END inside */
        uint64_t start = blo > LO + 2000 ? blo - 2000 : LO;
        for (uint64_t lo = start | 1; lo < bhi; lo += 2 * SEG) {
            uint64_t hi = lo + 2 * SEG; if (hi > bhi) hi = bhi;
            memset(seg, 0, SEG);
            for (int k = 0; k < nb; k++) {
                uint64_t q = bp[k]; if (q * q >= hi) break;
                uint64_t st = q * q; if (st < lo) { st = ((lo + q - 1) / q) * q; if (!(st & 1)) st += q; }
                for (uint64_t m = st; m < hi; m += 2 * q) seg[(m - lo) / 2] = 1;
            }
            for (uint64_t i = 0; i < SEG; i++) if (!seg[i]) {
                uint64_t p = lo + 2 * i; if (p >= hi) break;
                if (passes(p)) { if (run < 64) runp[run] = p; run++; }
                else {
                    if (run >= K && runp[0] >= blo - 0 && runp[0] < bhi) {   /* report runs starting in this block */
                        #pragma omp critical
                        { printf("run %d start %llu:", run, (unsigned long long)runp[0]);
                          for (int j = 0; j < run && j < 64; j++) printf(" %llu", (unsigned long long)runp[j]);
                          printf(" | next %llu\n", (unsigned long long)p); fflush(stdout); }
                    }
                    run = 0;
                }
                prevp = p;
            }
        }
        free(seg);
    }
    return 0;
}

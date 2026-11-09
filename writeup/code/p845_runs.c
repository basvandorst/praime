/* Puzzle 845 / OEIS A275508, A066793.  Omega(n) = number of prime factors with multiplicity.
   Segmented sieve of the parity of Omega(n) (Liouville) for n < N: for every prime power p^k <= N
   toggle the parity of its multiples and multiply a cofactor tracker by p; a leftover cofactor > 1
   is one more prime (it has no factor <= sqrt(N)), so it toggles once more.
   Then the runs of consecutive n with even Omega: the first run of exactly L, for every L, and the
   record runs.  Threads take contiguous blocks and recompute a margin so that runs crossing a block
   boundary are attributed to the block where they start.
   Usage: p845_runs N [threads]          gcc -O2 -fopenmp -o p845_runs p845_runs.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#define MAXL 128
#define MARGIN 256
int main(int argc, char **argv) {
    uint64_t N = strtoull(argv[1], 0, 10); int th = argc > 2 ? atoi(argv[2]) : 2;
    uint64_t R = 1; while (R * R < N) R++;
    char *comp = calloc(R + 2, 1); uint32_t *pr = malloc(4 * (R / 2 + 100)); int np = 0;
    for (uint64_t i = 2; i <= R; i++) if (!comp[i]) { pr[np++] = i; for (uint64_t j = i * i; j <= R; j += i) comp[j] = 1; }
    uint64_t first[MAXL]; memset(first, 0, sizeof first);   /* first start of a run of exactly L */
    uint64_t recstart[MAXL]; memset(recstart, 0, sizeof recstart);
    const uint64_t SEG = 1u << 22;
    omp_set_num_threads(th);
    #pragma omp parallel
    {
        int t = omp_get_thread_num(), T = omp_get_num_threads();
        uint64_t blo = (N / T) * t, bhi = (t == T - 1) ? N : (N / T) * (t + 1);
        if (blo < 1) blo = 1;
        uint64_t lo0 = blo > 1 ? blo - 1 : 1;              /* one before, to know whether a run started earlier */
        uint64_t hi0 = bhi + MARGIN; if (hi0 > N + MARGIN) hi0 = N + MARGIN;
        uint8_t *par = malloc(SEG + MARGIN); uint64_t *rem = malloc(8 * (SEG + MARGIN));
        uint64_t myfirst[MAXL]; memset(myfirst, 0, sizeof myfirst);
        int run = 0; uint64_t runstart = 0; int prev_even_before = 0;
        for (uint64_t lo = lo0; lo < hi0; lo += SEG) {
            uint64_t hi = lo + SEG; if (hi > hi0) hi = hi0; uint64_t len = hi - lo;
            memset(par, 0, len); for (uint64_t i = 0; i < len; i++) rem[i] = 1;
            for (int k = 0; k < np; k++) {
                uint64_t p = pr[k]; if (p * p > hi) { /* primes above sqrt(hi) only as single factors: handled by leftover */ break; }
                for (uint64_t q = p; q < hi; q *= p) {    /* each prime power p, p^2, ... */
                    uint64_t st = ((lo + q - 1) / q) * q;
                    for (uint64_t m = st; m < hi; m += q) { par[m - lo] ^= 1; rem[m - lo] *= p; }
                    if (q > hi / p) break;
                }
            }
            for (uint64_t i = 0; i < len; i++) {
                uint64_t n = lo + i; int even = !(par[i] ^ (rem[i] != n));   /* leftover prime toggles */
                if (n == lo0 && lo0 == blo - 1) { prev_even_before = even; continue; }   /* sentinel */
                if (even) { if (run == 0) runstart = n; run++; }
                else {
                    if (run > 0) {
                        int L = run < MAXL ? run : MAXL - 1;
                        int startedinside = (runstart >= blo) && !(runstart == blo && prev_even_before);
                        if (startedinside && runstart < bhi && !myfirst[L]) myfirst[L] = runstart;
                        run = 0;
                    }
                }
                if (n >= bhi + 1 && run == 0) break;        /* past the block and no open run */
            }
            if (lo + SEG >= hi0) break;
        }
        #pragma omp critical
        for (int L = 1; L < MAXL; L++) if (myfirst[L] && (!first[L] || myfirst[L] < first[L])) first[L] = myfirst[L];
        free(par); free(rem);
    }
    printf("N = %llu\n", (unsigned long long)N);
    uint64_t best = 0;
    for (int L = 1; L < MAXL; L++) if (first[L]) printf("exactly %d: %llu\n", L, (unsigned long long)first[L]);
    printf("record runs (start, length):");
    for (int L = 1; L < MAXL; L++) if (first[L] && (!best || first[L] < best)) { }
    /* records: scan all (start, L) pairs in order of start */
    for (;;) { uint64_t mn = 0; int mL = 0; for (int L = 1; L < MAXL; L++) if (first[L] && first[L] > best && (!mn || first[L] < mn)) { mn = first[L]; mL = L; }
        if (!mn) break; best = mn; static int rec = 0; if (mL > rec) { rec = mL; printf(" (%llu, %d)", (unsigned long long)mn, mL); } }
    printf("\n"); return 0;
}

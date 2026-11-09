/* Puzzle 1185 / OEIS A247148: primes p that divide LS(p), their look-and-say value.
   Segmented sieve of Eratosthenes (odd numbers only) up to N; for every prime p the
   look-and-say string is built as a 128-bit integer (runs of length >= 10 get a
   two-digit count, handled generically) and reduced mod p.
   Usage: p1185_search N [threads]      Output: "p LS(p)" per solution, then pi(N).
   gcc -O2 -fopenmp -o p1185_search p1185_search.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef unsigned __int128 u128;
static int ls_divides(uint64_t p, u128 *lsout) {
    int d[24], n = 0; uint64_t q = p;
    while (q) { d[n++] = q % 10; q /= 10; }
    u128 ls = 0;
    for (int i = n - 1; i >= 0;) {            /* most significant digit first */
        int j = i; while (j > 0 && d[j - 1] == d[i]) j--;
        int c = i - j + 1;                    /* run length */
        if (c >= 10) ls = ls * 100 + c; else ls = ls * 10 + c;
        ls = ls * 10 + d[i];
        i = j - 1;
    }
    *lsout = ls;
    return ls % p == 0;
}
static void print128(u128 x, char *buf) {
    char t[48]; int k = 0;
    if (!x) t[k++] = '0';
    while (x) { t[k++] = '0' + (int)(x % 10); x /= 10; }
    for (int i = 0; i < k; i++) buf[i] = t[k - 1 - i];
    buf[k] = 0;
}
int main(int argc, char **argv) {
    uint64_t N = strtoull(argv[1], 0, 10);
    int th = argc > 2 ? atoi(argv[2]) : 2;
    omp_set_num_threads(th);
    /* base primes to sqrt(N) */
    uint64_t R = 1; while (R * R < N) R++;
    char *small = calloc(R + 1, 1);
    uint32_t *bp = malloc(sizeof(uint32_t) * (R / 2 + 10)); int nb = 0;
    for (uint64_t i = 3; i <= R; i += 2) if (!small[i]) { bp[nb++] = i; for (uint64_t j = i * i; j <= R; j += 2 * i) small[j] = 1; }
    const uint64_t SEG = 1ull << 24;          /* odd numbers per segment: 2^24 -> span 2^25 */
    uint64_t nseg = (N / 2 + SEG - 1) / SEG;
    uint64_t total = 0; int found = 0;
    if (N >= 2) { printf("2 12\n"); total = 1; found = 1; }   /* the even prime, LS(2) = 12 */
    #pragma omp parallel
    {
        uint8_t *seg = malloc(SEG);
        #pragma omp for schedule(dynamic) reduction(+:total,found)
        for (uint64_t s = 0; s < nseg; s++) {
            uint64_t lo = 2 * s * SEG + 1;          /* first odd number of the segment */
            uint64_t hi = lo + 2 * SEG;             /* exclusive */
            memset(seg, 0, SEG);
            if (s == 0) seg[0] = 1;                 /* 1 is not prime */
            for (int k = 0; k < nb; k++) {
                uint64_t q = bp[k]; if (q * q >= hi) break;
                uint64_t st = q * q; if (st < lo) { st = ((lo + q - 1) / q) * q; if (!(st & 1)) st += q; }
                for (uint64_t m = st; m < hi; m += 2 * q) seg[(m - lo) / 2] = 1;
            }
            for (uint64_t i = 0; i < SEG; i++) if (!seg[i]) {
                uint64_t p = lo + 2 * i; if (p >= N) break;
                total++;
                u128 ls;
                if (ls_divides(p, &ls)) {
                    char buf[48]; print128(ls, buf);
                    #pragma omp critical
                    { printf("%llu %s\n", (unsigned long long)p, buf); fflush(stdout); }
                    found++;
                }
            }
        }
        free(seg);
    }
    printf("pi(%llu) = %llu, solutions = %d\n", (unsigned long long)N, (unsigned long long)total, found);
    return 0;
}

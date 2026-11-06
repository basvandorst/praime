/* Puzzle 990 Q1: all 4-term arithmetic progressions n, n+d, n+2d, n+3d (d > 0) with
   sigma equal on all four terms and n + 3d <= N.  sigma by a divisor sieve, then the
   indices sorted by sigma and every group of equal sigma scanned for 4-term APs.
   Usage: p990_ap4 N        (N <= 2e8; sigma fits in 32 bits there; memory 4N bytes plus the largest class mod 16) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static uint32_t *sig;
static int cmp(const void *a, const void *b) {
    uint32_t x = sig[*(const uint32_t *)a], y = sig[*(const uint32_t *)b];
    return x < y ? -1 : x > y ? 1 : (*(const uint32_t *)a < *(const uint32_t *)b ? -1 : 1);
}
int main(int argc, char **argv) {
    uint32_t N = (uint32_t)atoll(argv[1]);
    sig = calloc(N + 1, 4);
    for (uint32_t d = 1; d <= N; d++) for (uint32_t m = d; m <= N; m += d) sig[m] += d;
    /* 16 passes: the numbers with sigma = r (mod 16) are sorted and scanned together; an AP's
       four terms share their sigma, so they land in the same pass.  Memory: sigma + N/16 indices. */
    long found = 0;
    uint32_t per[16] = {0}, mx = 0;
    for (uint32_t i = 1; i <= N; i++) per[sig[i] & 15]++;
    for (int r = 0; r < 16; r++) if (per[r] > mx) mx = per[r];
    uint32_t *idx = malloc((size_t)mx * 4);   /* sigma is biased mod 16, so the classes differ in size */
    for (int r = 0; r < 16; r++) {
        uint32_t cnt = 0;
        for (uint32_t i = 1; i <= N; i++) if ((sig[i] & 15) == (uint32_t)r) idx[cnt++] = i;
        qsort(idx, cnt, 4, cmp);
        for (uint32_t a = 0; a < cnt;) {
            uint32_t b = a; while (b < cnt && sig[idx[b]] == sig[idx[a]]) b++;
            if (b - a >= 4) {
                for (uint32_t i = a; i < b; i++) for (uint32_t j = i + 1; j < b; j++) {
                    uint64_t d = idx[j] - idx[i];
                    uint64_t t3 = idx[i] + 3 * d; if (t3 > N) break;
                    uint32_t t2 = idx[i] + 2 * d, ok = 0;
                    uint32_t lo = j + 1, hi = b;
                    while (lo < hi) { uint32_t mid = (lo + hi) / 2; if (idx[mid] < t2) lo = mid + 1; else hi = mid; }
                    if (lo < b && idx[lo] == t2) {
                        uint32_t lo2 = lo + 1, hi2 = b;
                        while (lo2 < hi2) { uint32_t mid = (lo2 + hi2) / 2; if (idx[mid] < t3) lo2 = mid + 1; else hi2 = mid; }
                        if (lo2 < b && idx[lo2] == t3) ok = 1;
                    }
                    if (ok) { found++; printf("AP4: n = %u d = %llu sigma = %u : %u %u %u %llu\n", idx[i], (unsigned long long)d, sig[idx[i]], idx[i], idx[j], t2, (unsigned long long)t3); fflush(stdout); }
                }
            }
            a = b;
        }
    }
    printf("N = %u: %ld four-term APs with equal sigma and last term <= N\n", N, found);
    return 0;
}

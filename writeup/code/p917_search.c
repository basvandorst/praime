/* Puzzle 917.  P prime with 2d digits, A = first d digits, B = last d digits (B may have
   leading zeros), Q = P + A*B, halved if even.  Find prime Q reached by >= K primes P.
   Method: for every prime P with 2d digits (segmented sieve) compute Q and increment a byte
   counter for Q; the Q range [lo, hi) is processed in chunks that fit in memory.  Q with
   count >= K are tested for primality (deterministic Miller-Rabin, 64-bit) and their
   preimages are recovered by solving B = (T - A 10^d)/(A+1), T in {Q, 2Q}, over A.
   Also prints the histogram of counts >= 3 among prime Q (for the Poisson model).
   Usage: p917_search d K [chunk_log2=30]           gcc -O2 -o p917_search p917_search.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
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
static uint64_t POW10[20];
int main(int argc, char **argv) {
    int d = atoi(argv[1]), K = atoi(argv[2]); int CL = argc > 3 ? atoi(argv[3]) : 30;
    POW10[0] = 1; for (int i = 1; i < 20; i++) POW10[i] = POW10[i - 1] * 10;
    uint64_t M = POW10[d], Pmin = POW10[2 * d - 1], Pmax = POW10[2 * d];
    /* base primes up to sqrt(Pmax) */
    uint64_t R = 1; while (R * R < Pmax) R++;
    char *small = calloc(R + 1, 1); uint32_t *bp = malloc(4 * (R / 2 + 10)); int nb = 0;
    for (uint64_t i = 3; i <= R; i += 2) if (!small[i]) { bp[nb++] = i; for (uint64_t j = i * i; j <= R; j += 2 * i) small[j] = 1; }
    uint64_t Qlo = Pmin / 2, Qhi = Pmax + (M - 1) * (M - 1) + 1;   /* Q range */
    uint64_t CH = 1ull << CL;                                        /* odd Q per chunk */
    uint8_t *cnt = malloc(CH);
    const uint64_t SEG = 1ull << 22; uint8_t *seg = malloc(SEG);
    uint64_t hist[9] = {0}; uint64_t nP = 0;
    for (uint64_t clo = Qlo | 1; clo < Qhi; clo += 2 * CH) {
        uint64_t chi = clo + 2 * CH;                                 /* odd Q in [clo, chi) */
        memset(cnt, 0, CH); nP = 0;
        for (uint64_t lo = Pmin | 1; lo < Pmax; lo += 2 * SEG) {
            uint64_t hi = lo + 2 * SEG; if (hi > Pmax) hi = Pmax;
            memset(seg, 0, SEG);
            for (int k = 0; k < nb; k++) {
                uint64_t q = bp[k]; if (q * q >= hi) break;
                uint64_t st = q * q; if (st < lo) { st = ((lo + q - 1) / q) * q; if (!(st & 1)) st += q; }
                for (uint64_t m = st; m < hi; m += 2 * q) seg[(m - lo) / 2] = 1;
            }
            for (uint64_t i = 0; i < SEG; i++) if (!seg[i]) {
                uint64_t P = lo + 2 * i; if (P >= hi) break;
                nP++;
                uint64_t A = P / M, B = P % M, Q = P + A * B; if (!(Q & 1)) Q >>= 1;
                if (Q >= clo && Q < chi && (Q & 1)) { uint64_t j = (Q - clo) / 2; if (cnt[j] < 255) cnt[j]++; }
            }
        }
        for (uint64_t j = 0; j < CH; j++) if (cnt[j] >= 3) {
            uint64_t Q = clo + 2 * j; if (!isprime64(Q)) continue;
            int c = cnt[j] > 8 ? 8 : cnt[j]; hist[c]++;
            if (cnt[j] >= K) {
                printf("k=%d Q=%llu P:", cnt[j], (unsigned long long)Q);
                for (int t = 1; t <= 2; t++) {                     /* T = Q (odd branch) and 2Q (even branch) */
                    uint64_t T = t * Q;
                    uint64_t Amin = T / M > 0 ? (T / M) / 2 : 0; if (Amin < M / 10) Amin = M / 10;
                    for (uint64_t A = Amin; A <= T / M && A < M; A++) {
                        if ((T - A * M) % (A + 1)) continue;
                        uint64_t B = (T - A * M) / (A + 1); if (B >= M) continue;
                        uint64_t P = A * M + B; if (P < Pmin || P >= Pmax) continue;
                        uint64_t q = P + A * B; int even = !(q & 1); if (even) q >>= 1;
                        if (q != Q || !isprime64(P)) continue;
                        printf(" %llu%s", (unsigned long long)P, B < M / 10 ? "(B<10^(d-1))" : "");
                    }
                }
                printf("\n"); fflush(stdout);
            }
        }
        fprintf(stderr, "chunk [%llu, %llu) done, %llu primes P\n", (unsigned long long)clo, (unsigned long long)chi, (unsigned long long)nP);
    }
    printf("primes P with %d digits: %llu\n", 2 * d, (unsigned long long)nP);
    printf("histogram of prime Q by count (3..8+):"); for (int c = 3; c <= 8; c++) printf(" %llu", (unsigned long long)hist[c]); printf("\n");
    return 0;
}

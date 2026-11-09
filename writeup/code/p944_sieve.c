/* Puzzle 944.  Palindromes P with n digits (n odd); A = P spread with zeros (digit i at 10^(2i)),
   B = 10 (100^(n-1) - 1)/99; the intertwined numbers are A + D B.  For D = 0..K-1 to be prime,
   A must avoid the residues -D B mod q for every small prime q not dividing B.  Depth-first over the
   (n+1)/2 free digits with A mod q carried along; survivors are tested for P prime (64-bit
   Miller-Rabin, n <= 19) and printed; the AP itself is tested in Python with gmpy2.
   Usage: p944_sieve n K        gcc -O2 -o p944_sieve p944_sieve.c */
#include <stdio.h>
#include <stdlib.h>
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
#define NQ 40
static int n, h, K, nq; static uint32_t Q[NQ]; static uint64_t forb[NQ][4];     /* bitmask of forbidden residues (q < 256) */
static uint32_t w[NQ][32];      /* w[j][i] = weight of digit position i (from the top, i = 0..n-1) in A mod Q[j] */
static int dig[32]; static uint32_t amod[32][NQ]; static long long leaves = 0, surv = 0, primes = 0;
static uint64_t Pval(void) { uint64_t v = 0; for (int i = 0; i < n; i++) v = v * 10 + dig[i]; return v; }
static void dfs(int i) {           /* choose digit i (0-based from the top) and its mirror n-1-i */
    if (i == h) {
        leaves++;
        for (int j = 0; j < nq; j++) { uint32_t r = amod[h][j]; if (forb[j][r >> 6] >> (r & 63) & 1) return; }
        surv++;
        uint64_t P = Pval(); if (!isprime64(P)) return;
        primes++; printf("%llu\n", (unsigned long long)P);
        return;
    }
    int mirror = n - 1 - i;
    for (int d = (i == 0 ? 1 : 0); d < 10; d++) {
        if (i == 0 && (d == 2 || d == 4 || d == 5 || d == 6 || d == 8)) continue;
        dig[i] = d; dig[mirror] = d;
        for (int j = 0; j < nq; j++) {
            uint32_t add = (uint64_t)d * (w[j][i] + (mirror != i ? w[j][mirror] : 0)) % Q[j];
            amod[i + 1][j] = (amod[i][j] + add) % Q[j];
        }
        dfs(i + 1);
    }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); K = atoi(argv[2]); h = (n + 1) / 2;
    /* B mod q and the weights 100^(n-1-i) mod q */
    static const uint32_t cand[] = {11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199};
    nq = 0;
    for (unsigned c = 0; c < sizeof cand / sizeof cand[0] && nq < NQ; c++) {
        uint32_t q = cand[c];
        uint64_t Bq = 0, pw = 1;                          /* B = 10 * sum_{j=0}^{n-2} 100^j */
        for (int j = 0; j < n - 1; j++) { Bq = (Bq + pw) % q; pw = pw * 100 % q; }
        Bq = Bq * 10 % q;
        if (Bq == 0) continue;                            /* q | B: every term = A mod q, no sieve */
        Q[nq] = q; forb[nq][0] = forb[nq][1] = forb[nq][2] = forb[nq][3] = 0;
        for (int D = 0; D < K; D++) { uint32_t r = (q - (uint64_t)D * Bq % q) % q; forb[nq][r >> 6] |= 1ull << (r & 63); }
        uint64_t pwr = 1;
        for (int i = n - 1; i >= 0; i--) { w[nq][i] = pwr; pwr = pwr * 100 % q; }
        nq++;
    }
    fprintf(stderr, "n=%d K=%d sieve primes (not dividing B): %d, first %u last %u\n", n, K, nq, Q[0], Q[nq - 1]);
    for (int j = 0; j < nq; j++) amod[0][j] = 0;
    dfs(0);
    fprintf(stderr, "palindromes %lld, sieve survivors %lld, palprimes among them %lld\n", leaves, surv, primes);
    return 0;
}

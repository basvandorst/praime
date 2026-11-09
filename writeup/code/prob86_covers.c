// Problem 86: Sierpinski numbers divisible by 3.
//
// k odd, 3 | k, and k*2^n + 1 composite for all n >= 1 because every n is
// covered by a prime p of a covering set S (3 cannot be in S).  Prime p with
// m = ord_p(2) covers n == a (mod m) exactly when k == -2^(-a) (mod p).
//
// This program enumerates every "canonical" cover of Z/L by residue classes
// a_p (mod m_p) drawn from a given list of primes (one class per prime,
// each prime optional).  A cover is canonical when every used prime is, at
// the moment it is chosen, the first prime in list order that covers the
// smallest still uncovered point; every cover contains a canonical sub-cover,
// and the sub-cover's k is not larger (fewer congruences).  So the minimum
// over canonical covers of the least k == 3 (mod 6) in the CRT class is the
// least Sierpinski number divisible by 3 whose covering set lies in the list
// (with period dividing L).
//
// For every cover we compute the CRT solution by Garner's algorithm and stop
// as soon as it exceeds BOUND (then the cover cannot give k < BOUND).
// Output: every cover with k < BOUND, the number of covers, and
// sum over covers of 1/prod(p), the expected density of such k.
//
// Build: gcc -O2 -o prob86_covers prob86_covers.c -lm
// Run:   ./prob86_covers L BOUND p1 p2 p3 ...        (primes with ord_p(2) | L, no 3)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;
#define MAXL 2520
#define W ((MAXL + 63) / 64)
#define MAXP 64

static int L, NP;
static uint64_t P[MAXP]; static int M[MAXP];        // prime, order of 2
static uint64_t *classmask[MAXP];                    // classmask[i][c*W .. ] : points n == c (mod M[i])
static uint64_t invpow2[MAXP][MAXL];                 // invpow2[i][a] = -2^(-a) mod P[i]  (k residue)
static int cls[MAXP];                                // assigned class or -1
static uint64_t forb[MAXP][W];                       // forbidden points per prime
static double BOUND;
static uint64_t KBOUND;
static long long ncov = 0, nsmall = 0; static double density = 0; static long long nodes = 0;
static uint64_t bestk = ~0ULL;

static inline int popc(const uint64_t *m) { int s = 0; for (int w = 0; w < W; w++) s += __builtin_popcountll(m[w]); return s; }

static uint64_t powmod(uint64_t b, uint64_t e, uint64_t m) {
    uint64_t r = 1; b %= m;
    while (e) { if (e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; }
    return r;
}
static uint64_t inv(uint64_t a, uint64_t m) { // modular inverse via Fermat for prime m
    return powmod(a % m, m - 2, m);
}

// Garner: k == r_i (mod P[i]) for assigned i, plus k == 1 (mod 2), k == 0 (mod 3).
// Returns k if k < KBOUND, else 0 (meaning "too big").  Also sets *prodlog.
static uint64_t garner(double *prodlog) {
    // moduli list
    uint64_t mods[MAXP + 2]; uint64_t res[MAXP + 2]; int n = 0;
    mods[n] = 2; res[n] = 1; n++;
    mods[n] = 3; res[n] = 0; n++;
    *prodlog = log(6.0);
    for (int i = 0; i < NP; i++) if (cls[i] >= 0) { mods[n] = P[i]; res[n] = invpow2[i][cls[i]]; n++; *prodlog += log((double)P[i]); }
    // mixed radix: k = v0 + v1*m0 + v2*m0*m1 + ...
    u128 k = 0, prod = 1;
    for (int i = 0; i < n; i++) {
        // v_i = (res_i - k) * inv(prod, m_i) mod m_i
        uint64_t mi = mods[i];
        uint64_t kk = (uint64_t)(k % mi);
        uint64_t t = (res[i] + mi - kk) % mi;
        uint64_t pr = (uint64_t)(prod % mi);
        uint64_t v;
        if (mi == 2) v = t % 2; else if (mi == 3) v = (t * inv(pr, 3)) % 3; else v = (u128)t * inv(pr, mi) % mi;
        if (v) {
            if (prod > (u128)KBOUND) return 0;      // k >= prod > bound
            k += (u128)v * prod;
            if (k >= (u128)KBOUND) return 0;
        }
        prod *= mi;
        if (prod > ((u128)1 << 120)) { // product too large to continue safely; if all later v are 0, k stays; otherwise too big.
            // remaining congruences must be satisfied by the current k
            for (int j = i + 1; j < n; j++) if ((uint64_t)(k % mods[j]) != res[j]) return 0;
            return (uint64_t)k;
        }
    }
    return (uint64_t)k;
}

static void leaf(void) {
    ncov++;
    double pl; uint64_t k = garner(&pl);
    density += exp(-pl) * 6.0;   // 6/prod: density of k == 3 mod 6 (per unit k) is 1/(6 prod) ... keep 6/prod for E[#k<X] = X * density / 36? see below
    if (k) {
        nsmall++;
        if (k < bestk) bestk = k;
        printf("K %llu :", (unsigned long long)k);
        for (int i = 0; i < NP; i++) if (cls[i] >= 0) printf(" %llu^%d", (unsigned long long)P[i], cls[i]);
        printf("\n"); fflush(stdout);
    }
}

static void dfs(uint64_t *cov, int uncovered) {
    nodes++;
    if (uncovered == 0) { leaf(); return; }
    // smallest uncovered point
    int n = -1;
    for (int w = 0; w < W; w++) { uint64_t x = ~cov[w]; if (w == W - 1 && (L & 63)) x &= (1ULL << (L & 63)) - 1; if (x) { n = w * 64 + __builtin_ctzll(x); break; } }
    // capacity pruning: remaining primes can cover at most sum L/m points
    int cap = 0; for (int i = 0; i < NP; i++) if (cls[i] < 0) cap += L / M[i];
    if (cap < uncovered) return;
    uint64_t newcov[W];
    for (int i = 0; i < NP; i++) {
        if (cls[i] >= 0) continue;
        int c = n % M[i];
        const uint64_t *cm = classmask[i] + (size_t)c * W;
        int ok = 1; for (int w = 0; w < W; w++) if (cm[w] & forb[i][w]) { ok = 0; break; }
        if (!ok) continue;
        int newunc = uncovered;
        for (int w = 0; w < W; w++) { newcov[w] = cov[w] | cm[w]; newunc -= __builtin_popcountll(cm[w] & ~cov[w]); }
        cls[i] = c;
        // CRT pruning at every node: the least k of the partial class is a lower bound for every completion
        { double pl; if (garner(&pl) == 0) { cls[i] = -1; continue; } }
        // primes before i (unassigned) must never cover n
        for (int j = 0; j < i; j++) if (cls[j] < 0) forb[j][n >> 6] |= 1ULL << (n & 63);
        dfs(newcov, newunc);
        for (int j = 0; j < i; j++) if (cls[j] < 0) forb[j][n >> 6] &= ~(1ULL << (n & 63));
        cls[i] = -1;
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: L BOUND p1 p2 ...\n"); return 1; }
    L = atoi(argv[1]); BOUND = atof(argv[2]); KBOUND = (uint64_t)BOUND;
    NP = argc - 3;
    for (int i = 0; i < NP; i++) {
        P[i] = strtoull(argv[3 + i], 0, 10);
        // order of 2 mod p
        int m = 1; uint64_t x = 2 % P[i]; while (x != 1) { x = x * 2 % P[i]; m++; if (m > MAXL) { fprintf(stderr, "order too large for %llu\n", (unsigned long long)P[i]); return 1; } }
        M[i] = m;
        if (L % m) { fprintf(stderr, "order %d of %llu does not divide L\n", m, (unsigned long long)P[i]); return 1; }
        classmask[i] = calloc((size_t)m * W, 8);
        for (int c = 0; c < m; c++) for (int n = c; n < L; n += m) classmask[i][(size_t)c * W + (n >> 6)] |= 1ULL << (n & 63);
        uint64_t p = P[i]; uint64_t i2 = inv(2, p); uint64_t t = 1;  // t = 2^(-a)
        for (int a = 0; a < m; a++) { invpow2[i][a] = (p - t) % p; t = (u128)t * i2 % p; }
        cls[i] = -1; memset(forb[i], 0, sizeof forb[i]);
    }
    uint64_t cov[W]; memset(cov, 0, sizeof cov);
    dfs(cov, L);
    // density = sum over covers T of 1/M_T.  Each cover contributes one class
    // mod 6 M_T, so about X/(6 M_T) values k < X; E(X) = X * density / 6, and
    // the least k is expected near 6/density.
    printf("DONE L=%d primes=%d nodes=%lld covers=%lld covers_with_k_below_%.3g=%lld bestk=%llu sum_1overM=%.6e expected_k_min=%.4g\n",
           L, NP, nodes, ncov, BOUND, nsmall, (unsigned long long)bestk, density, 6.0 / density);
    return 0;
}

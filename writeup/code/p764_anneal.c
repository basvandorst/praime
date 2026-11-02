/* Puzzle 764 (writeup/764.md): heuristic search for rows of N distinct integers whose neighbour gcds are 1..N-1,
 * with small sum. The search runs over orders of the labels 1..N-1; an order is scored by a greedy value assignment:
 * cell 0 = smallest unused multiple of g_1; cell i = smallest unused multiple of lcm(g_i, g_{i+1}) whose gcd with the
 * previous cell is exactly g_i; last cell likewise with g_{N-1}. (For consecutive labels a, b, c the cells need
 * a/gcd(a,b) coprime to c/gcd(b,c), or no multiplier works: such orders are infeasible and get a penalty.)
 * Simulated annealing with segment reversal, swap and insertion moves. Seeds: the row given in a file (its label
 * order), or random orders.
 * Usage: ./p764_anneal N iterations seed [start_row.txt]  -> prints the best row found and its sum
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV 4000000
static int N, K;
static unsigned long long rs = 0x9E3779B97F4A7C15ULL;
static inline unsigned long long rng(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rng() >> 11) * (1.0 / 9007199254740992.0); }
static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static inline long lcm(int a, int b) { return (long)a / gcd(a, b) * b; }
static unsigned char *used;
static long stamp_used[1]; /* unused */

static long evaluate(const int *G, int *V) {       /* the sum, plus 100000 for every cell that cannot be set */
    static int touched[512]; int nt = 0;
    long sum = 0; int bad = 0;
    for (int i = 0; i <= K; i++) {
        long L; int need = -1, prev = (i > 0) ? V[i - 1] : 0, g = (i > 0) ? G[i - 1] : G[0];
        if (i == 0) L = G[0]; else if (i == K) L = G[K - 1]; else L = lcm(G[i - 1], G[i]);
        if (i > 0) need = G[i - 1];
        /* the next cell must be able to reach an exact gcd G[i]: the cofactor w / G[i] must avoid the primes of
           G[i+1] / gcd(G[i], G[i+1]) */
        int nextc = (i < K) ? G[i] / gcd(G[i], (i + 1 < K) ? G[i + 1] : 1) : 1;
        if (i < K && i + 1 < K) nextc = G[i + 1] / gcd(G[i], G[i + 1]); else nextc = 1;
        int found = 0;
        for (int t = 1; t <= 60; t++) {
            long w = L * t;
            if (w >= MAXV) break;
            if (used[w]) continue;
            if (i > 0 && gcd(prev, (int)w) != need) continue;
            if (i < K && gcd((int)(w / G[i]), nextc) != 1) continue;
            V[i] = (int)w; used[w] = 1; touched[nt++] = (int)w; sum += w; found = 1; break;
        }
        if (!found) {                                 /* keep going with a plain unused multiple, and count a failure */
            bad++;
            for (int t = 1; ; t++) { long w = L * t; if (w < MAXV && !used[w]) { V[i] = (int)w; used[w] = 1; touched[nt++] = (int)w; sum += w; break; } if (w >= MAXV) { V[i] = (int)L; sum += L; break; } }
        }
    }
    for (int j = 0; j < nt; j++) used[touched[j]] = 0;
    return sum + 100000L * bad;
}

int main(int argc, char **argv) {
    N = atoi(argv[1]); K = N - 1; long iters = atol(argv[2]); rs = 0x9E3779B97F4A7C15ULL ^ ((unsigned long long)atol(argv[3]) * 0xD1B54A32D192ED03ULL + 0x2545F4914F6CDD1DULL); if (!rs) rs = 1;
    used = calloc(MAXV, 1);
    int *G = malloc(sizeof(int) * K), *best = malloc(sizeof(int) * K), *cand = malloc(sizeof(int) * K);
    int *V = malloc(sizeof(int) * N);
    for (int i = 0; i < K; i++) G[i] = i + 1;
    if (argc > 4) {                                   /* seed with the label order of a given row */
        FILE *f = fopen(argv[4], "r"); int *row = malloc(sizeof(int) * N);
        for (int i = 0; i < N; i++) if (fscanf(f, "%d", &row[i]) != 1) return 1;
        fclose(f);
        for (int i = 0; i < K; i++) G[i] = gcd(row[i], row[i + 1]);
    } else for (int i = K - 1; i > 0; i--) { int j = rng() % (i + 1); int t = G[i]; G[i] = G[j]; G[j] = t; }
    long cur = evaluate(G, V), bestc = cur; memcpy(best, G, sizeof(int) * K);
    double T0 = 0.02 * (cur % 100000), T1 = 0.0005 * (cur % 100000); if (T0 < 1) { T0 = 20; T1 = 0.5; }
    for (long it = 0; it < iters; it++) {
        double T = T0 * pow(T1 / T0, (double)it / iters);
        memcpy(cand, G, sizeof(int) * K);
        int a = rng() % K, b = rng() % K, mv = rng() % 3;
        if (a == b) continue;
        if (mv == 0) { int t = cand[a]; cand[a] = cand[b]; cand[b] = t; }
        else if (mv == 1) { if (a > b) { int t = a; a = b; b = t; } for (int i = a, j = b; i < j; i++, j--) { int t = cand[i]; cand[i] = cand[j]; cand[j] = t; } }
        else { int x = cand[a]; if (a < b) { memmove(cand + a, cand + a + 1, sizeof(int) * (b - a)); cand[b] = x; } else { memmove(cand + b + 1, cand + b, sizeof(int) * (a - b)); cand[b] = x; } }
        long c = evaluate(cand, V);
        if (c <= cur || urand() < exp((cur - c) / T)) {
            memcpy(G, cand, sizeof(int) * K); cur = c;
            if (c < bestc) { bestc = c; memcpy(best, G, sizeof(int) * K); }
        }
    }
    evaluate(best, V);
    printf("N = %d: best sum %ld%s\n", N, bestc, bestc >= 100000 ? " (INFEASIBLE)" : "");
    for (int i = 0; i < N; i++) printf("%d%c", V[i], i == N - 1 ? '\n' : ' ');
    return 0;
}

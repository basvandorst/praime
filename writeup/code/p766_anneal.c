/* Puzzles 766, 767 (writeup/764.md): heuristic search for an N x N grid or N x N x N cube of distinct integers whose
 * neighbour gcds are exactly 1..E, with small sum. State: an assignment of the labels 1..E to the edges. Score: cells
 * are filled in a fixed order; a cell gets the smallest unused multiple of the lcm of its edge labels such that the
 * gcd with every already filled neighbour is exactly the edge label, and its cofactor avoids the primes that a not
 * yet filled neighbour will need (otherwise that neighbour can never reach an exact gcd). Cells that cannot be filled
 * count as failures (penalty). Simulated annealing with swaps of two edge labels.
 * Usage: ./p766_anneal grid|cube N iterations seed [start.txt]   (start.txt: a known arrangement, row by row)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV 40000000
static int n, E;
static int eu[2048], ev[2048];               /* edge endpoints */
static int deg[1024], adj[1024][6], adje[1024][6];
static unsigned char *used;
static unsigned long long rs;
static inline unsigned long long rng(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rng() >> 11) * (1.0 / 9007199254740992.0); }
static long gcdl(long a, long b) { while (b) { long t = a % b; a = b; b = t; } return a; }
static long lcml(long a, long b) { return a / gcdl(a, b) * b; }

static long evaluate(const int *lab, long *V) {
    static int touched[1024]; int nt = 0; long sum = 0; int bad = 0;
    for (int c = 0; c < n; c++) {
        long L = 1;
        for (int j = 0; j < deg[c]; j++) L = lcml(L, lab[adje[c][j]]);
        if (L >= MAXV) { bad++; V[c] = L; sum += L; continue; }
        int found = 0;
        for (long w = L; w < MAXV && w <= 200L * L; w += L) {
            if (used[w]) continue;
            int ok = 1;
            for (int j = 0; j < deg[c] && ok; j++) {
                int u = adj[c][j], g = lab[adje[c][j]];
                if (u < c) { if (gcdl(w, V[u]) != g) ok = 0; }
                else {                                  /* future neighbour: its lcm over g must be coprime to w / g */
                    long Lu = 1;
                    for (int k = 0; k < deg[u]; k++) Lu = lcml(Lu, lab[adje[u][k]]);
                    if (gcdl(w / g, Lu / g) != 1) ok = 0;
                }
            }
            if (!ok) continue;
            V[c] = w; used[w] = 1; touched[nt++] = (int)w; sum += w; found = 1; break;
        }
        if (!found) { bad++; long w = L; while (w < MAXV && used[w]) w += L; V[c] = w; if (w < MAXV) { used[w] = 1; touched[nt++] = (int)w; } sum += w; }
    }
    for (int j = 0; j < nt; j++) used[touched[j]] = 0;
    return sum + 10000000L * bad;
}

int main(int argc, char **argv) {
    const char *kind = argv[1]; int N = atoi(argv[2]); long iters = atol(argv[3]);
    rs = 0x9E3779B97F4A7C15ULL ^ ((unsigned long long)atol(argv[4]) * 0xD1B54A32D192ED03ULL + 0x2545F4914F6CDD1DULL); if (!rs) rs = 1;
    int dim = strcmp(kind, "grid") == 0 ? 2 : 3;
    int na = N, nb = N, nc = dim == 3 ? N : 1; n = na * nb * nc; E = 0;
    #define ID(a, b, c) (((a) * nb + (b)) * nc + (c))
    for (int a = 0; a < na; a++) for (int b = 0; b < nb; b++) for (int c = 0; c < nc; c++) {
        int x = ID(a, b, c);
        if (c + 1 < nc) { eu[E] = x; ev[E] = ID(a, b, c + 1); E++; }
        if (b + 1 < nb) { eu[E] = x; ev[E] = ID(a, b + 1, c); E++; }
        if (a + 1 < na) { eu[E] = x; ev[E] = ID(a + 1, b, c); E++; }
    }
    for (int e = 0; e < E; e++) { adj[eu[e]][deg[eu[e]]] = ev[e]; adje[eu[e]][deg[eu[e]]++] = e; adj[ev[e]][deg[ev[e]]] = eu[e]; adje[ev[e]][deg[ev[e]]++] = e; }
    used = calloc(MAXV, 1);
    int *lab = malloc(sizeof(int) * E), *best = malloc(sizeof(int) * E); long *V = malloc(sizeof(long) * n);
    for (int e = 0; e < E; e++) lab[e] = e + 1;
    if (argc > 5) {
        FILE *f = fopen(argv[5], "r"); long *row = malloc(sizeof(long) * n);
        for (int i = 0; i < n; i++) if (fscanf(f, "%ld", &row[i]) != 1) return 1;
        fclose(f);
        for (int e = 0; e < E; e++) lab[e] = (int)gcdl(row[eu[e]], row[ev[e]]);
    } else for (int i = E - 1; i > 0; i--) { int j = rng() % (i + 1); int t = lab[i]; lab[i] = lab[j]; lab[j] = t; }
    long cur = evaluate(lab, V), bestc = cur; memcpy(best, lab, sizeof(int) * E);
    double base = (double)(cur % 10000000L); if (base < 10) base = 1000;
    double T0 = 0.02 * base, T1 = 0.0003 * base;
    for (long it = 0; it < iters; it++) {
        double T = T0 * pow(T1 / T0, (double)it / iters);
        int a = rng() % E, b = rng() % E; if (a == b) continue;
        int t = lab[a]; lab[a] = lab[b]; lab[b] = t;
        long c = evaluate(lab, V);
        if (c <= cur || urand() < exp((cur - c) / T)) { cur = c; if (c < bestc) { bestc = c; memcpy(best, lab, sizeof(int) * E); } }
        else { t = lab[a]; lab[a] = lab[b]; lab[b] = t; }
    }
    evaluate(best, V);
    printf("%s N = %d: best sum %ld%s\n", kind, N, bestc, bestc >= 10000000L ? " (INFEASIBLE)" : "");
    for (int i = 0; i < n; i++) printf("%ld%c", V[i], (i % N == N - 1) ? '\n' : ' ');
    return 0;
}

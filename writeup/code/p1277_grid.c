/* Puzzle 1277 (writeup/1277.md): arrange 1..n^2 in an n x n grid so that few edge-neighbours are coprime.
 *
 * Start (structured): 1 and the primes p > n^2/2 are coprime to every other entry ("lonely"), so they go on the
 * outer ring (degree 3) and, for the rest, on the second ring. The primes n^2/4 < p <= n^2/2 have only one or two
 * multiples; they fill the rest of the second ring with 2p directly inside them. Everything else is filled in by
 * smallest prime factor (evens first, then multiples of 3, 5, ...).
 * Then local improvement: simulated annealing with swaps of two cells, the count of coprime edges as energy.
 *
 * The outer two rings and the 2p cells (the skeleton) stay fixed during the annealing.
 * Usage: ./p1277_grid n iterations seed out.txt [T0 [start.txt mask.txt]]   (prints the best count; writes the grid)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int n, NN;
static int *g;                       /* g[r*n+c] */
static char *frozen;                 /* skeleton cells that the annealing may not move */
static unsigned long long rng_s;
static inline unsigned long long rng(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }
static inline double urand(void) { return (rng() >> 11) * (1.0 / 9007199254740992.0); }
static inline int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static inline int cop(int a, int b) { return gcd(a, b) == 1; }

static int cell_cost(int idx, int v) {          /* coprime edges of value v placed at idx (current neighbours) */
    int r = idx / n, c = idx % n, s = 0;
    if (r > 0) s += cop(v, g[idx - n]);
    if (r < n - 1) s += cop(v, g[idx + n]);
    if (c > 0) s += cop(v, g[idx - 1]);
    if (c < n - 1) s += cop(v, g[idx + 1]);
    return s;
}
static long total(void) {
    long s = 0;
    for (int r = 0; r < n; r++) for (int c = 0; c < n; c++) {
        if (c + 1 < n) s += cop(g[r*n+c], g[r*n+c+1]);
        if (r + 1 < n) s += cop(g[r*n+c], g[(r+1)*n+c]);
    }
    return s;
}
static int delta_swap(int a, int b) {
    int va = g[a], vb = g[b];
    int before = cell_cost(a, va) + cell_cost(b, vb);
    g[a] = vb; g[b] = va;
    int after = cell_cost(a, vb) + cell_cost(b, va);
    g[a] = va; g[b] = vb;
    int adj = (abs(a - b) == n) || (abs(a - b) == 1 && a / n == b / n);
    /* an a-b edge is counted twice in both sums and does not change: no correction needed */
    (void)adj;
    return after - before;
}

int main(int argc, char **argv) {
    n = atoi(argv[1]); long iters = atol(argv[2]); rng_s = 88172645463325252ULL ^ (unsigned long long)atol(argv[3]);
    NN = n * n; g = malloc(sizeof(int) * NN); frozen = calloc(NN, 1);
    char *isp = calloc(NN + 1, 1); for (int i = 2; i <= NN; i++) isp[i] = 1;
    for (int i = 2; (long)i * i <= NN; i++) if (isp[i]) for (int j = i * i; j <= NN; j += i) isp[j] = 0;
    int *used = calloc(NN + 1, sizeof(int)); for (int i = 0; i < NN; i++) g[i] = 0;
    /* lonely numbers: 1 and primes > NN/2, largest first */
    int *lon = malloc(sizeof(int) * NN), nl = 0;
    for (int v = NN; v >= 2; v--) if (isp[v] && 2 * v > NN) lon[nl++] = v;
    lon[nl++] = 1;
    /* semi-lonely: primes in (NN/4, NN/2], largest first */
    int *semi = malloc(sizeof(int) * NN), ns = 0;
    for (int v = NN / 2; v >= 2; v--) if (isp[v] && 4 * v > NN) semi[ns++] = v;
    /* ring cells in cyclic order */
    int *ring = malloc(sizeof(int) * 4 * n), il = 0, is = 0;
    for (int k = 0; k < 2 && n - 2 * k >= 2; k++) {
        int m = 0, lo = k, hi = n - 1 - k;
        for (int c = lo; c < hi; c++) ring[m++] = lo * n + c;
        for (int r = lo; r < hi; r++) ring[m++] = r * n + hi;
        for (int c = hi; c > lo; c--) ring[m++] = hi * n + c;
        for (int r = hi; r > lo; r--) ring[m++] = r * n + lo;
        for (int i = 0; i < m; i++) {
            int idx = ring[i], r = idx / n, c = idx % n;
            int corner = (r == lo || r == hi) && (c == lo || c == hi);
            int inner = (k == 1) ? -1 : 0;
            if (k == 1 && !corner) {    /* inward neighbour */
                if (r == lo) inner = idx + n; else if (r == hi) inner = idx - n; else if (c == lo) inner = idx + 1; else inner = idx - 1;
            }
            int remaining_cells = m - i;
            if (il < nl && (k == 0 || corner || nl - il >= remaining_cells || is >= ns)) { g[idx] = lon[il]; used[lon[il++]] = 1; frozen[idx] = 1; }
            else if (k == 1 && is < ns && 2 * semi[is] <= NN && inner >= 0 && !g[inner]) {
                g[idx] = semi[is]; used[semi[is]] = 1; g[inner] = 2 * semi[is]; used[2 * semi[is]] = 1; is++; frozen[idx] = frozen[inner] = 1;
            }
        }
    }
    /* fill the rest by smallest prime factor: evens, then 3, 5, 7, ... (primes after their multiples) */
    int *spf = calloc(NN + 1, sizeof(int));
    for (int v = 2; v <= NN; v++) if (!spf[v]) for (int j = v; j <= NN; j += v) if (!spf[j]) spf[j] = v;
    int *rest = malloc(sizeof(int) * NN), nr = 0;
    for (int p = 2; p <= NN; p++) if (isp[p]) for (int v = p; v <= NN; v += p) if (!used[v] && spf[v] == p) { rest[nr++] = v; used[v] = 1; }
    for (int v = 1; v <= NN; v++) if (!used[v]) { rest[nr++] = v; used[v] = 1; }
    int k = 0;
    for (int r = 0; r < n; r++) for (int c = 0; c < n; c++) {
        int cc = (r % 2) ? n - 1 - c : c;            /* snake order keeps runs adjacent */
        if (!g[r * n + cc]) g[r * n + cc] = rest[k++];
    }
    if (argc > 7) {                      /* start grid and frozen mask from files (p1277_build.py) */
        FILE *fs = fopen(argv[6], "r"), *fm = fopen(argv[7], "r");
        for (int i = 0; i < NN; i++) { int m; if (fscanf(fs, "%d", &g[i]) != 1 || fscanf(fm, "%d", &m) != 1) return 1; frozen[i] = (char)m; }
        fclose(fs); fclose(fm);
    }
    long cur = total(), best = cur;
    fprintf(stderr, "start %ld\n", cur);
    int *bestg = malloc(sizeof(int) * NN); memcpy(bestg, g, sizeof(int) * NN);
    double T0 = argc > 5 ? atof(argv[5]) : 0.6, T1 = 0.03;
    for (long it = 0; it < iters; it++) {
        double T = T0 * pow(T1 / T0, (double)it / iters);
        int a = rng() % NN, b;
        if (urand() < 0.5) b = rng() % NN;
        else { int r = a / n + (int)(rng() % 11) - 5, c = a % n + (int)(rng() % 11) - 5;
               if (r < 0 || r >= n || c < 0 || c >= n) continue; b = r * n + c; }
        if (a == b || frozen[a] || frozen[b]) continue;
        int d = delta_swap(a, b);
        if (d <= 0 || urand() < exp(-d / T)) {
            int t = g[a]; g[a] = g[b]; g[b] = t; cur += d;
            if (cur < best) { best = cur; memcpy(bestg, g, sizeof(int) * NN); }
        }
    }
    memcpy(g, bestg, sizeof(int) * NN);
    printf("n = %d: %ld coprime neighbour pairs out of %d (check %ld)\n", n, best, 2 * n * (n - 1), total());
    FILE *f = fopen(argv[4], "w");
    for (int r = 0; r < n; r++) { for (int c = 0; c < n; c++) fprintf(f, c ? " %d" : "%d", g[r * n + c]); fprintf(f, "\n"); }
    fclose(f);
    return 0;
}

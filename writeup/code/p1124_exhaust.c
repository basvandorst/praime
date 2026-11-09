/* Puzzles 1123 / 1124.  Triangle with R rows (R = 4: tetractys with 0..9, R = 5: pentactys with 0..14).
   Along every straight line (rows and the two diagonal directions) every contiguous segment is read forwards
   and backwards; the numbers are concatenated in decimal; count the DISTINCT primes.  Exhaustive DFS over all
   arrangements, symmetry broken by requiring the top corner value < the other two corners and left corner <
   right corner (the dihedral group acts as S3 on the corners).  Prime ids by value (different sequences can give
   the same number, e.g. 1|13 and 11|3), looked up in a table over value sequences of length <= R.
   Prints every arrangement with >= T distinct primes and the maximum.  Usage: p1124_exhaust R T [threads]
   gcc -O2 -fopenmp -o p1124_exhaust p1124_exhaust.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
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
static int R, NC, NV, T;                  /* rows, cells, values (= NC), threshold */
static int cell[5][5];                    /* cell index of (row, i) */
#define MAXSEG 512
static int segs[MAXSEG][6], seglen[MAXSEG], nseg = 0;   /* ordered cell lists (one per reading direction) */
static int *pid[6];                       /* pid[len][mixed-radix index over values] = prime id or -1 */
static int npr = 0;                       /* number of distinct primes */
static uint64_t *prval; static int prcap = 0;
static int find_or_add(uint64_t v) {      /* linear-probing hash on value */
    static uint64_t *hk = 0; static int *hv = 0; static int hs = 1 << 20;
    if (!hk) { hk = calloc(hs, 8); hv = malloc(hs * sizeof(int)); }
    uint64_t h = (v * 11400714819323198485ull) >> 44;
    while (hk[h]) { if (hk[h] == v) return hv[h]; h = (h + 1) & (hs - 1); }
    hk[h] = v; hv[h] = npr; return npr++;
}
static void build_tables(void) {
    int idx[6]; 
    for (int len = 1; len <= R; len++) {
        long sz = 1; for (int i = 0; i < len; i++) sz *= NV;
        pid[len] = malloc(sz * sizeof(int));
        for (long code = 0; code < sz; code++) {
            long c = code; int ok = 1; uint64_t v = 0;
            for (int i = len - 1; i >= 0; i--) { idx[i] = c % NV; c /= NV; }
            for (int i = 0; i < len && ok; i++) for (int j = 0; j < i; j++) if (idx[i] == idx[j]) { ok = 0; break; }
            if (!ok) { pid[len][code] = -1; continue; }
            for (int i = 0; i < len; i++) { int x = idx[i]; v = v * (x >= 10 ? 100 : 10) + x; }
            pid[len][code] = isprime64(v) ? find_or_add(v) : -1;
        }
    }
}
static void add_line(int *cells, int len) {           /* all contiguous segments, both directions */
    for (int a = 0; a < len; a++) for (int b = a; b < len; b++) {
        int L = b - a + 1;
        for (int i = 0; i < L; i++) segs[nseg][i] = cells[a + i]; seglen[nseg] = L; nseg++;
        if (L > 1) { for (int i = 0; i < L; i++) segs[nseg][i] = cells[b - i]; seglen[nseg] = L; nseg++; }
    }
}
static int order[25], depthof[25];        /* fill order and depth of each cell */
static int *compl_at[25]; static int ncompl[25];   /* segments completing at each depth */
static long long found = 0; static int best = 0;
static void dfs(int d, int *val, int *used, int *cnt, int distinct) {
    if (d == NC) {
        #pragma omp critical
        { if (distinct > best) best = distinct;
          if (distinct >= T) { found++; printf("%d primes:", distinct); for (int i = 0; i < NC; i++) printf(" %d", val[i]); printf("\n"); fflush(stdout); } }
        return;
    }
    int c = order[d];
    for (int v = 0; v < NV; v++) {
        if (used[v]) continue;
        /* symmetry: top corner (cell 0) < left corner cell[R-1][0] < right corner cell[R-1][R-1] */
        if (c == cell[R-1][0] && v < val[0]) continue;
        if (c == cell[R-1][R-1] && (v < val[0] || v < val[cell[R-1][0]])) continue;
        val[c] = v; used[v] = 1;
        int added = 0; int touched[64]; int nt = 0;
        for (int k = 0; k < ncompl[d]; k++) {
            int s = compl_at[d][k]; long code = 0;
            for (int i = 0; i < seglen[s]; i++) code = code * NV + val[segs[s][i]];
            int id = pid[seglen[s]][code];
            if (id >= 0) { if (cnt[id]++ == 0) added++; touched[nt++] = id; }
        }
        dfs(d + 1, val, used, cnt, distinct + added);
        for (int k = 0; k < nt; k++) cnt[touched[k]]--;
        used[v] = 0;
    }
    val[c] = -1;
}
int main(int argc, char **argv) {
    R = atoi(argv[1]); T = atoi(argv[2]); int th = argc > 3 ? atoi(argv[3]) : 2;
    NC = R * (R + 1) / 2; NV = NC;
    int k = 0; for (int r = 0; r < R; r++) for (int i = 0; i <= r; i++) cell[r][i] = k++;
    build_tables();
    int line[6];
    for (int r = 0; r < R; r++) { for (int i = 0; i <= r; i++) line[i] = cell[r][i]; add_line(line, r + 1); }           /* rows */
    for (int i = 0; i < R; i++) { int n = 0; for (int r = i; r < R; r++) line[n++] = cell[r][i]; add_line(line, n); }   /* left edges */
    for (int j = 0; j < R; j++) { int n = 0; for (int r = j; r < R; r++) line[n++] = cell[r][r - j]; add_line(line, n); } /* right edges */
    /* fill order: top corner, then the two bottom corners early (symmetry pruning), then the rest row by row */
    int o = 0; order[o++] = cell[0][0]; order[o++] = cell[R-1][0]; order[o++] = cell[R-1][R-1];
    for (int c = 0; c < NC; c++) { int seen = 0; for (int i = 0; i < o; i++) if (order[i] == c) seen = 1; if (!seen) order[o++] = c; }
    for (int d = 0; d < NC; d++) depthof[order[d]] = d;
    for (int d = 0; d < NC; d++) { compl_at[d] = malloc(nseg * sizeof(int)); ncompl[d] = 0; }
    for (int s = 0; s < nseg; s++) { int m = 0; for (int i = 0; i < seglen[s]; i++) if (depthof[segs[s][i]] > m) m = depthof[segs[s][i]]; compl_at[m][ncompl[m]++] = s; }
    fprintf(stderr, "R=%d: %d cells, %d segment readings, %d distinct primes available\n", R, NC, nseg, npr);
    if (argc > 4 && argc >= 4 + NC) {            /* evaluation mode: p1124_exhaust R T threads v0 v1 ... (row by row) */
        int val[25], *cnt = calloc(npr + 1, sizeof(int)); int distinct = 0;
        for (int i = 0; i < NC; i++) val[i] = atoi(argv[4 + i]);
        for (int s = 0; s < nseg; s++) { long code = 0; for (int i = 0; i < seglen[s]; i++) code = code * NV + val[segs[s][i]];
            int id = pid[seglen[s]][code]; if (id >= 0 && cnt[id]++ == 0) distinct++; }
        printf("arrangement has %d distinct primes\n", distinct); return 0;
    }
    omp_set_num_threads(th);
    /* parallel over the value of the top corner and the left corner */
    #pragma omp parallel for schedule(dynamic) collapse(2)
    for (int v0 = 0; v0 < NV; v0++) for (int v1 = 0; v1 < NV; v1++) {
        if (v1 <= v0) continue;
        int val[25], used[25], *cnt = calloc(npr + 1, sizeof(int));
        for (int i = 0; i < 25; i++) { val[i] = -1; used[i] = 0; }
        val[order[0]] = v0; used[v0] = 1; val[order[1]] = v1; used[v1] = 1;
        int distinct = 0;
        for (int d = 0; d < 2; d++) for (int k2 = 0; k2 < ncompl[d]; k2++) {
            int s = compl_at[d][k2]; long code = 0;
            for (int i = 0; i < seglen[s]; i++) code = code * NV + val[segs[s][i]];
            int id = pid[seglen[s]][code]; if (id >= 0 && cnt[id]++ == 0) distinct++;
        }
        dfs(2, val, used, cnt, distinct);
        free(cnt);
    }
    printf("maximum %d, arrangements with >= %d primes: %lld\n", best, T, found);
    return 0;
}

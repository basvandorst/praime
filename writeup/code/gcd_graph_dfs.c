/* Puzzles 764, 766, 767 (writeup/764.md): independent check of the CP-SAT results of gcd_graph_cpsat.py.
 * Exhaustive depth-first search: vertices get distinct values in [1, B]; whenever both ends of an edge are set, the
 * gcd must be a label in 1..E not used yet; the sum of all values must be <= S. Prunes when the current sum plus the
 * smallest unused values for the remaining vertices exceeds S. Prints the first solution found, or "none".
 * Usage: ./gcd_graph_dfs path|grid|cube N S        (B = S - n(n-1)/2, see the Python model)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n, E, B, S;
static int adj[64][8], deg[64];
static int val[64], used[4096], lab[64];
static long nodes;
static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

static int dfs(int i, int sum) {
    nodes++;
    if (i == n) return 1;
    /* lower bound: remaining n-i vertices need distinct unused values */
    int need = 0, cnt = 0, top = 0;
    for (int v = 1; v <= B && cnt < n - i; v++) if (!used[v]) { need += v; cnt++; top = v; }
    if (cnt < n - i || sum + need > S) return 0;
    for (int v = 1; v <= B; v++) {
        /* the other n-i-1 vertices need the smallest unused values other than v */
        if (used[v] || sum + v + (v <= top ? need - v : need - top) > S) continue;
        int ok = 1, got[8], k = 0;
        for (int j = 0; j < deg[i]; j++) {
            int u = adj[i][j];
            if (u < i) {
                int g = gcd(v, val[u]);
                if (g > E || lab[g]) { ok = 0; break; }
                for (int t = 0; t < k; t++) if (got[t] == g) { ok = 0; break; }
                if (!ok) break;
                got[k++] = g;
            }
        }
        if (!ok) continue;
        val[i] = v; used[v] = 1;
        for (int t = 0; t < k; t++) lab[got[t]] = 1;
        if (dfs(i + 1, sum + v)) return 1;
        for (int t = 0; t < k; t++) lab[got[t]] = 0;
        used[v] = 0;
    }
    return 0;
}

int main(int argc, char **argv) {
    const char *kind = argv[1]; int N = atoi(argv[2]); S = atoi(argv[3]);
    int id[8][8][8]; n = 0; E = 0;
    int dim = strcmp(kind, "path") == 0 ? 1 : strcmp(kind, "grid") == 0 ? 2 : 3;
    int na = dim >= 1 ? N : 1, nb = dim >= 2 ? N : 1, nc = dim >= 3 ? N : 1;
    for (int a = 0; a < na; a++) for (int b = 0; b < nb; b++) for (int c = 0; c < nc; c++) id[a][b][c] = n++;
    for (int a = 0; a < na; a++) for (int b = 0; b < nb; b++) for (int c = 0; c < nc; c++) {
        int x = id[a][b][c];
        if (a + 1 < na) { int y = id[a + 1][b][c]; adj[x][deg[x]++] = y; adj[y][deg[y]++] = x; E++; }
        if (b + 1 < nb) { int y = id[a][b + 1][c]; adj[x][deg[x]++] = y; adj[y][deg[y]++] = x; E++; }
        if (c + 1 < nc) { int y = id[a][b][c + 1]; adj[x][deg[x]++] = y; adj[y][deg[y]++] = x; E++; }
    }
    B = S - n * (n - 1) / 2;
    printf("%s N=%d: %d vertices, %d edges, sum <= %d, values <= %d\n", kind, N, n, E, S, B);
    if (dfs(0, 0)) {
        int sum = 0; printf("solution:"); for (int i = 0; i < n; i++) { printf(" %d", val[i]); sum += val[i]; }
        printf("  sum %d  (%ld nodes)\n", sum, nodes);
    } else printf("none: no assignment with sum <= %d  (%ld nodes)\n", S, nodes);
    return 0;
}

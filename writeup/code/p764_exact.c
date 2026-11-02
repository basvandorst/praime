/* Puzzle 764 (writeup/764.md): exact search for a row of N distinct positive integers whose neighbour gcds are
 * 1..N-1, with sum <= S. Two ingredients:
 *  1. Relaxed completion cost. If the cell just placed has right label g and the labels in 'mask' are still unused,
 *     the remaining cells cost at least rem[mask][g] = min over orders g1..gk of mask of
 *     lcm(g,g1) + lcm(g1,g2) + ... + lcm(g_{k-1},g_k)   (a cell with labels a, b is a multiple of lcm(a, b);
 *     the last cell, with one label, is at least that label: it is counted as lcm(g_{k-1}, g_k) >= g_k here
 *     only through its left label, so the true last cell costs lcm(g_{k-1},g_k)... see below).
 *     Precisely: cells i+1..N-1 are counted as lcm(label_left, label_right) for the interior ones and label_left
 *     for the last one; rem[0][g] = 0 is the case "cell i is the last cell".
 *  2. Depth-first search over cells: choose the next label g among the unused ones and the next value w <= B,
 *     a multiple of g, unused, with gcd(previous value, w) = g; prune when sum + w + rem[mask - g][g] > S.
 * Every value is <= B = S - N(N-1)/2 (the other values are distinct). Prints the first row found, or "none".
 * Usage: ./p764_exact N S
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int N, m, S, B;
static int *rem;                     /* rem[mask * m + (g-1)] */
static unsigned char *usedv;
static int val[64], lab[64];
static long nodes;
static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static int lcm(int a, int b) { return a / gcd(a, b) * b; }

static int dfs(int i, long mask, int sum) {        /* cell i to place; cell i-1 placed with right label lab[i-1] */
    nodes++;
    if (i == N) return 1;
    int g = lab[i - 1], prev = val[i - 1];
    if (i == N - 1) {                              /* last cell: only its left label */
        for (int w = g; w <= B && sum + w <= S; w += g)
            if (!usedv[w] && gcd(prev, w) == g) { val[i] = w; return 1; }
        return 0;
    }
    for (int h = 1; h <= m; h++) {                 /* right label of cell i */
        if (!(mask >> (h - 1) & 1)) continue;
        long nm = mask & ~(1L << (h - 1));
        int L = lcm(g, h), r = rem[nm * m + (h - 1)];
        for (int w = L; w <= B && sum + w + r <= S; w += L) {
            if (usedv[w] || gcd(prev, w) != g) continue;
            val[i] = w; lab[i] = h; usedv[w] = 1;
            if (dfs(i + 1, nm, sum + w)) return 1;
            usedv[w] = 0;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    N = atoi(argv[1]); S = atoi(argv[2]); m = N - 1; B = S - N * (N - 1) / 2;
    long full = 1L << m;
    rem = malloc(sizeof(int) * full * m);
    for (int g = 1; g <= m; g++) rem[0 * m + g - 1] = 0;
    for (long mask = 1; mask < full; mask++)
        for (int g = 1; g <= m; g++) {
            int best = 1 << 29;
            for (int h = 1; h <= m; h++) if (mask >> (h - 1) & 1) {
                int c = lcm(g, h) + rem[(mask & ~(1L << (h - 1))) * m + h - 1];
                if (c < best) best = c;
            }
            rem[mask * m + g - 1] = best;
        }
    usedv = calloc(B + 1, 1);
    printf("N = %d, sum <= %d, values <= %d\n", N, S, B);
    /* first cell: value v0 = multiple of its right label g1 */
    for (int g = 1; g <= m; g++) {
        long nm = (full - 1) & ~(1L << (g - 1));
        for (int v0 = g; v0 <= B && v0 + rem[nm * m + g - 1] <= S; v0 += g) {
            val[0] = v0; lab[0] = g; usedv[v0] = 1;
            if (dfs(1, nm, v0)) {
                int sum = 0; printf("row:"); for (int i = 0; i < N; i++) { printf(" %d", val[i]); sum += val[i]; }
                printf("  sum %d  (%ld nodes)\n", sum, nodes); return 0;
            }
            usedv[v0] = 0;
        }
    }
    printf("none: no row with sum <= %d  (%ld nodes)\n", S, nodes);
    return 0;
}

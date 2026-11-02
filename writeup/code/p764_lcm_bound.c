/* Puzzle 764 (writeup/764.md): lower bound by relaxation. A cell with edge labels a, b holds a multiple of lcm(a, b),
 * an end cell with label a holds a multiple of a. So the sum is at least
 *     min over orders (g_1, ..., g_{N-1}) of the labels 1..N-1 of  g_1 + lcm(g_1,g_2) + ... + lcm(g_{N-2},g_{N-1}) + g_{N-1},
 * a shortest Hamiltonian path, computed exactly by dynamic programming over subsets (Held-Karp).
 * The relaxation ignores that values must be distinct and that gcds must be exact, so it is a lower bound only.
 * Usage: ./p764_lcm_bound N      (N - 1 labels; N <= 24 fits in memory)
 */
#include <stdio.h>
#include <stdlib.h>
static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
int main(int argc, char **argv) {
    int N = atoi(argv[1]), m = N - 1;
    long full = 1L << m;
    int *dp = malloc(sizeof(int) * full * m);
    int *par = malloc(sizeof(int) * full * m);
    for (long i = 0; i < full * m; i++) dp[i] = 1 << 29;
    for (int g = 1; g <= m; g++) dp[(1L << (g - 1)) * m + (g - 1)] = g;
    for (long mask = 1; mask < full; mask++)
        for (int last = 0; last < m; last++) {
            int c = dp[mask * m + last];
            if (c >= 1 << 29 || !(mask >> last & 1)) continue;
            for (int nx = 0; nx < m; nx++) {
                if (mask >> nx & 1) continue;
                int a = last + 1, b = nx + 1, l = a / gcd(a, b) * b;
                long nm = mask | (1L << nx);
                if (c + l < dp[nm * m + nx]) { dp[nm * m + nx] = c + l; par[nm * m + nx] = last; }
            }
        }
    int best = 1 << 29, bl = -1;
    for (int last = 0; last < m; last++) if (dp[(full - 1) * m + last] + last + 1 < best) { best = dp[(full - 1) * m + last] + last + 1; bl = last; }
    printf("N = %d: lower bound %d  (label order:", N, best);
    int order[64], k = 0; long mask = full - 1; int last = bl;
    while (k < m) { order[k++] = last + 1; int p = par[mask * m + last]; mask &= ~(1L << last); last = p; }
    for (int i = m - 1; i >= 0; i--) printf(" %d", order[i]);
    printf(")\n");
    return 0;
}

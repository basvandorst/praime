// Puzzle 11: first N pairwise distinct consecutive gaps, first N strictly increasing / decreasing gaps.
// Records runs whose START prime lies in [LO, HI). Sieve starts at LO - OVL to rebuild the window state.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef uint64_t u64;
#define SEG (1u<<24)
#define NMAX 64
#define WMAX 4096
int main(int argc, char **argv) {
  u64 LO = strtoull(argv[1], 0, 10), HI = strtoull(argv[2], 0, 10);
  const u64 OVL = 2000000; u64 S0 = LO > OVL + 3 ? LO - OVL : 3; if (!(S0 & 1)) S0++;
  u64 END = HI + OVL;
  u64 r = (u64)sqrtl((long double)END) + 2; char *bs = calloc(r + 1, 1); u64 *bp = malloc(sizeof(u64) * r / 2); int nb = 0;
  for (u64 i = 2; i <= r; i++) if (!bs[i]) { bp[nb++] = i; for (u64 j = i * i; j <= r; j += i) bs[j] = 1; }
  u64 *nxt = malloc(sizeof(u64) * nb);
  for (int k = 1; k < nb; k++) { u64 p = bp[k]; u64 m = (S0 + p - 1) / p * p; if (m < p * p) m = p * p; if (!(m & 1)) m += p; nxt[k] = m; }
  // ring buffer of recent primes and gaps
  static u64 P[WMAX]; static int G[WMAX]; long np = 0;          // P[t % WMAX]
  static int cnt[2048]; long wi = 0;                              // distinct window over gap indices [wi, np-1)
  static u64 dfirst[NMAX], ifirst[NMAX], rfirst[NMAX];
  long inc = 0, dec = 0;                                          // current run lengths (gaps) ending at last gap
  unsigned char *seg = malloc(SEG);
  for (u64 lo = S0; lo < END; lo += 2ull * SEG) {
    memset(seg, 1, SEG); u64 hi = lo + 2ull * SEG;
    for (int k = 1; k < nb; k++) { u64 p = bp[k]; if (p * p >= hi) break; u64 s = nxt[k];
      for (; s < hi; s += 2 * p) seg[(s - lo) >> 1] = 0; nxt[k] = s; }
    for (u64 i = 0; i < SEG; i++) if (seg[i]) {
      u64 n = lo + 2 * i; if (n >= END) break;
      P[np % WMAX] = n;
      if (np >= 1) {
        long t = np - 1;                                          // gap index t: P[t] -> P[t+1]
        int g = (int)(n - P[t % WMAX]); if (g > 2047) g = 2047;
        G[t % WMAX] = g;
        // distinct window: gaps wi..t
        while (cnt[g]) { cnt[G[wi % WMAX]]--; wi++; }
        cnt[g]++;
        long L = t - wi + 1;                                      // number of distinct gaps ending at t
        for (long N = 2; N <= L && N < NMAX; N++) {
          u64 st = P[(t - N + 1) % WMAX];
          if (!dfirst[N] && st >= LO && st < HI) dfirst[N] = st;
        }
        if (t >= 1) {
          int gp = G[(t - 1) % WMAX];
          inc = (g > gp) ? inc + 1 : 1; dec = (g < gp) ? dec + 1 : 1;
        } else { inc = dec = 1; }
        for (long N = 2; N <= inc && N < NMAX; N++) { u64 st = P[(t - N + 1) % WMAX]; if (!ifirst[N] && st >= LO && st < HI) ifirst[N] = st; }
        for (long N = 2; N <= dec && N < NMAX; N++) { u64 st = P[(t - N + 1) % WMAX]; if (!rfirst[N] && st >= LO && st < HI) rfirst[N] = st; }
        if (t - wi + 1 > WMAX - 8) { fprintf(stderr, "window overflow\n"); return 1; }
      }
      np++;
    }
  }
  for (int N = 2; N < NMAX; N++) {
    if (dfirst[N]) printf("D %d %llu\n", N, (unsigned long long)dfirst[N]);
    if (ifirst[N]) printf("I %d %llu\n", N, (unsigned long long)ifirst[N]);
    if (rfirst[N]) printf("R %d %llu\n", N, (unsigned long long)rfirst[N]);
  }
  return 0;
}

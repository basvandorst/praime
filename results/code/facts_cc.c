// facts_cc.c: least start of a complete Cunningham chain of exactly k primes, both kinds, for FACTS.md section 6B.
// Build: gcc -O3 -o facts_cc facts_cc.c -lm. Run: ./facts_cc 30000000000 (about 6 minutes). Only p = 29 (mod 30) and p = 1 (mod 30) are tried (fact 6201),
// so chains starting at 2, 3 or 5 are not reported. Chain members are tested with deterministic Miller-Rabin (12 bases, exact below 3.18e23).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 mulm(u64 a, u64 b, u64 m) { return (u128)a * b % m; }
static u64 powm(u64 a, u64 e, u64 m) { u64 r = 1; a %= m; while (e) { if (e & 1) r = mulm(r, a, m); a = mulm(a, a, m); e >>= 1; } return r; }
static int isprime(u64 n) {                     // deterministic Miller-Rabin for n < 3.3e24 (bases <= 37 suffice below 2^64)
  static const int sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
  if (n < 2) return 0;
  for (int i = 0; i < 25; i++) { if (n == (u64)sp[i]) return 1; if (n % sp[i] == 0) return 0; }
  u64 d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
  for (int i = 0; i < 12; i++) { u64 x = powm(sp[i], d, n); if (x == 1 || x == n - 1) continue;
    int ok = 0; for (int r = 1; r < s; r++) { x = mulm(x, x, n); if (x == n - 1) { ok = 1; break; } } if (!ok) return 0; }
  return 1; }
static u64 first[2][16];
static void chain(u64 p, int kind) {             // p prime, p = 29 or 1 (mod 30)
  long d = kind ? -1 : 1; u64 q = p; int L = 1;
  while (1) { q = 2 * q + d; if (!isprime(q)) break; L++; }
  if (L < 4) return;
  u64 pre = kind ? (p + 1) / 2 : (p - 1) / 2; if (isprime(pre) && 2 * pre + d == p) return;
  if (L < 16 && !first[kind][L]) { first[kind][L] = p; printf("kind %d length %d start %llu\n", kind + 1, L, (unsigned long long)p); fflush(stdout); }
}
int main(int argc, char **argv) {
  u64 X = strtoull(argv[1], 0, 10), R = (u64)sqrt((double)X) + 2;
  unsigned char *bs = calloc(R + 1, 1); u64 *bp = malloc(sizeof(u64) * R); long nb = 0;
  for (u64 i = 2; i <= R; i++) if (!bs[i]) { if (i > 2) bp[nb++] = i; for (u64 j = i * i; j <= R; j += i) bs[j] = 1; }
  const u64 S = 1 << 22; unsigned char *seg = malloc(S);
  for (u64 base = 3; base <= X; base += 2 * S) {
    memset(seg, 0, S); u64 top = base + 2 * S;
    for (long t = 0; t < nb; t++) { u64 p = bp[t]; if (p * p >= top) break;
      u64 st = p * p; if (st < base) { st = ((base + p - 1) / p) * p; if (!(st & 1)) st += p; }
      for (u64 j = st; j < top; j += 2 * p) seg[(j - base) / 2] = 1; }
    for (u64 i = 0; i < S; i++) { u64 p = base + 2 * i; if (p > X) break;
      if (seg[i]) continue; u64 r = p % 30;
      if (r == 29) chain(p, 0); else if (r == 1) chain(p, 1); }
  }
  for (int k = 0; k < 2; k++) { printf("KIND %d:", k + 1); for (int L = 4; L < 16; L++) printf(" %d:%llu", L, (unsigned long long)first[k][L]); printf("\n"); }
}

/* Puzzle 795.  F(a, m): reverse the digit strings between successive occurrences of the digit a in m.
   Wanted: the smallest prime p that uses exactly D distinct digits, such that F(a, p) is a prime
   different from p for every digit a of p.  A digit used once is a fixed point, so every used digit
   occurs at least twice.  This program enumerates, in increasing numeric order, all numbers with
   exactly L digits that use exactly D distinct digits each at least twice (depth-first over the
   digits, most significant first), and tests p and its D transforms with a deterministic
   Miller-Rabin (64-bit, L <= 19).  Filters: last digit in {1,3,7,9} (shared by all transforms),
   digit sum not divisible by 3 (shared too).  Stops at the first solution.
   Usage: p795_search L D [count]      (count: report the first `count` solutions)
   gcc -O2 -o p795_search p795_search.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
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
static int L, D, want, found = 0; static int dig[24], cnt[10], used; static long long tested = 0;
static uint64_t value(const int *d) { uint64_t v = 0; for (int i = 0; i < L; i++) v = v * 10 + d[i]; return v; }
static int check(void) {
    uint64_t p = value(dig); tested++;
    static const uint32_t tp[] = {7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
    for (int i = 0; i < 22; i++) if (p % tp[i] == 0) return 0;   /* cheap filter before Miller-Rabin */
    if (!isprime64(p)) return 0;
    int t[24];
    for (int a = 0; a < 10; a++) if (cnt[a]) {
        memcpy(t, dig, sizeof(int) * L);
        int prev = -1;
        for (int i = 0; i < L; i++) if (dig[i] == a) {
            if (prev >= 0) { int x = prev + 1, y = i - 1; while (x < y) { int tmp = t[x]; t[x] = t[y]; t[y] = tmp; x++; y--; } }
            prev = i;
        }
        uint64_t f = value(t);
        if (f == p || !isprime64(f)) return 0;
    }
    printf("L=%d D=%d p=%llu transforms:", L, D, (unsigned long long)p);
    for (int a = 0; a < 10; a++) if (cnt[a]) {
        memcpy(t, dig, sizeof(int) * L); int prev = -1;
        for (int i = 0; i < L; i++) if (dig[i] == a) { if (prev >= 0) { int x = prev + 1, y = i - 1; while (x < y) { int tmp = t[x]; t[x] = t[y]; t[y] = tmp; x++; y--; } } prev = i; }
        printf(" F(%d)=%llu", a, (unsigned long long)value(t));
    }
    printf("  (candidates tested: %lld)\n", tested); fflush(stdout);
    return 1;
}
static void dfs(int pos, int need) {      /* need = positions still required to bring every used digit to 2 and to introduce the missing digits */
    if (found >= want) return;
    int left = L - pos;
    if (need > left) return;
    if (pos == L) {
        int s = 0; for (int a = 0; a < 10; a++) s += a * cnt[a];
        if (s % 3 == 0) return;
        if (check()) found++;
        return;
    }
    for (int a = (pos == 0 ? 1 : 0); a < 10; a++) {
        if (pos == L - 1 && !(a == 1 || a == 3 || a == 7 || a == 9)) continue;
        int newdig = (cnt[a] == 0);
        if (newdig && used == D) continue;
        int dneed = cnt[a] < 2 ? -1 : 0;                  /* first and second copies of a digit each supply one required position */
        cnt[a]++; used += newdig; dig[pos] = a;
        dfs(pos + 1, need + dneed);
        cnt[a]--; used -= newdig;
        if (found >= want) return;
    }
}
int main(int argc, char **argv) {
    L = atoi(argv[1]); D = atoi(argv[2]); want = argc > 3 ? atoi(argv[3]) : 1;
    if (L > 19 || 2 * D > L) { fprintf(stderr, "bad parameters\n"); return 1; }
    dfs(0, 2 * D);
    if (!found) printf("L=%d D=%d: no solution (candidates tested: %lld)\n", L, D, tested);
    return 0;
}

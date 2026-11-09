// Puzzle 823.  Strings over {1,3,7,9} of length L; K = number of distinct primes among the substrings read
// left to right.  Simulated annealing over single-digit changes (plus occasional swaps), full re-evaluation with
// a memo of substring primality (64-bit Miller-Rabin up to 19 digits, GMP mpz_probab_prime_p beyond; GMP is
// linked by hand-written prototypes since the header is not installed).  Seeds: random or a given string.
// Usage: p823_anneal L iterations seed [start_string]      g++ -O2 -o p823_anneal p823_anneal.cpp -l:libgmp.so.10
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <random>
extern "C" {
typedef struct { int _mp_alloc; int _mp_size; void *_mp_d; } __mpz_struct; typedef __mpz_struct mpz_t[1];
void __gmpz_init(mpz_t); void __gmpz_clear(mpz_t); int __gmpz_set_str(mpz_t, const char *, int); int __gmpz_probab_prime_p(const mpz_t, int);
}
typedef unsigned __int128 u128;
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) { return (u128)a * b % m; }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) { uint64_t r = 1; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static bool isprime64(uint64_t n) {
    if (n < 2) return false;
    static const uint64_t sp[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) { if (n == sp[i]) return true; if (n % sp[i] == 0) return false; }
    uint64_t d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        uint64_t x = powmod(sp[i], d, n); if (x == 1 || x == n - 1) continue;
        bool ok = false; for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { ok = true; break; } }
        if (!ok) return false;
    }
    return true;
}
static std::unordered_map<std::string, bool> memo;
static mpz_t z;
static bool isprime_str(const std::string &s) {
    auto it = memo.find(s); if (it != memo.end()) return it->second;
    bool r;
    if (s.size() <= 19) r = isprime64(strtoull(s.c_str(), 0, 10));
    else { __gmpz_set_str(z, s.c_str(), 10); r = __gmpz_probab_prime_p(z, 20) > 0; }
    if (memo.size() < 60000000) memo[s] = r;
    return r;
}
static int evaluate(const std::string &s) {
    std::unordered_set<std::string> seen; int L = s.size();
    for (int i = 0; i < L; i++) for (int j = i + 1; j <= L; j++) {
        std::string t = s.substr(i, j - i);
        if (isprime_str(t)) seen.insert(t);
    }
    return seen.size();
}
int main(int argc, char **argv) {
    int L = atoi(argv[1]); long iters = atol(argv[2]); unsigned seed = atoi(argv[3]);
    __gmpz_init(z);
    std::mt19937_64 rng(seed); const char D[4] = {'1', '3', '7', '9'};
    std::string s;
    if (argc > 4) s = argv[4]; else { s.resize(L); for (auto &c : s) c = D[rng() % 4]; }
    int cur = evaluate(s), best = cur; std::string bests = s;
    double T0 = 2.0, T1 = 0.2;
    for (long it = 0; it < iters; it++) {
        double T = T0 * pow(T1 / T0, (double)it / iters);
        std::string t = s;
        int mv = rng() % 10;
        if (mv < 8) { int k = rng() % L; char c; do c = D[rng() % 4]; while (c == t[k]); t[k] = c; }
        else if (mv == 8) { int a = rng() % L, b = rng() % L; std::swap(t[a], t[b]); }
        else { int k = rng() % L, m = rng() % L; t[k] = D[rng() % 4]; t[m] = D[rng() % 4]; }
        int val = evaluate(t);
        if (val >= cur || exp((val - cur) / T) > (double)(rng() % 1000000) / 1e6) {
            s = t; cur = val;
            if (cur > best) { best = cur; bests = s; printf("it %ld: K=%d %s\n", it, best, bests.c_str()); fflush(stdout); }
        }
    }
    printf("final best K=%d %s\n", best, bests.c_str());
    return 0;
}

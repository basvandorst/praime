# Results in short

PRP = probable prime. Full explanations and verification scripts: [results/](../results).

## Records beaten

- **65:** With p = 19765897565467, the 12 consecutive primes p+{0,16,22,54,60,76} and p+{4,6,36,40,70,72} have equal sums of k-th powers for k = 1 to 5.
- **87 Q1:** [53 251 197 / 311 167 23 / 137 83 281] + [257 11 269 / 191 179 167 / 89 347 101] + 1 = [311 263 467 / 503 347 191 / 227 431 383], total 1041.
- **87 Q2:** [41 113 59 / 89 71 53 / 83 29 101] + [29 167 71 / 131 89 47 / 107 11 149] + [67 79 163 / 199 103 7 / 43 127 139] = [137 359 293 / 419 263 107 / 233 167 389], total 789.
- **92 (n = 5):** bottom layers 23 | 61 127 | 41 29 71 | 139 43 7 103 | 17 101 137 89 37 (each ball above is the sum of the three below), top 4999, total 16561.
- **99:** p = 10^3000 + 10^1889 + 10^1500 + 10^1111 + 1 satisfies p² = R(R(p)²) (3001 digits, PRP).
- **134 Q3:** a 201-digit carrousel prime 10·B + {1,3,7,9}, B in [134.md](../results/134.md).
- **227:** t = a/2^2760 (⌊t⌋ = 3001, a has 835 digits): ⌊t^j⌋ is prime for j = 1 to 239. a in [data/227.txt](../results/data/227.txt).
- **249:** [c−b, c+a+b, c−a / c−a+b, c, c+a−b / c+a, c−a−b, c+b] with c = 320515678943, a = 3000000, b = 4554 (magic sum 961547036829).
- **307:** 7·(2^136279841 − 1)·10^2192 is a Smith number (41,026,513 digits).
- **324:** 99989897979696968585858484848473737373737373626262626262515140405151404051 (74 digits).
- **340:** p = ⌈10^3923·(1 + √22001)/2⌉ (3925 digits, PRP) is contained in p².
- **342:** opposite faces (4, 64), (36, 144), (49, 361), total 658 (minimal).
- **356:** max 661: 2 3 5 29 79 97 107 193 239 311 349 401 461 541 577 599 607 641 647 661. Max 1297: 2 3 13 31 103 151 281 307 347 389 433 449 577 599 673 761 797 983 1033 1063 1097 1117 1289 1291 1297. Max 2213: 2 3 7 23 89 139 193 211 467 503 547 631 659 719 877 997 1021 1301 1439 1451 1549 1613 1811 1873 1907 2053 2083 2161 2207 2213.
- **359:** the first 202 primes in a circle (202 is the maximum possible): [data/359.txt](../results/data/359.txt).
- **373:** 2^82589933 − 1 is a self number.
- **398 and 438:** N = A·10^900 + 1 (1057 digits, proven prime), A = 1057447767801982899810936078986870322855097083184870799512294142679453394807877910422539690997221686644824565309308629748242833059104863529225434120688108130.
- **419:** N = product of 336 primes with SOPF = 3^42·5^25 gives σ(N) = (SOPF·τ(N))^47. Primes in [data/419.txt](../results/data/419.txt).
- **450:** the sum of the twin primes 2996863034895·2^1290000 ± 1 has 1,290,006 prime factors.
- **476:** a 458-digit seed with 30 consecutive prime odd terms (PRP): [data/476.txt](../results/data/476.txt).
- **522:** x1 = 3355845769, y1 = 5130974053 give k = 15.
- **553:** p = (523·541·…·3001)/2593 − 1296, q = p + 2592 (1060 digits): p + … + q = 523·541·…·3001 (333 primes).
- **597:** x² − 9999x + 25242757 gives 9743.
- **607:** (10^1888529 − 10^944264 − 1)^7 has 13,219,703 digits and no digit 6.
- **793:** K = 6: 1012031234062461463. K = 8: 10203123405142534560716273657467. K = 10: 101203123405142534560716273647567809182938495869789 (PRP).
- **808:** 108864424468801 is palindromic in 227 bases.
- **830:** K = 2012: N = 3520·21·4·2^11·3^10 = 35757277839360. K = 2013: N = 60·54·2^13·3^8 = 174142586880.
- **905:** P(x) = x + d + 18·(x − a0)(x − a1)…(x − a25), with PrimeGrid's AP-27 (a0 = 224584605939537911, d = 81292139·223092870), gives 28 primes by iteration (last one PRP).
- **908:** p = 68516818704324257201795329467591133865835633315118443095966860117, q = p² − p + 199 give 34 remainders 199, 197, 193, …, 5, 2, 1 (PRP).
- **957:** 1-2+3\*(4\*5^6)^(7\*89) (2989 digits), -(1+2)+34^(5\*6\*7+890) (1685), -9+8\*7^(6\*54\*3^2)\*1 (2466), -9+8\*76^(54-3+2^10) (2023), all PRP.
- **1001:** P = 391177608521960466267478281835624430526821956454882813 gives Q/P = 5.5201061257.
- **1028:** (35968n mod 35969) + (9481023n mod 94764426) + (204182n mod 4083624) + (1230n mod 36886) + (13256n mod 530202) gives 49 primes for n = 1 to 49.
- **1097:** 9699690·(n mod 10) + Σ e_j·((⌊n/10⌋ + j) mod 13), e = [933516151, −104434997, 23700127, −128911085, 56766199, 48442819, −1599797, −152799521, 32343931, −71698241, 45728731, 44720503, 43956415]/1014, gives 130 distinct primes.
- **1099:** H(n) = 35969 + 9699690n + Σ c_j·⌊n/(10j)⌋ gives 570 primes for n = 0 to 569. Coefficients in [data/1099.txt](../results/data/1099.txt).
- **1135:** 5n^4 − 330n^3 + 9395n^2 − 130510n + 709337 gives 30 primes (n = 0 to 29) whose digits never alternate parity.
- **1190:** 9999999967 = ((9·9+9)/9)^(9+9/9) − 6·7 + 9.
- **1209:** −3x² + 2799x + 22433 gives 34 unlucky primes (x = 0 to 33). The listed 45 contains five lucky numbers.
- **1217:** with h = (p(n+1) − p(n))/2, h² − 79h + 1601 is prime for 1,319,943 successive n; |S(h − 43)| with the Dress–Landreau sextic S gives 103,518 n (53 distinct primes).
- **1230:** 5n^4 − 330n^3 + 9395n^2 − 130510n + 709337 gives 30 primes (n = 0 to 29), all with an even-length period of 1/p.
- **1251:** b = 518918399 gives 477 cycles below 10^9: [data/1251.txt](../results/data/1251.txt).

## First answers to open questions

- **75 Q2:** the smallest palprime with digit sum 2000 has 223 digits: 10^223 − 1 − (10^221 + 10^220 + 10^124 + 10^111 + 10^98 + 10^2 + 10) (PRP).
- **77:** a 10-row prime pine, bottom row 3548161 12418561 54405121 351859201 368417281 3437575681 8200980481 18681062401 47753502721 809036659201 (each row above is the average of the two below: 55 distinct primes).
- **81 c):** squares [c+x, c−x−y, c+y / c−x+y, c, c+x−y / c−y, c+x+y, c−x]. For p, 2p+1, 4p+3: c = 48353471, x = 38813580, y = 8220240 (magic sum 145060413). For p, 2p−1, 4p−3: c = 47722789, x = 39390450, y = 2801400 (magic sum 143168367).
- **87 Q3:** 4×4 versions with constants 548 and 570, squares in [87.md](../results/87.md).
- **92 (n = 6, 7):** bottom layers n = 6: 43 | 137 19 | 61 59 79 | 173 37 197 211 | 223 53 13 191 17 | 5 229 71 67 73 283. n = 7: 409 | 457 131 | 61 83 139 | 47 281 313 71 | 241 331 79 41 379 | 7 59 11 431 181 593 | 13 557 127 211 19 347 109.
- **106:** for 2000, k = 4: [17 89 71 / 113 59 5 / 47 29 101] + [101 29 47 / 5 59 113 / 71 89 17] + [269 1493 1061 / 1733 941 149 / 821 389 1613] + [1613 389 821 / 149 941 1733 / 1061 1493 269]; k = 2 is impossible, k = 26 is the maximum. For 2001: at least 3, at most 27 ([106.md](../results/106.md)).
- **134 Q4:** 9773317917193993331131917339173977131919931391113371131103919313311737119331717179393197317197911119.
- **170 Q10:** impossible. Since n+1 is even, the pseudoprime condition gives a nontrivial divisor of n.
- **172:** K = 9, L = 113: 1.1.1.1.1.1.1.3.19.17.7.23.2.13.5.11.1.1.1.1.1.1. K = 10, L = 144: 1.1.1.1.1.1.23.1.11.17.3.13.19.29.5.7.1.1.1.2.1.1.1.1.1. Both minimal (exhaustive search).
- **179 Q1:** TROIS + ONZE = SEPT + SEPT gives 12689 + 6053 = 9371 + 9371, all prime. Also CINQ + CINQ + ONZE = 3·SEPT gives 2879 + 2879 + 6701 = 3·4153.
- **189:** paths exist for 52, 53, 56, 58 and every open n above; none for 44, 45, 46, 47, 49, 51, 54 and 57, so n0 = 57. Paths in [data/189.txt](../results/data/189.txt).
- **299 Q2:** 35014031359.
- **325 Q4:** p1 = 12487, p2 = (10^21 + 1245)·7723# + 1, A = (p2 − 1)/12486 − 1, B = 12487 − A, p3 = A·p2 + B; Z = p1·p2·p3 has 10,008 digits (p3 PRP).
- **435:** m = 60493590969525342720 = 2^9·3^5·5·7^3·13^2·31^2·41·42577 satisfies φ(m) + σ(m) = 5m.
- **458:** n = 678395112906914180368298499070649370.
- **490:** no other odd n exists, and k = 3 always (Ramanujan–Nagell for k = 3, Beukers for other k).
- **493:** yes, for every odd prime p and every x + y = p; a and b come from expanding (√x + i√y)^p.
- **516 Q1:** 59 is the only solution. No solution has more than 7 digits, and a check below 10^7 finds only 59.
- **631:** yes, by Shiu's theorem on strings of congruent primes.
- **639:** n = P·(P + 6m + 2) whenever P, P+2, P+6m, P+6m+2 are all prime. Checked for every m up to 3.3·10^7; holds for all m under the prime k-tuples conjecture.
- **677:** the tried values are impossible (neighbours must alternate between 1 and 3 mod 4); the possible range is narrowed.
- **890 Q1:** impossible (proof via parity and mod 4).
- **1080:** σ3(54608832726) = σ3(56689652415) = σ3(57381413815) = 190592992522906181161863557566464, exactly three solutions, so a(3) ≤ 1.9·10^32.
- **1145:** the covering set behind Wesolowski's k is {5, 7, 11, 31, 73, 97, 151, 241, 631, 673, 1321, 23311}. His a = 153479820268467961 is minimal for it, and no cheaper set has period ≤ 1440.
- **1204:** Hardy–Littlewood predicts 51.7%, because 30030 = 2·3·5·7·11·13.

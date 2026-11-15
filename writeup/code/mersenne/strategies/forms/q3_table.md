| P | h(-4P) | criterion (method A, qf_extra.py) | primes checked | agree | exceptions | method B (qf_extra_check.py): class polynomial H_{-4P} | method B: simple polynomial |
|---:|---:|---|---:|---:|---|---|---|
| 27 | 3 | p = 1 (mod 3) and 2^((p-1)/3) = 1 (mod p) | 41,538 | 41,538 | 0 | 41,525 of 41,525 agree (skipped p dividing 4P*disc H: [2, 3, 5, 11, 17, 23, 29, 41, 47, 59, 71, 89, 101]) | p = 1 mod 3 and x^3 - 2 has a root: 41,538 of 41,538 agree; same field as polredabs(H) = x^3 - 2: True |
| 64 | 4 | p = 1 (mod 4) and 2^((p-1)/4) = 1 (mod p) | 41,538 | 41,538 | 0 | 41,517 of 41,517 agree (skipped p dividing 4P*disc H: [2, 3, 7, 11, 19, 23, 31, 43, 47, 59, 71, 83, 103, 107, 131, 139, 191, 211, 227, 239, 251]) | p = 1 mod 4 and x^4 - 2 has a root: 41,538 of 41,538 agree; same field as polredabs(H) = x^4 - 2: True |
| 14 | 4 | (-14/p) = 1 and (x^2+1)^2 = 8 (mod p) solvable | 41,536 | 41,536 | 0 (set aside: p=2: drawn=no, criterion=no; p=7: drawn=no, criterion=no) | 41,527 of 41,527 agree (skipped p dividing 4P*disc H: [2, 7, 11, 17, 29, 31, 37, 41, 43, 47, 53]) | (-14/p) = 1 and (x^2+1)^2 - 8 has a root: 41,538 of 41,538 agree; same field as polredabs(H) = x^4 - 2*x^3 + x^2 - 2*x + 1: True |
| 11 | 3 | x^3 - x^2 + x + 1 has three roots mod p | 41,536 | 41,536 | 0 (set aside: p=2: drawn=no, criterion=no; p=11: drawn=yes, criterion=no) | 41,530 of 41,530 agree (skipped p dividing 4P*disc H: [2, 7, 11, 13, 17, 19, 29, 41]) | x^3 - x^2 + x + 1 has three roots: 41,537 of 41,538 agree, exceptions [11]; same field as polredabs(H) = x^3 - x^2 + x + 1: True |
| 23 | 3 | x^3 - x - 1 has three roots mod p | 41,536 | 41,536 | 0 (set aside: p=2: drawn=no, criterion=no; p=23: drawn=yes, criterion=no) | 41,526 of 41,526 agree (skipped p dividing 4P*disc H: [2, 5, 7, 11, 17, 19, 23, 53, 61, 79, 83, 89]) | x^3 - x - 1 has three roots: 41,537 of 41,538 agree, exceptions [23]; same field as polredabs(H) = x^3 - x^2 + 1: True |

Further consistency rows of qf_extra.py:

| statement | primes checked | agree | disagree |
|---|---:|---:|---:|
| P=11 (-11/p)=1 and x^3-x^2+x+1 has a root | 41,536 | 41,536 | 0 |
| P=27 x^3-2 three roots | 41,536 | 41,536 | 0 |
| P=64 (p=1 mod 4) and x^4-2 has a root | 41,537 | 41,537 | 0 |
| P=14 x^4+2x^2-7 has a root and (-14/p)=1 | 41,536 | 41,536 | 0 |
| x^3-x-1 root-count distribution (p != 23) | {'0': 13809, '1': 20818, '3': 6910} | | |

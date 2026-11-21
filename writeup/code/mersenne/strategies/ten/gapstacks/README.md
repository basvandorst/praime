# Prime gap stacks (prime40.html)

Falling blocks with vertical bars only: the gap from p_n to p_(n+1) is a bar of max(1, gap/2) cells, dropped in order into a box of W columns and H rows; a full row clears, and the game ends when the next bar does not fit. The column is the lowest one (leftmost on ties), a random one (fixed LCG, seed 20261121) or p_n mod W. The page sieves the primes up to 20,000,000 itself. (An earlier version of the page dropped round balls; it was replaced because frozen balls let small balls overlap.)

- `ref40game.py` → `ref40game.json`: an independent replay of the rules in Python (numpy sieve) for five settings.
- `probe2.js`: plays the same five games on the page at full speed; `compare40.txt`: page against replay, 0 mismatches. For example 12 × 40 with the lowest column: 14,356 bars and 6,492 cleared rows, ended by the gap 86 after 155,921 (43 cells, taller than the box); 12 × 40 with p mod 12: 47 bars, 0 rows, because columns 0, 4, 6, 8, 9, 10 never get a bar (a column sharing a factor with W gets at most the bar of a prime dividing W).
- `ref40.py` → `ref40.json`: the first gaps and their counts (sympy).
- `t40.js`: the headless runner.

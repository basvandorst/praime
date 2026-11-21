# Prime gap balls (prime40.html)

One ball per gap between consecutive primes (diameter = the gap, in order 1, 2, 2, 4, 2, 4, …) falls into a box 100 units high; the page lets you set the width (10 to 300 units) and the speed. The page sieves the primes up to 2,000,000 itself.

- `ref40.py` → `ref40.json`: the first 200 gaps and the gap counts, recomputed with sympy; the page's gaps agree (first 30 checked in headless Chrome, and the legend after a full box of width 60: gaps 1, 2, 4, 6, 8, 10, 12, 14, 18 occur 1, 30, 32, 39, 12, 13, 6, 4, 1 times among the first 138, the primes up to 797).
- `probe.js`: fast-forwards the simulation until the box is full and reports balls, steps, the largest overlap (relative to the smaller radius) and balls outside the box. Widths 60, 200, 300: full after about 141, 315, 406 balls, 0 outside, overlap at most 0.47 (four solver passes per step). Width 10: stopped by the gap 14 after 113.
- `t40.js`: the headless runner (`node t40.js $PWD/prime40.html '[{"q":"","tag":"x","wait":6000}]'`).

The physics (Verlet, gravity, walls, pairwise pushes weighted by r², a little damping) makes a picture, not a packing measurement.

# Problem 88, Q3: a magic square with 7 squares and two triangular numbers? (attempt)

> **In short:** Scaling the square by K² makes both non-squares 222121 and 360721 triangular exactly when 8·222121·K² + 1 and 8·360721·K² + 1 are both squares. This is a pair of simultaneous Pell equations. The page's "no" (Jan van Delden) subtracts the two equations as if they shared the same x, which they need not. Every solution has K = k_n(a) = k_m(b), a common term of the two Pell sequences. That forces |n·ln ε_a − m·ln ε_b − ½·ln(a/b)| < 2·10⁻²⁷. A scan of every n ≤ 10⁶ finds the smallest value of that linear form to be 4.8·10⁻⁵. So there is **no K > 0 below about 10^(1.37·10⁷)**. A proof for all K needs a lower bound for linear forms in two logarithms (Laurent) plus a continued-fraction reduction, which we did not carry out.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/problems/prob_088.htm)

**Status:** open · no solution with K < 10^(1.37·10⁷) (proof); all K needs a Baker-type bound

## The question

Multiply every entry of Bremner's magic square (7 squares plus 222121 and 360721) by K². The squares stay squares. Q3 asks whether some K makes both remaining entries triangular.

## Reduction

s·K² = t(t+1)/2 exactly when (2t + 1)² − 8s·K² = 1. So we need, with a = 222121 and b = 360721,

x² − 8a·K² = 1  and  y² − 8b·K² = 1  with the same K > 0.

van Delden's argument subtracts the two equations as if x = y, but in general x ≠ y. Eliminating K gives b·x² − a·y² = b − a. This always has the solution x = y = 1 (K = 0), so congruences cannot rule out the problem.

## Linear forms in logarithms

The fundamental solutions are (h, k) = (22632419659503, 16978181966) for 8a, and (86286108053753141137230600854902563, 50793737929709259543442603308701) for 8b. We checked both with `diop_DN`.

All solutions are k_n = (εⁿ − ε⁻ⁿ)/(2√(8s)), with ε = h + k√(8s). So K = k_n(a) = k_m(b) with n, m ≥ 1. Taking logarithms gives

|n·ln ε_a − m·ln ε_b − ½·ln(a/b)| ≤ 2(ε_a⁻² + ε_b⁻²) < 2·10⁻²⁷.

For each n ≤ 10⁶, take the nearest m and evaluate the left side at 80-digit precision. Its minimum is 4.8·10⁻⁵, at (n, m) = (321961, 124773). So no solution has n ≤ 10⁶, which means none has K < k_(10⁶+1)(a) ≈ 10^(1.37·10⁷).

## What would finish it

A lower bound for |n·ln ε_a − m·ln ε_b − c|, from Laurent's theorem on linear forms in two logarithms, caps n at some explicit value, probably around 10²⁵. A Baker–Davenport reduction with the continued fraction of ln ε_a / ln ε_b then lowers the cap below 10⁶, and the scan above covers the rest. Results of Bennett and of Bennett, Cipu, Mignotte and Okazaki also limit such systems to at most two solutions.

## Verify it yourself

`writeup/code/prob88_pell.py` checks the fundamental solutions and scans n ≤ 10⁶. It takes about 20 minutes.

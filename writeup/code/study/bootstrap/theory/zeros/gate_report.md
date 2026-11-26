# Prior-art search for the zero traps (26 November 2026)

Report of a literature agent, as handed back (condensed; wording of verdicts unchanged). Not checked line by line by a second reader.

## The statements under test

(A) Trap lemma: a feasible function gt of cost kt and a function h of the same finite-dimensional space with eta*hhat <= gthat and eta*W(h) > kt force, at every higher degree, a zero of the optimal function in the set {eta*h > gt}.
(B) The dual optimal measure is a Gauss-Radau type quadrature formula (K free nodes, one fixed node at 0) for the prime-power measure.
(C) The double zeros of the best polynomial converge to the prime powers.

## Verdicts

**(A) NOT FOUND, provisional.** Karlin-Studden, Krein-Nudelman and Krein's 1951 survey were not read (no text access), so the verdict is void until they are checked. Closest:
1. Cohn-Kumar 2009 (arXiv math/0403263), Prop 4.3 and Lemma 4.4: weak duality localises the vectors of a lattice (a dual point of the full problem) near the roots of f; Lemma 4.4 uses a Hermite-type function (single root at one norm, double roots at the others) to force at least 196560 minimal vectors. Differences: full-problem dual objects, not truncated dual optima; no constraint on the transform side.
2. Johns, Part I, Prop 5.5 and Cor 5.3: window-mass bounds for every admissible pair (upper bounds; nothing about finite-degree optima).
3. Chebyshev-Markov-Stieltjes separation (quoted by Denisov-Simon, arXiv math/0209329, from Freud): statements about Gauss nodes, not a trap built from a feasible pair. The generic step "a positive rule exact on V has a node where some f in V is positive" is folklore.

**(B) CLOSEST PRIOR ART.** The principle is classical and in print for analogous linear programmes; no instance for the span of Xi^2 times polynomials was found.
- Cohn-Goncalves 2019 (arXiv 1712.04438), section 5, Prop 5.1 and Cor 5.2: a stationary optimum with free double roots gives a summation formula with k+2 nodes, including the fixed node 0, on the (2k+2)-dimensional space ("reminiscent of Gauss-Jacobi quadrature").
- Afkhami-Jeddi, Cohn, Hartman, de Laat, Tajdini 2020 (arXiv 2006.02560), eq. (2.23).
- Markov-Krein maximal-mass theorem and the principal representation through an endpoint (Karlin-Studden Thm 4.1, Thm 3.1), seen only through Owhadi-Scovel-Sullivan (arXiv 1304.7046). Caveat: the classical theorem needs a T-system on [0, oo) and an inactive constraint on the transform side.

**(C) NOT FOUND for the explicit formula; no proof of root convergence was found in any field.** The analogues are numerical or conjectural: Cohn-Goncalves Table 4.2 and Conj 4.2; Hartman-Mazac-Rastelli section 2.5 ("zeros ... appear to converge toward the non-negative integers"); Afkhami-Jeddi et al. Conj 4.1 ("no proof"); Cohn-Miller section 3 ("a proof is completely lacking"); Cohn-Elkies Tables 1-2. Warning from two of these papers: at fixed truncation the top third of the roots does not track the target set, so (C) must be stated per fixed index as the degree grows.

## The Johns preprints

- Versions found: Part I v1.4 and Part II v1.3, both repositories last pushed 2026-10-08; no third part on GitHub. The Zenodo pages could not be opened (403, then a timeout).
- A companion draft "Hermite interpolation at the prime powers" lies in the repository of Part I (supplementary folder). It works with fixed nodes only. A remark there says that the near-extremal certificates do not vanish at 6, 10, 12.
- Not in the corpus as read: optimisation with free nodes; a trap lemma; any convergence claim. Part II, section 1.5: "we claim nothing about the convergence of the exact Hermite family of Part I or of other interpolation schemes"; section 9.4(3) lists it as open.
- Part I, Thm 2.12(a, c): weak limits of minimising sequences that are magic vanish on the support of every admissible measure (so at every prime power under RH). This is about zeros of the limit, not about the nodes at finite degree.
- Part II's magic function vanishes doubly at all integers n >= 2, so a limit of the free-node optima need not be that function.

## Read in full text

Johns Part I v1.4 (sections 1.3-1.6, 2.4-2.5, Thm 3.6, Cor 3.9, all of 4, 5.2-5.4, all of 6, App. B.2.3); Part II v1.3 (1.5, 1.6, 8(6), all of 9); the companion draft (parts); Cohn-Elkies sections 7-8; Cohn-Kumar section 4; Cohn-Miller sections 1 and 3; Cohn-Goncalves sections 4-5; Afkhami-Jeddi et al. 2.4 and 4.1; Hartman-Mazac-Rastelli 2.5; El-Showk-Paulos section 3; Mazac 2016 3.1 and 8; Denisov-Simon pp. 1-2; Owhadi-Scovel-Sullivan sections 4-5.

Keyword sweeps only: Cohn-Triantafillou; de Laat-Leijenhorst; Carneiro-Milinovich-Soundararajan; Chirre-Goncalves-de Laat; Odlyzko 1990; Poitou 1977 and the Bourbaki talk; Mestre; Booker; Miller. None states that optimal explicit-formula test functions vanish at or near log n.

Named but not read: Karlin-Studden; Krein-Nudelman; Krein 1951; Freud; Szego; Afkhami-Jeddi-Hartman-Tajdini 2019; El-Showk-Paulos 2016; the Zenodo pages.

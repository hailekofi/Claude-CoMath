# Verified bibliography for the skeleton/core rebuild

Literature pass 2026-06-10 (web-verified against journals/arXiv/DOI). **All 22 items
EXIST.** Use only these keys in the .tex. Corrections vs. my tentative attributions are
flagged ⚠. Per-quantity mapping is in `ws_memo_nonrigidity_transcendence.md` (Part II).

## Physics: Type-1 / MLZ
- **[BE1993]** S. Brundobler, V. Elser, "S-matrix for generalized Landau–Zener problem,"
  *J. Phys. A* **26** (1993) 1211–1227. DOI 10.1088/0305-4470/26/5/037. — the BE law.
- **[SC2017]** N. A. Sinitsyn, V. Y. Chernyak, "The quest for solvable multistate
  Landau–Zener models," *J. Phys. A* **50** (2017) 255203. arXiv:1701.01870. — canonical
  modern review of solvable MLZ + integrability conditions. (Bow-tie origin: Demkov–
  Ostrovsky 2001 / Ostrovsky–Nakamura 1997, cite if the bow-tie itself is invoked.)
- **[OWY2009]** H. K. Owusu, K. Wagh, E. A. Yuzbashyan, "The link between integrability,
  level crossings and exact solution in quantum models," *J. Phys. A* **42** (2009)
  035206. arXiv:0807.0259. ⚠ **2009, not 2010/11**; the integrability↔level-crossing link
  (the orientation/selector ground), **not** the Type classification.
- **[OY2011]** H. K. Owusu, E. A. Yuzbashyan, "Classification of parameter-dependent
  quantum integrable models...," *J. Phys. A* **44** (2011) 395302. arXiv:1106.1831. ⚠
  **this** is the Type-M definition/classification paper (N−M commuting operators linear
  in coupling; Types 1,2,3).
- **[Y2018]** E. A. Yuzbashyan, "Integrable time-dependent Hamiltonians, solvable
  Landau–Zener models and Gaudin magnets," *Ann. Phys.* **392** (2018) 323–339.
  arXiv:1804.04926. — the cleanest cite for the **Type-1 ↔ Gaudin-magnet** map.
- **[VO2004]** M. V. Volkov, V. N. Ostrovsky, "Exact results for survival probability in
  the multistate Landau–Zener model," *J. Phys. B* **37** (2004) 4069–4084 (corr. **38**
  (2005) 907). — first analytic proof of BE (extremal state).
- **[DS2006]** B. E. Dobrescu, N. A. Sinitsyn, Comment on the above, *J. Phys. B* **39**
  (2006) 1253–1255. ⚠ a Comment correcting VO2004's change-of-variables; cite as the pair
  [VO2004, DS2006].

## Irregular connections / Stokes / wild geometry
- **[Wasow1965]** W. Wasow, *Asymptotic Expansions for ODEs*, Wiley 1965 (Dover repr.). —
  formal/Thome solutions, Stokes, turning points.
- **[JMU1981]** M. Jimbo, T. Miwa, K. Ueno, "Monodromy preserving deformation of linear
  ODEs with rational coefficients. I," *Physica D* **2** (1981) 306–352. — isomonodromy +
  τ-function, incl. irregular points.
- **[LR2016]** M. Loday-Richaud, *Divergent Series, Summability and Resurgence II*, LNM
  **2154**, Springer 2016. — summability / Stokes / Newton polygon.
- **[Boalch2001]** P. Boalch, "Symplectic manifolds and isomonodromic deformations,"
  *Adv. Math.* **163** (2001) 137–205. ⚠ **the** cite for "wild isomonodromy on these
  moduli = Painlevé/Garnier."
- **[Boalch2007]** P. Boalch, "Quasi-Hamiltonian geometry of meromorphic connections,"
  *Duke Math. J.* **139** (2007) 369–405. — symplectic/quasi-Ham structure of Stokes data.
- **[Boalch2014]** P. Boalch, "Geometry and braiding of Stokes data; fission and wild
  character varieties," *Ann. of Math.* **179** (2014) 301–365. arXiv:1111.6228. — wild
  character varieties (the dim-6 count).
- **[Boalch2010]** P. Boalch, "Towards a nonlinear Schwarz's list," in *The Many Facets of
  Geometry* (OUP 2010) 210–236. arXiv:0707.3375. — Painlevé/Schwarz-list survey (optional).

## Rigidity / middle convolution
- **[Katz1996]** N. M. Katz, *Rigid Local Systems*, Ann. Math. Studies **139**, Princeton
  1996. — rigidity index, Katz algorithm.
- **[DR2000]** M. Dettweiler, S. Reiter, "An algorithm of Katz and its application to the
  inverse Galois problem," *J. Symbolic Comput.* **30** (2000) 761–798.
- **[DR2007]** M. Dettweiler, S. Reiter, "Middle convolution of Fuchsian systems and the
  construction of rigid differential systems," *J. Algebra* **318** (2007) 1–24.
  — together: rigid ⇒ explicit (Euler/Laplace) ⇒ Γ-ratio connection data.

## Differential Galois / Liouvillian / transcendence
- **[Kolchin1973]** E. R. Kolchin, *Differential Algebra and Algebraic Groups*, Academic
  Press 1973. — Liouvillian ⇔ solvable differential Galois group.
- **[vdPS2003]** M. van der Put, M. F. Singer, *Galois Theory of Linear Differential
  Equations*, Grundlehren **328**, Springer 2003. — local theory; **Ramis density Thm 8.10**.
- **[Ramis1985]** J.-P. Ramis, "Phénomène de Stokes et filtration Gevrey sur le groupe de
  Picard–Vessiot," *C. R. Acad. Sci. Paris I* **301** (1985) 165–167 (full proof:
  Martinet–Ramis, *Ann. IHP* **54** (1991) 331–401). — cite via [vdPS2003] for a clean
  theorem statement.
- **[CS2007]** P. J. Cassidy, M. F. Singer, "Galois theory of parameterized differential
  equations and linear differential algebraic groups," IRMA Lect. **9**, EMS 2007.
  arXiv:math/0502396. ⚠ **this** (not Hardouin–Singer 2008, which is the *difference* case)
  is parameterized PV theory over a **differential** parameter field — the right frame for
  "σ non-Liouvillian as a function of (ε,γ,a)."

## Isomonodromy irreducibility (the transcendence engine)
- **[Nishioka1988]** K. Nishioka (Keiji), "A note on the transcendency of Painlevé's first
  transcendent," *Nagoya Math. J.* **109** (1988) 63–67.
- **[Umemura1990]** H. Umemura, "Second proof of the irreducibility of the first
  differential equation of Painlevé," *Nagoya Math. J.* **117** (1990) 125–171 (machinery:
  "Birational automorphism groups and differential equations," *Nagoya* **119** (1990) 1–80).
- **[Casale2008]** G. Casale, "Le groupoïde de Galois de P₁ et son irréductibilité,"
  *Comment. Math. Helv.* **83** (2008) 471–519. arXiv:math/0510657.
- **[CL2009]** S. Cantat, F. Loray, "Dynamics on character varieties and Malgrange
  irreducibility of Painlevé VI equation," *Ann. Inst. Fourier* **59** (2009) 2927–2978.
  — irreducibility of PVI **via dynamics on the character variety** = directly our frame.
- **[Ronveaux1995]** A. Ronveaux (ed.), *Heun's Differential Equations*, OUP 1995. — Heun /
  confluent Heun (the de-confluence target).

## Notes for the .tex
- Drop the merged/duplicated keys; the existing proof doc's `\bibitem`s for Wasow, JMU,
  Katz, Boalch, OWY, BE, Sinitsyn map to [Wasow1965],[JMU1981],[Katz1996],[Boalch2014],
  [OWY2009] (⚠ retag: this is the link paper; add [OY2011] for classification),[BE1993],
  [SC2017].
- For "rigid ⇒ Γ-ratio": [Katz1996]+[DR2000,DR2007]. For "non-rigid ⇒ Painlevé/Garnier":
  [JMU1981]+[Boalch2001]. For "Painlevé non-Liouvillian": [Nishioka1988],[Umemura1990],
  [Casale2008],[CL2009]. For "Liouvillian = solvable Galois / Ramis density":
  [Kolchin1973],[vdPS2003]. For "σ non-Liouvillian over the parameter field": [CS2007].

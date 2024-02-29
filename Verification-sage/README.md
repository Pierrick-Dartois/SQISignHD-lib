# Theta-dim4 sagemath library

This library is a `python/sagemath` implementation of the dimension 4 isogeny computations with the Theta model. Only chains of 2-isogenies and Theta models of level 2 are implemented. The library is suited to the computation of isogenies between products of elliptic curves given by Kani's lemma, especially in the context of SQISignHD verification. We recover isogenies in dimension 1 given their evaluation on torsion points. Hence, the library could also be used to test SIDH attacks.

Note that this code is experimental. It has not been completely tested and optimised.

## Requirement

You should have `python` and `sagemath` installed on your computer. Instructions to download `sagemath` can be found here https://doc.sagemath.org/html/en/installation/index.html.

## How to run the code?

Launch sagemath in the terminal (type command `sage`) and type 
`load("path to 'file' from your current directory")`
where `file` is the file containing the code you want to run.

## Generic tests

The script in `Tests.py` computes several instanciations of dimension 4 2-isogeny chains derived from Kani's lemma given as endomorphisms:
$$F:=\begin{array}{cccc} a_1 & a_2 & \widehat{\sigma} & 0 \\ -a_2 & a_1 & 0 & \widehat{\sigma} \\ -\sigma & 0 & a_1 & -a_2 \\ 0 & -\sigma & a_2 & a_1 \end{array}\in End(E_1^2\times E_2^2),$$
where:
- $E_1$ is a "random" supersingular elliptic curve, generated as the codomain of a random isogeny walk of degree $\simeq p$ starting from the supersingular elliptic curve $E_0$ of $j$-invariat $1728$;
- $\sigma: E_1\longrightarrow E_2$ is an isogeny of degree $\ell_B^{e_B}$ with $\ell_B=3$ or $7$; 
- Integers $a_1, a_2, e_A$ have been chosen such that:
$$a_1^2+a_2^2+\ell_B^{e_B}=2^{e_A};$$
- The characteristic $p$ is of the form $p:=c\cdot 2^{e_A+2}\cdot 3^{e_B}-1$, so that we have enough accessible torsion defined over $\mathbb{F}_{p^2}$ to compute $\sigma$ and evaluate it on $E_1[2^{e_A+2}]$. 

Different set of parameters $a_1, a_2, e_A, e_B, p$... can be found in `parameters/parameters.txt` (for $\ell_B=3$) and `parameters/parameters_7.txt` (for $\ell_B=7$). They all have been generated with the function `find_prime_gen` of `parameters/parameter_generation.py`.

This script includes two different ways of computing such a dimension 4 isogeny:
- the function `test_kani_endomorphism` tests instances when we have the full torsion available.
- the function `test_kani_endomorphism_half` tests instances when only half the torsion is available (as in SQISignHD, see Section 4.3 of the paper).

### Chain with the full torsion

When we can access the full $2^{e_A+2}$-torsion on $E_1$ and $E_2$, the whole isogeny chain $F: E_1^2\times E_2^2\longrightarrow E_1^2\times E_2^2$ can be computed directly. The function `test_kani_endomorphism` tests this direct chain computation for specified sets of parameters. 

INPUT:
- index: specifies the index of the set of parameters in the list extracted from the file `parameters/parameters.txt`
	if $\ell_B=3$ or `parameters/parameters_7.txt` if $\ell_B=7$.
- $\ell_B$: prime specifying the degree of the embedded isogeny $\sigma$ ($\deg(\sigma)=\ell_B^{e_B}$), $\ell_B=3$ or $7$ ($7$ by default).

OUTPUT: an object of the class `KaniEndo` (imported from `isogenies/Kani_endomorphism.py`) representing a dimension 4 $2$-isogeny chain derived from Kani's lemma.

EXAMPLES:

```sage: load("Tests.py")
sage: F=test_kani_endomorphism(1)
Testing KaniEndo with parameters:
 - Prime characteristic p = 1 * 2**21 * 7**5 - 1
 - Degree of the embedded isogeny sigma q = 7**5
 - a1 = 37
 - a2 = -336
 - Length of the dimension 4 2-isogeny = 17
Parameter import: 0.18054890632629395 s
Random walk: 0.0447688102722168 s
Generation of sigma: 0.005819082260131836 s
Generation and evaluation of the torsion basis: 0.007441997528076172 s
Dimension 4 endomorphism: 0.20083093643188477 s
Is evaluation correct?
True
Time evaluation: 0.008864164352416992 s
```

```sage: F=test_kani_endomorphism(1,3)
Testing KaniEndo with parameters:
 - Prime characteristic p = 1 * 2**19 * 3**9 - 1
 - Degree of the embedded isogeny sigma q = 3**9
 - a1 = 213
 - a2 = 22
 - Length of the dimension 4 2-isogeny = 16
Parameter import: 0.009882926940917969 s
Random walk: 0.018561124801635742 s
Generation of sigma: 0.010208845138549805 s
Generation and evaluation of the torsion basis: 0.004849910736083984 s
Dimension 4 endomorphism: 0.16254210472106934 s
Is evaluation correct?
True
Time evaluation: 0.007250070571899414 s
```

### Chain with half the necessary torsion

### Implementation restrictions


## Verifying real SQISignHD signatures



## License

SQIsign is licensed under Apache-2.0. See LICENSE and NOTICE in the root directory.

Third party code is used in this directory (Verification-sage):

- `Verification/montgomery_isogenies`; MIT: "Copyright (c) 2023 Giacomo Pope"
- `Verification/utilities`; "Copyright (c) 2023 Andrea Basso, Luciano Maino and Giacomo Pope"
- `Verification/isogenies_dim2`; MIT: "Copyright (c) 2023 Pierrick Dartois, Luciano Maino, Giacomo Pope and Damien Robert" 

# Old

## Project structure

Most files are located in the root directory, except the code for supersingular elliptic curves and isogenies in the Montgomery model due to Giacomo Pope. The files that can be run in the root directory are:
- `Kani_isogeny_diamond.py`;
- `benchmark.py`.

They test the computation of 2-isogeny chains between products of 4 elliptic curves, in a context unrelated to SQIsignHD. The other files contain auxiliary code.

The file `benchmark.py` realizes a profiling of the code `Kani_isogeny_diamond.py`. Results are saved in the file `benchmarks/isogeny.cProfile`.

Files to test the computation of two 2-isogeny chains $F_1$ and $\widetilde{F}_2$ defined over a product of elliptic curves and having the same codomain are located in directory `SQISignHD`. The files to run are:
- `tests_dividing_two.py`;
- `Verify.py`.

The first file performs the computation described above, but not in the context of SQISignHD. The second performs the verification of real SQISignHD signatures that have been saved in the file `SQISignHD_executions.txt`.

The parameters used to run `Kani_isogeny_diamond.py`, `benchmark.py` and `tests_dividing_two.py` are located in directory `parameters`. They are stored in `parameters.txt` and have been generated with `parameter_generation.py`. 

## How to run the code?

Launch sagemath in the terminal (command `sage`) and launch 
`load("path to your file from the current directory")`

## Which file to run?

### Kani_isogeny_diamond.py - single 2-isogeny chain

When one runs `Kani_isogeny_diamond.py`, a "random" supersingular elliptic curve $E_1$ is computed, along with an isogeny $\phi: E_1\longrightarrow E_2$ of degree $3^{e_B}$. The integer $e_B$ has been chosen, along with integers $x, y, e_A$ such that
$$x^2+y^2+3^{e_B}=2^{e_A},$$
and the characteristic $p$ is of the form $p:=c\cdot 2^{e_A+2}\cdot 3^{e_B}-1$, so that we have enough accessible torsion defined over $\mathbb{F}_{p^2}$ to compute $\phi$ and evaluate it on $E_1[2^{e_A+2}]$. The parameters $x, y, e_A, e_B, p$... come from the file `parameters/parameters.txt`. 

This way, by Kani's lemma, we can "embed" $\phi$ in a $2^{e_A}$-isogeny $F\in \mbox{End}(E_1^2\times E_2^2)$. Kani_isogeny_diamond.py computes this isogeny as a chain of 2-isogenies in dimension 4 and verifies that $F$ represents $\phi$ by evaluating $F$ on $E_1[2^{e_A+2}]$.

### SQISignHD/tests_dividing_two.py - two 2-isogeny chains

This file runs with the same parameters as the previous one (`parameters/parameters.txt`), except that $F\in \mbox{End}(E_1^2\times E_2^2)$ is decomposed in two parts $F:=F_2\circ F_1$, as in Section 4.3 of the SQIsignHD paper. After computing $F_1$ and $\widetilde{F}_2$, we check that both isogenies have the same codomain, compute $F_2$ as the dual of $\widetilde{F}_2$, evaluate $F_2\circ F_1$ on a torsion point and verify that the image is of the right form.

### SQISignHD/Verify.py - SQISignHD verification

This file runs the verification of SQISignHD (Alogorithms 4 and 5 of the SQISignHD paper) with real SQISignHD parameters and signatures (at NIST-1 level). The signatures have been saved in the file `SQISignHD_executions.txt`. Each line of this file is a signature. Since we have not already implemented an interface between the C signature implementation and the verification, we simply copy pasted the signatures we obtained.

#### Implementation restrictions

Note that the implementation is slightly different from the presentation of the article because:
- We have to recompute the challenge $\varphi: E_A\longrightarrow E_2$ given its kernel.
- For convenience (since we have to recompute $\varphi$ and push points through $\varphi$), we "embed" in dimension 4 the dual of the signature isogeny $\widehat{\sigma}: E_2\longrightarrow E_1$ instead of $\sigma$ itself.
- Our code restricts to signature isogenies of degree $q=\deg(\sigma)$ such that $2^e-q\equiv 5 \mod 8$. In that case, the first 2-isogenies of the chains $F_1$ and $\widetilde{F}_2$ are easier to compute. We leave the implementation of the general case $2^e-q\equiv 1 \mod 4$ to future works.
- Since the code is still experimental, note that signatures are not compressed. In particular, we have more torsion than needed and we have to multiply some torsion points by a power of 2.

#### What does this file do?

Given the kernel of the challenge $\varphi: E_A\longrightarrow E_2$, we compute and evaluate it on the basis $(P_A,Q_A)$ of $E_A[2^f]$ obtained from the signature data. The result $(P_2,Q_2):=(\varphi(P_A),\varphi(Q_A))$ is a basis of $E_2[2^f]$ and we already know the basis $(R_1,R_2):=(\widehat{\sigma}(P_2),\widehat{\sigma}(Q_2))$ from the signature. We can then compute integers $a_1, a_2$ such that $a_1^2+a_2^2+q=2^e$ and compute the embedding $F\in \mbox{End}(E_2^2\times E_1^2)$ of $\widehat{\sigma}$ given by Kani's lemma by decomposing it in two parts $F:=F_2\circ F_1$ as previously. We then evaluate $F(Q,0,0,0)$, where $Q$ is a point of $2^f 3^{f'}$-torsion and check that $F(Q,0,0,0)=([a_1]Q,-[a_2]Q,*,0)$ (up to signs on each component).



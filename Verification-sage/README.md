# Theta-dim4 sagemath library

This library is a `python/sagemath` implementation of the dimension 4 isogeny computations with the Theta model. Only chains of 2-isogenies and Theta models of level 2 are implemented. The library is suited to the computation of isogenies between products of elliptic curves given by Kani's lemma, especially in the context of SQISignHD verification. We recover isogenies in dimension 1 given their evaluation on torsion points. Hence, the library could also be used to test SIDH attacks.

## Disclaimer

This library is experimental. Further tests and optimizations are still needed. 

Due to already implemented optimizations, the code is significantly different from what is described in Appendix F of the SQISignHD paper [1]. Another paper dedicated to dimension 4 $2$-isogeny computations in the Theta model is in preparation to describe the algorithms implemented in this library [2]. 

Note that this library heavily relies on dimension 2 isogeny computations in the Theta-model as described in [3].

## Requirement

You should have `python` and `sagemath` installed on your computer. Instructions to download `sagemath` can be found here https://doc.sagemath.org/html/en/installation/index.html.

## How to run the code?

Launch sagemath in the terminal (type command `sage`) and type 
`load("path to 'file' from your current directory")`
where `file` is the file containing the code you want to run.

## <a name="Generic"></a> Generic tests

The script in `Tests.py` computes several instanciations of dimension 4 2-isogeny chains derived from Kani's lemma given as endomorphisms:

$$F:=\left(\begin{matrix} a_1 & a_2 & \widehat{\sigma} & 0 \\ -a_2 & a_1 & 0 & \widehat{\sigma} \\ -\sigma & 0 & a_1 & -a_2 \\ 0 & -\sigma & a_2 & a_1 \end{matrix}\right)\in End(E_1^2\times E_2^2),$$

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

### Chain with the full torsion (as in the SIDH attacks)

When we can access the full $2^{e_A+2}$-torsion on $E_1$ and $E_2$, the whole isogeny chain $F: E_1^2\times E_2^2\longrightarrow E_1^2\times E_2^2$ can be computed directly. The function `test_kani_endomorphism` tests this direct chain computation for specified sets of parameters. 

INPUT:
- index: specifies the index of the set of parameters in the list extracted from the file `parameters/parameters.txt`
	if $\ell_B=3$ or `parameters/parameters_7.txt` if $\ell_B=7$.
- $\ell_B$: prime specifying the degree of the embedded isogeny $\sigma$ ($\deg(\sigma)=\ell_B^{e_B}$), $\ell_B=3$ or $7$ ($7$ by default).

OUTPUT: an object of the class `KaniEndo` (imported from `isogenies/Kani_endomorphism.py`) representing a dimension 4 $2$-isogeny chain derived from Kani's lemma.

EXAMPLE:
```
sage: load("Tests.py")
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
sage: F=test_kani_endomorphism(1,3)
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

To automatically test all sets of parameters in `parameters/parameters.txt` (except the set of index 0), set: 
```python
test_endomorphism_3=True
```
in the file `Tests.py`. To automatically test all sets of parameters in `parameters/parameters_7.txt`, set: 
```python
test_endomorphism_7=True
```
in the file `Tests.py`.

### Chain with half the necessary torsion (as in SQISignHD)

When we cannot access the full $2^{e_A+2}$-torsion on $E_1$ and $E_2$, the computation of the isogeny chain $F: E_1^2\times E_2^2\longrightarrow E_1^2\times E_2^2$ can has to be divided in two, as in Section 4.3 of the SQISignHD paper. Namely, we compute two isogeny chains $F_1: E_1^2\times E_2^2\longrightarrow C$ and $\widetilde{F_2}: E1^2\times E2^2\longrightarrow C$ such that $F=F_2\circ F_1$.

The function `test_kani_endomorphism_half` tests this chain computation for specified sets of parameters. Note that, in practice, the full torsion is available with the tested sets of parameters but only half is used.

INPUT:
- index: specifies the index of the set of parameters in the list extracted from the file `parameters/parameters.txt`
	if $\ell_B=3$ or `parameters/parameters_7.txt` if $\ell_B=7$.
- $\ell_B$: prime specifying the degree of the embedded isogeny $\sigma$ ($\deg(\sigma)=\ell_B^{e_B}$), $\ell_B=3$ or $7$ ($7$ by default).

OUTPUT: an object of the class `KaniEndoHalf` (imported from `isogenies/Kani_endomorphism.py`) representing a dimension 4 $2$-isogeny chain derived from Kani's lemma.

EXAMPLE:
```
sage: load("Tests.py")
sage: F=test_kani_endomorphism_half(2)
Testing KaniEndoHalf with parameters:
 - Prime characteristic p = 9 * 2**35 * 7**5 - 1
 - Degree of the embedded isogeny sigma q = 7**5
 - a1 = 35037
 - a2 = 85804
 - Length of the dimension 4 2-isogeny = 33
 - Used available torsion = 2**19
Parameter import: 0.03811001777648926 s
Random walk: 0.024472951889038086 s
Generation of sigma: 0.0059468746185302734 s
Generation and evaluation of the torsion basis: 0.0032341480255126953 s
Dimension 4 endomorphism: 0.28946805000305176 s
Is evaluation correct?
True
Time evaluation: 0.008486032485961914 s
sage: F=test_kani_endomorphism_half(2,3)
Testing KaniEndoHalf with parameters:
 - Prime characteristic p = 1 * 2**34 * 3**13 - 1
 - Degree of the embedded isogeny sigma q = 3**13
 - a1 = -52562
 - a2 = 39123
 - Length of the dimension 4 2-isogeny = 32
 - Used available torsion = 2**18
Parameter import: 0.00571894645690918 s
Random walk: 0.03211402893066406 s
Generation of sigma: 0.008935213088989258 s
Generation and evaluation of the torsion basis: 0.006306886672973633 s
Dimension 4 endomorphism: 0.24866199493408203 s
Is evaluation correct?
True
Time evaluation: 0.00876617431640625 s
```

To automatically test all sets of parameters in `parameters/parameters.txt` (except the set of index 0), set: 
```python
test_half_endomorphism_3=True
```
in the file `Tests.py`. To automatically test all sets of parameters in `parameters/parameters_7.txt`, set: 
```python
test_half_endomorphism_7=True
```
in the file `Tests.py`.

### Implementation restrictions

When computing dimension 4 isogenies with the Theta model, extra care is needed to compute gluing isogenies (isogenies starting from a product of abelian varieties of dimension less than 4). In most cases, when we compute an isogeny chain derived from Kani's lemma, we only encounter gluings at the beginning. However, it may happen in small characteristic that an isogeny splits into a product of abelian varieties in the middle of the chain. In this case, we have to compute a gluing afterwards and this step generally fails since we can only compute gluings when we know their location in advance in the chain. Our code then returns a NotImplementedError.    

In particular, the set of parametres of index 0 in `parameters/parameters.txt` always strikes a failure.

```
sage: F=test_kani_endomorphism(0,3)
Testing KaniEndo with parameters:
 - Prime characteristic p = 1 * 2**18 * 3**5 - 1
 - Degree of the embedded isogeny sigma q = 3**5
 - a1 = 238
 - a2 = -93
 - Length of the dimension 4 2-isogeny = 16
Parameter import: 0.0007309913635253906 s
Random walk: 0.015373945236206055 s
Generation of sigma: 0.004636049270629883 s
Generation and evaluation of the torsion basis: 0.0027260780334472656 s
Zero dual theta constants.
Random symplectic base change is being used for duplication.
Doublings are more costly than expected.

[...]

NotImplementedError: The codomain of this 2-isogeny could not be computed.
We may have encountered a product of abelian varieties
somewhere unexpected along the chain.
This is exceptionnal and should not happen in larger characteristic.
sage: F=test_kani_endomorphism_half(0,3)
Testing KaniEndoHalf with parameters:
 - Prime characteristic p = 1 * 2**18 * 3**5 - 1
 - Degree of the embedded isogeny sigma q = 3**5
 - a1 = 238
 - a2 = -93
 - Length of the dimension 4 2-isogeny = 16
 - Used available torsion = 2**10
Parameter import: 0.13568615913391113 s
Random walk: 0.015301942825317383 s
Generation of sigma: 0.005235910415649414 s
Generation and evaluation of the torsion basis: 0.0023212432861328125 s

[...]

NotImplementedError: The codomain of this 2-isogeny could not be computed.
We may have encountered a product of abelian varieties
somewhere unexpected along the chain.
This is exceptionnal and should not happen in larger characteristic.
```

## Verifying real SQISignHD signatures

The file `Verify_SQISignHD.py` runs the verification of SQISignHD (Alogorithms 4 and 5 of the SQISignHD paper [1]) with real SQISignHD parameters and signatures (at NIST-1 level). Ten signatures have been saved in the file `SQIsignHD_data/SQISignHD_executions.txt`. Each line of this file is a signature. Since we have not already implemented an interface between the C signature implementation and the verification, we simply copy pasted the signatures we obtained. You can test new signatures by running new signatures and copy pasting them in the file `SQIsignHD_data/SQISignHD_executions.txt`.

To run all the signatures contained in `SQIsignHD_data/SQISignHD_executions.txt`, simply type:
```python
load("Verify_SQISignHD.py")
```

### Verifying a specific signature

Once `Verify_SQISignHD.py` is loaded, one can also verify a specific signature with the function `verify` by specifying its index in the list `L_exec` (between 0 and 9).

```
sage: load("Verify_SQISignHD.py")
===========================================================
Testing 10 instances of SQISignHD verification parameters:
 - Prime characteristic p = 13 * 2**126 * 3**78 - 1
 - Length of the dimension 4 2-isogeny = 142
 - Used available torsion = 2**73
===========================================================

Test 1
Signature data parsing time 0.19607090950012207 s
Challenge computation time 0.18125009536743164 s
Endomorphism F=F2*F1 computation time 1.4973571300506592 s
Do F1 and F2_dual have the same codomain?
True
Codomain matching verification time 0.00028967857360839844 s
Time to find a point of order 2^f*3^fp 0.012469053268432617 s
Time point evaluation 0.05667710304260254 s
Is the point evaluation correct ?
True
Total verification time 1.754267930984497 s

[...]

sage: data=verify(L_exec[0])
Signature data parsing time 0.08100771903991699 s
Challenge computation time 0.17457008361816406 s
Endomorphism F=F2*F1 computation time 1.55623197555542 s
Do F1 and F2_dual have the same codomain?
True
Codomain matching verification time 0.00023794174194335938 s
Time to find a point of order 2^f*3^fp 0.01227712631225586 s
Time point evaluation 0.0759129524230957 s
Is the point evaluation correct ?
True
Total verification time 1.825559139251709 s
```

### What does this file do?

Given the kernel of the challenge $\varphi: E_A\longrightarrow E_2$, we compute and evaluate it on the basis $(P_A,Q_A)$ of $E_A[2^f]$ obtained from the signature data. The result $(P_2,Q_2):=(\varphi(P_A),\varphi(Q_A))$ is a basis of $E_2[2^f]$ and we already know the basis $(R_1,R_2):=(\widehat{\sigma}(P_2),\widehat{\sigma}(Q_2))$ from the signature. We can then compute integers $a_1, a_2$ such that $a_1^2+a_2^2+q=2^e$ and compute the embedding $F\in \mbox{End}(E_2^2\times E_1^2)$ of $\widehat{\sigma}$ given by Kani's lemma by decomposing it in two parts $F:=F_2\circ F_1$ as previously. We then evaluate $F(Q,0,0,0)$, where $Q$ is a point of $2^f 3^{f'}$-torsion and check that $F(Q,0,0,0)=([a_1]Q,-[a_2]Q,*,0)$ (up to signs on each component).

### Implementation restrictions

Note that the implementation is slightly different from the presentation of the article because:
- We have to recompute the challenge $\varphi: E_A\longrightarrow E_2$ given its kernel.
- For convenience (since we have to recompute $\varphi$ and push points through $\varphi$), we "embed" in dimension 4 the dual of the signature isogeny $\widehat{\sigma}: E_2\longrightarrow E_1$ instead of $\sigma$ itself.
- For now, only isogenies of degree $q=\deg(\sigma)$ such that $2^e-q\equiv 5 \mod 8$ are outputted by the signature protocol but our verification procedure should work with any degree $q$ such that $2^e-q\equiv 1 \mod 4$, as specified in the SQISignHD paper [1].
- Since the code is still experimental, note that signatures are not compressed. In particular, we have more torsion than needed and we have to multiply some torsion points by a power of 2.

## Organization of the library

The main test files `Tests.py` and `Verification_SQIsignHD.py` are located in the main directory `Verification-sage` of this library.

`Verification-sage` contains several subdirectories:
- `basis_change` contains code for changing level 2 Theta structures in dimension 2 and 4 useful for the computation of gluing and splitting isogenies (in the beginning and in the end of the chain).
- `isogenies` contains code to compute 2-isogenies and chains of 2-isogenies in dimension 4 in the Theta model.
- `isogenies_dim2` contains code to compute 2-isogenies and chains of 2-isogenies in dimension 2 in the Theta model, based on the implementation of [3]. This code is useful for the first steps of dimension 4 2-isogeny chains (involving gluings).
- `montgomery_isogenies` contains code to compute dimension 1 isogeny on the Kummer line faster than `sagemath` using $x$-only arithmetic. These files are due to Giacomo Pope.
- `parameters` contains the parameters mentionned in [Generic tests](#Generic).
- `SQISignHD_data` contains real SQISignHD signatures.
- `theta_structures` contains code for different models in dimension 1, 2 and 4 and translations between those models (Montgomery model on elliptic curves (dimension 1) or product of elliptic curves, level 2 Theta models in dimensions 1, 2 and 4).
- `utilities` contains several useful functions to work on supersingular elliptic curves faster than with the `sagemath` generic functions (discrete logarithms, Weil pairings...).

## License

SQIsignHD is licensed under Apache-2.0. See LICENSE and NOTICE in the root directory.

Third party code is used in this directory (`Verification-sage`):

- `montgomery_isogenies`; MIT: "Copyright (c) 2023 Giacomo Pope"
- `utilities`; "Copyright (c) 2023 Andrea Basso, Luciano Maino and Giacomo Pope"
- `isogenies_dim2`; MIT: "Copyright (c) 2023 Pierrick Dartois, Luciano Maino, Giacomo Pope and Damien Robert" 

## References

[1] P. Dartois, A. Leroux, D. Robert, and B. Wesolowski. SQIsignHD: New Dimensions in Cryptography. Cryptology ePrint Archive, Paper 2023/436. 2023. url: https://eprint.iacr.org/2023/436.

[2] P. Dartois. Fast computation of 2-isogenies in dimension 4 with the Theta model and cryptographic applications. In preparation. 2024.

[3] P. Dartois, L. Maino, G. Pope, and D. Robert. An Algorithmic Approach to (2,2)-isogenies in the Theta Model and Applications to Isogeny-based Cryptography. Cryptology ePrint Archive, Paper 2023/1747. 2023. url: https://eprint.iacr.org/2023/1747.


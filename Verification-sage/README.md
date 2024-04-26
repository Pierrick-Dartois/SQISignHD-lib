# Theta-dim4 sagemath library

This library is a `python/sagemath` implementation of the dimension 4 isogeny computations with the Theta model. Only chains of 2-isogenies and Theta models of level 2 are implemented. The library is suited to the computation of isogenies between products of elliptic curves given by Kani's lemma, especially in the context of SQISignHD verification. We recover isogenies in dimension 1 given their evaluation on torsion points. Hence, the library could also be used to test SIDH attacks.

## Disclaimer

This library is experimental. Further tests and optimizations are still needed. 

Due to already implemented optimizations, the code is significantly different from what is described in Appendix F of the SQISignHD paper [1]. Another paper dedicated to dimension 4 $2$-isogeny computations in the Theta model is in preparation to describe the algorithms implemented in this library [2]. 

Note that this library heavily relies on dimension 2 isogeny computations in the Theta-model as described in [3].

## Requirement

You should have `python 3` and `sagemath` (version 10 at least) installed on your computer. Instructions to download `sagemath` can be found here https://doc.sagemath.org/html/en/installation/index.html.

## How to run the code?

Start a terminal and use the `cd` command to locate in the `Verification-sage` directory. From there, you can access two command line interfaces (CLI):
- `Tests.py` running generic tests on 4-dimensional 2-isogeny chains derived from Kani's lemma;
- `Verify_SQISignHD.py` to verify real SQISignHD signatures.

To use the CLI, type:

`sage <Tests.py/Verify_SQISignHD.py> <arguments>`

More details on the `<arguments>` of each CLI are provided in the following. 

## <a name="Generic"></a> Generic tests

The script in `Tests.py` computes instanciations of 4-dimensional 2-isogeny chains derived from Kani's lemma given as endomorphisms:

$$F:=\left(\begin{matrix} a_1 & a_2 & \widehat{\sigma} & 0 \\ 
-a_2 & a_1 & 0 & \widehat{\sigma} \\ 
-\sigma & 0 & a_1 & -a_2 \\ 
0 & -\sigma & a_2 & a_1 \end{matrix}\right)\in \mathrm{End}(E_1^2\times E_2^2),$$

where:
- $E_1$ is a "random" supersingular elliptic curve, generated as the codomain of a random isogeny walk of degree $\simeq p$ starting from the supersingular elliptic curve $E_0$ of $j$-invariat $1728$;
- $\sigma: E_1\longrightarrow E_2$ is an isogeny of degree $\ell_B^{e_B}$ where $\ell_B$ is an odd prime ($\ell_B=3$ or $7$ here); 
- Integers $a_1, a_2, e_A$ have been chosen such that:
$$a_1^2+a_2^2+\ell_B^{e_B}=2^{e_A};$$
- The characteristic $p$ is of the form $p:=f\cdot 2^{f_A}\cdot \ell_B^{f_B}-1$, with $f_A\geq e_A+2$ and $f_B\geq e_B$ so that we have enough accessible torsion defined over $\mathbb{F}_{p^2}$ to compute $\sigma$ and evaluate it on $E_1[2^{e_A+2}]$.

This script includes two different ways of computing such a dimension 4 isogeny:
- the command line instruction `--KaniEndo` tests instances when we have the full torsion available.
- the function `--KaniEndoHalf` tests instances when only half the torsion is available (as in SQISignHD, see Section 4.3 of the paper).

### Parameters 

Different set of parameters can be found in the subdirectory `parameters`. The files labeled `parameters/parameters_3.txt` and `parameters/parameters_7.txt` contain lists of parameters for the values $\ell_B=3$ and $\ell_B=7$ respectively. Every line of those files contains a list:

$$[e_A, e_B, a_1, a_2, f, f_A, f_B, p(, m)]$$

such that:
- $p:=f\cdot 2^{f_A}\cdot \ell_B^{f_B}-1$;
- $f_A\geq e_A+2$;
- $f_B\geq e_B$;
- $m=\max(v_2(a_1),v_2(a_2))$ (this parameter is not present in the list when $\ell_B=3$ since it is always $1$ in this case).

These parameters have been generated with the functions `find_param` and `find_param_gen` of `parameters/parameter_generation.py`. 

#### Displaying parameters 

You can display available parameters with the command:

```
sage Tests.py --display
```

or in short:

```
sage Tests.py -d
```

This command will display all saved parameters with their index in the list. Those indices will be used as reference to execute tests on a specific set of parameters.

If you want to restrict the displayed lists of parameters, type:

```
sage Tests.py -d -l_B=<value>
```

to display all parameters with a given `<value>` of $\ell_B$ or type:

```
sage Tests.py -d -l_B=<value l_B> -i=<value index>
```

to display the list of parameters with a given `<value l_B>` of $\ell_B$ indexed with `<value index>`.

#### Creating new parameters

If you want to create new lists of parameters, you can type:

```sage --add_params -l_B=<value l_B> -e_A=<value e_A>```

this command will create lists of parameters with $\ell_B=$`<value l_B>`, $e_A=$`<value e_A>` and $e_B$ varying in the range $2^{e_A/2}\leq \ell_B^{e_B}\leq 2^{e_A}$. Depending on the specified values, multiple parameters or no parameter can be found. The new parameters are saved in an existing file (after the previously saved parameters) or in a newly created file `parameters/parameters_<value l_B>.txt`. The newly created parameters are displayed along with their index in the file.

EXAMPLE:
```
% sage Tests.py --add_params -l_B=11 -e_A=66
===========================================
New parameters with second prime l_B=11.
===========================================

 - Index in the list of parameters = 0
 - Prime characteristic p = 1 * 2**70 * 11**11 - 1
 - Degree of the embedded isogeny sigma q = 11**11
 - a1 = 1440875693
 - a2 = -8468226098
 - m = max(v_2(a1),v_2(a2)) = 1
 - Length of the 4-dimensional 2-isogeny = 66
 ```

The parameter generation is subexponential and may take a while. We do not recommend to select `<value e_A>` bigger than $200$.

### Chain with the full torsion available (as in the SIDH attacks)

When we can access the full $2^{e_A+2}$-torsion on $E_1$ and $E_2$, the whole $2$-isogeny chain $F: E_1^2\times E_2^2\longrightarrow E_1^2\times E_2^2$ can be computed directly. The command line instruction `--KaniEndo` tests this direct chain computation for specified sets of parameters. 

Type:

```
sage Tests.py --KaniEndo
```

to test all saved lists of parameters. 

Type:

```
sage Tests.py --KaniEndo -l_B=<value l_B>
```

to test all saved lists of parameters with $\ell_B=$`<value l_B>`.

Type:

```
sage Tests.py --KaniEndo -l_B=<value l_B> -i=<value index>
```

to test the list of parameters with $\ell_B=$`<value l_B>` indexed by `<value index>`.

#### For experimented users only

You can test a specific set of parameters chosen manually. If the set of parameters is not well chosen, this could lead to errors.

```
sage Tests.py --KaniEndo -l_B=<value l_B> -e_A=<value e_A> -e_B=<value e_B> -a1=<value a1> -a2=<value a2> -f=<value f> -f_A=<value f_A> -f_B=<value f_B> -p=<value p> [optional: -m=<value m>]
```

EXAMPLE:
```
% sage Tests.py --KaniEndo -l_B=7 -e_A=17 -e_B=3 -a1=123 -a2=-340 -f=3 -f_A=20 -f_B=3 -p=1078984703
Testing KaniEndo with parameters:
 - Prime characteristic p = 3 * 2**20 * 7**3 - 1
 - Degree of the embedded isogeny sigma q = 7**3
 - a1 = 123
 - a2 = -340
 - m = max(v_2(a1),v_2(a2)) = 2
 - Length of the dimension 4 2-isogeny = 17
Setup: 0.04501986503601074 s
Random walk: 0.015614986419677734 s
Generation of sigma: 0.004041910171508789 s
Generation and evaluation of the torsion basis: 0.002375364303588867 s
Strategy computation: 0.0004138946533203125 s
Dimension 4 endomorphism: 0.16916584968566895 s
Is evaluation correct?
True
```

### Chain with half the necessary torsion (as in SQISignHD)

When we cannot access the full $2^{e_A+2}$-torsion on $E_1$ and $E_2$, the computation of the 2-isogeny chain $F: E_1^2\times E_2^2\longrightarrow E_1^2\times E_2^2$ can has to be divided in two, as in Section 4.3 of the SQISignHD paper. Namely, we compute two isogeny chains $F_1: E_1^2\times E_2^2\longrightarrow C$ and $\widetilde{F_2}: E1^2\times E2^2\longrightarrow C$ such that $F=F_2\circ F_1$. The command line instruction `--KaniEndoHalf` tests this computation for specified sets of parameters. 

Type:

```
sage Tests.py --KaniEndoHalf
```

to test all saved lists of parameters. 

Type:

```
sage Tests.py --KaniEndoHalf -l_B=<value l_B>
```

to test all saved lists of parameters with $\ell_B=$`<value l_B>`.

Type:

```
sage Tests.py --KaniEndoHalf -l_B=<value l_B> -i=<value index>
```

to test the list of parameters with $\ell_B=$`<value l_B>` indexed by `<value index>`.

#### For experimented users only

You can test a specific set of parameters chosen manually. If the set of parameters is not well chosen, this could lead to errors.

```
sage Tests.py --KaniEndoHalf -l_B=<value l_B> -e_A=<value e_A> -e_B=<value e_B> -a1=<value a1> -a2=<value a2> -f=<value f> -f_A=<value f_A> -f_B=<value f_B> -p=<value p> [optional: -m=<value m>]
```

EXAMPLE:
```
% sage Tests.py --KaniEndoHalf -l_B=7 -e_A=17 -e_B=3 -a1=123 -a2=-340 -f=3 -f_A=20 -f_B=3 -p=1078984703
Testing KaniEndoHalf with parameters:
 - Prime characteristic p = 3 * 2**20 * 7**3 - 1
 - Degree of the embedded isogeny sigma q = 7**3
 - a1 = 123
 - a2 = -340
 - m = max(v_2(a1),v_2(a2)) = 2
 - Length of the dimension 4 2-isogeny = 17
 - Used available torsion = 2**11
Setup: 0.04306387901306152 s
Random walk: 0.015376091003417969 s
Generation of sigma: 0.004205942153930664 s
Generation and evaluation of the torsion basis: 0.002997159957885742 s
Computation of strategies: 0.00025773048400878906 s
Dimension 4 endomorphism: 0.20006227493286133 s
Is evaluation correct?
True
Time evaluation: 0.00693202018737793 s
```

### Implementation restrictions

When computing dimension 4 isogenies with the Theta model, extra care is needed to compute gluing isogenies (isogenies starting from a product of abelian varieties of dimension less than 4). In most cases, when we compute an isogeny chain derived from Kani's lemma, we only encounter gluings at the beginning. However, it may happen in small characteristic that an isogeny splits into a product of abelian varieties in the middle of the chain. In this case, we have to compute a gluing afterwards and this step generally fails since we can only compute gluings when we know their location in advance in the chain. Our code then returns a NotImplementedError.    

In particular, the set of parameters with $\ell_B=3$ and index 0 always strikes a failure.

```
% sage Tests.py --KaniEndo -l_B=3 -i=0
Testing KaniEndo with parameters:
 - Prime characteristic p = 1 * 2**18 * 3**5 - 1
 - Degree of the embedded isogeny sigma q = 3**5
 - a1 = 238
 - a2 = -93
 - m = max(v_2(a1),v_2(a2)) = 1
 - Length of the dimension 4 2-isogeny = 16
Setup: 0.04303693771362305 s
Random walk: 0.018531084060668945 s
Generation of sigma: 0.004489898681640625 s
Generation and evaluation of the torsion basis: 0.002763986587524414 s
Strategy computation: 0.0003800392150878906 s
[...]
NotImplementedError: The codomain of this 2-isogeny could not be computed.
We may have encountered a product of abelian varieties
somewhere unexpected along the chain.
This is exceptionnal and should not happen in larger characteristic.
```

## Verifying real SQISignHD signatures

The file `Verify_SQISignHD.py` runs the verification of SQISignHD (Alogorithms 4 and 5 of the SQISignHD paper [1]) with real SQISignHD parameters and signatures (at NIST-1 level). Ten signatures have already been saved in the file `SQIsignHD_data/SQISignHD_executions.txt`.

### What does this script do?

Given the kernel of the challenge $\varphi: E_A\longrightarrow E_2$, we compute and evaluate it on the basis $(P_A,Q_A)$ of $E_A[2^f]$ obtained from the data from the signature file. The result $(P_2,Q_2):=(\varphi(P_A),\varphi(Q_A))$ is a basis of $E_2[2^f]$ and we already know the basis $(R_1,R_2):=(\widehat{\sigma}(P_2),\widehat{\sigma}(Q_2))$ from the signature. We can then compute integers $a_1, a_2$ such that $a_1^2+a_2^2+q=2^e$ and compute the embedding $F\in \mbox{End}(E_2^2\times E_1^2)$ of $\widehat{\sigma}$ given by Kani's lemma by decomposing it in two parts $F:=F_2\circ F_1$ as previously. We then evaluate $F(Q,0,0,0)$, where $Q$ is a point of $2^f 3^{f'}$-torsion and check that $F(Q,0,0,0)=([a_1]Q,-[a_2]Q,*,0)$ (up to signs on each component).

### Verifying signatures in the sample

To test all 10 signatures contained in the sample `SQIsignHD_data/SQISignHD_executions.txt`, type:

```
sage Verify_SQISignHD.py --verify_sample
```

or in short format:

```
sage Verify_SQISignHD.py -vs
```

You can also verify only one signature of index `<value index>` between 0 and 9:

```
sage Verify_SQISignHD.py -vs -i=0
```

EXAMPLE:
```
% sage Verify_SQISignHD.py --verify_sample
===============================================================
Testing 10 instances of SQISignHD verification with parameters:
 - Prime characteristic p = 13 * 2**126 * 3**78 - 1
 - Length of the dimension 4 2-isogeny = 142
 - Used available torsion = 2**73
===============================================================

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

### Verifying a generated signature

Instructions to generate signatures may be found in the `README.md` file of the `Signature` library. We strongly advise to save the generated signature in the `Verification-sage/SQIsignHD_data` subdirectory as instructed in the `Signature/README.md`:

```
./src/sqisignhd/ref/lvl1/test/sqisign_test_sqisignhd_lvl1 > ../../Verification-sage/SQISignHD_data/<specify file name>.txt
```

Otherwise, you have to specify a relative path from the `Verification-sage` to verify your signature.

Once, you have generated a signature, you may verify it from `Verification-sage` with the command:

```
sage Verify_SQISignHD.py --verify_one_signature <path to signature file>
```

or in short format:

```
sage Verify_SQISignHD.py -vo <path to signature file>
```

EXAMPLE:
```
% sage Verify_SQISignHD.py --verify_one_signature SQISignHD_data/signature.txt
===========================================================
Verifying one SQISignHD signature with parameters:
 - Prime characteristic p = 13 * 2**126 * 3**78 - 1
 - Length of the 4-dimensional 2-isogeny chain = 142
 - Used available torsion = 2**73
===========================================================

Signature data parsing time 0.20653986930847168 s
Challenge computation time 0.16781020164489746 s
Time to compute the strategies 0.0027136802673339844 s
Endomorphism F=F2*F1 computation time 1.514930248260498 s
Do F1 and F2_dual have the same codomain?
True
Codomain matching verification time 0.00026869773864746094 s
Time to find a point of order 2^f*3^fp 0.006500244140625 s
Time point evaluation 0.06669378280639648 s
Is the point evaluation correct ?
True
Total verification time 1.7653708457946777 s
```

### Implementation restrictions

Note that the implementation is slightly different from the presentation of the article because:
- We have to recompute the challenge $\varphi: E_A\longrightarrow E_2$ given its kernel.
- For convenience (since we have to recompute $\varphi$ and push points through $\varphi$), we "embed" in dimension 4 the dual of the signature isogeny $\widehat{\sigma}: E_2\longrightarrow E_1$ instead of $\sigma$ itself.
- For now, only isogenies of degree $q=\deg(\sigma)$ such that $2^e-q\equiv 5 \mod 8$ are outputted by the signature protocol but our verification procedure should work with any degree $q$ such that $2^e-q\equiv 1 \mod 4$, as specified in the SQISignHD paper [1].
- Since the code is still experimental, note that signatures are not compressed. In particular, the deterministic basis $(P_A,Q_A)$ of $E_A[2^f]$ is part of the signature and we have more torsion than needed (we have to multiply torsion points by a power of 2).

## Organization of the library

The main test files `Tests.py` and `Verify_SQIsignHD.py` are located in the main directory `Verification-sage` of this library.

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


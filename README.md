# SQISignHD library (version 2.0)

This library contains an implementation of F-SQISignHD (Fast Version of the Short Quaternion Isogeny Signature in Higher Dimension), as described in Pierrick Dartois' PhD thesis [1, Chapter 3]. 

Copyright (c) 2025, Pierrick Dartois, Antonin Leroux, Damien Robert and Benjamin Wesolowski.

Note that this library differs significantly from the first version of FastSQIsignHD described in the SQISignHD paper [2]. In particular, the new signature implementation follows from ideal-to-isogeny translation algorithms introduced for the SQIsign2D-West protocol [3] which is now the reference SQIsign NIST candidate [4].

This library contains two independent sublibraries:

- `Signature`: a `C` implementation of the signature which is a light modification of the SQISign2D-West implementation <https://github.com/SQISign/the-sqisign>. The original work has been implemented by Luca De Feo, Antonin Leroux and Giacomo Pope. Antonin Leroux in particular proposed an new implementation of SQIsignHD from which this work is derived.
- `Verification`: that includes a script for the F-SQIsignHD verification written in `sagemath` which imports functions from the `Theta_dim4` submodule <https://github.com/Pierrick-Dartois/Theta_dim4>, an experimental `sagemath` implementation of dimension 4 $2$-isogeny chains in the Theta model of level $2$ described in [5].

## Running the signature

### Requirements

- CMake (version 3.5 or later)
- C99-compatible compiler
- GMP (version 6.1.2 or later)

### Build and check

For the development build, go the `Signature` subfolder and type the following instructions in your terminal

- `mkdir -p build`
- `cd build`
- `cmake -DSQISIGN_BUILD_TYPE=ref ..`
- `make`
- `make test` (optional, for testing)

These commands create 3 executables `build/src/sqisignhd/ref/lvl[1,3,5]/test/sqisign_test_sqisignhd_lvl[1,3,5]` (one for each NIST security level 1, 3 or 5). To run, one of these executables, type

`./src/sqisignhd/ref/lvl[1,3,5]/test/sqisign_test_sqisignhd_lvl[1,3,5]`

Note that this instruction will only work while in the `build` folder.

Each one of these executables tests (for a given NIST security level) the key generation and signature algorithms and runs a benchmark based on 100 key generation-signature loop. The resulting public key and signature pairs will be printed out in the files `Verification/Data/Public_keys_lvl[1,3,5].txt` and `Verification/Data/Signatures_lvl[1,3,5].txt` (where paths are indicated from the root of the SQIsignHD library). Public keys and signatures are ordered by pairs in the files. The signatures from these files can then be checked with the [verification script](#Verif).

### Build options

CMake build options can be specified with `-D<BUILD_OPTION>=<VALUE>`.

An optimised executable can be built by running
`cmake -DCMAKE_BUILD_TYPE=Release ..`

#### ENABLE_GMP_BUILD

If set to `OFF` (by default), the gmp library on the system is dynamically linked.
If set to `ON`, a custom gmp library is linked, which is built as part of the overall build process. 

In the latter case, the following further options are available:
- `ENABLE_GMP_STATIC`: Does static linking against gmp. The default is `OFF`.
- `GMP_BUILD_CONFIG_ARGS`: Provides additional config arguments for the gmp build (for example `--disable-assembly`). By default, no config arguments are provided.

#### CMAKE_BUILD_TYPE

Can be used to specify special build types. The options are:

- `Release`: Builds with optimizations enabled and assertions disabled.
- `Debug`: Builds with debug symbols.
- `ASAN`: Builds with AddressSanitizer memory error detector.
- `MSAN`: Builds with MemorySanitizer detector for uninitialized reads.
- `LSAN`: Builds with LeakSanitizer for run-time memory leak detection.
- `UBSAN`: Builds with UndefinedBehaviorSanitizer for undefined behavior detection.

The default build type uses the flags `-O3 -Wstrict-prototypes -Wno-error=strict-prototypes -fvisibility=hidden -Wno-error=implicit-function-declaration -Wno-error=attributes`. (Notice that assertions remain enabled in this configuration, which harms performance.)

### For developpers

Further instructions addressed to developpers of the `C` signature code base may be found in `Signature/README.md`. This file has been written for the SQIsign2D-West implementation so some instructions may be irrelevant to SQIsignHD or outdated. We decline all responsability on the content of this file.

## <a name="Verif"></a> Running the verification

### Requirements

- Python 3
- SageMath (version 10.0 or later)

### Testing and benchmarking the verification

To test the verification implemented in `sagemath`, go to the `Verification` subfolder and type in your terminal

`sage Verify.py [-lvl/--level] [-test] [-bench] [-i/--index] [-n/--n_samples]`

with the following options (in brackets):

- `-lvl/--level` (compulsory) specifies the NIST security level to be tested (1, 3 or 5). For level 1, type `-lvl=1` or `--level=1`.
- `-test` or `-bench` (you must type one and only one of them). With the `-test` option, all execution details are printed so this option is suitable for testing. With the `-bench` option, minimal details are printed. This is meant to run clean benchmarks.
- `-i/--index` (optional) specifies an index (in the range $\{0,\cdots, 99\}$) of a single public key and signature pair to verify in the files `Verification/Data/Public_keys_lvl[1,3,5].txt` and `Verification/Data/Signatures_lvl[1,3,5].txt`. Type `-i=0` or `--index=0` to choose the pair of index 0.
- `-n/--n_samples` (optional) specifies the number of single public key and signature pairs to verify (in the range $\{0,\cdots, 100\}$) in the files `Verification/Data/Public_keys_lvl[1,3,5].txt` and `Verification/Data/Signatures_lvl[1,3,5].txt`. Type `-n=10` or `--n_samples=10` to run 10 samples.

Note that by default, the script verifies the 100 public key and signature pairs from the files `Verification/Data/Public_keys_lvl[1,3,5].txt` and `Verification/Data/Signatures_lvl[1,3,5].txt`. This may take a while so the options `-i/--index` and `-n/--n_samples` can be used to restrict the number of public key and signature pairs to test.

### Examples

```
% sage Verify.py -lvl=1 -test -n=10

##########################################
# Testing level 1 SQIsignHD verification #
# p = 5*2**248 - 1                       #
##########################################

Verifying response number 0

Starting verification:
Commitment and public key basis generation time 0.04948616027832031 s
Challenge and challenge basis computation time 0.11446309089660645 s
Response image recovery time 0.07320880889892578 s
4-dimensional isogeny computation time 1.5313968658447266 s
Do F1 and F2_dual have the same codomain?
True
Codomain matching verification time 0.00016999244689941406 s
Is the point evaluation correct ?
True
Time image verification 0.1070868968963623 s

End verification. Total verification time 1.8760039806365967 s

[...]

Verifying response number 9

Starting verification:
Commitment and public key basis generation time 0.006487846374511719 s
Challenge and challenge basis computation time 0.10206913948059082 s
Response image recovery time 0.0702669620513916 s
4-dimensional isogeny computation time 1.5113840103149414 s
Do F1 and F2_dual have the same codomain?
True
Codomain matching verification time 0.0001761913299560547 s
Is the point evaluation correct ?
True
Time image verification 0.1116790771484375 s

End verification. Total verification time 1.802258014678955 s


All tests passed.

Average verification time 1.8279810905456544 s
```

```
% sage Verify.py -lvl=1 -test -i=98

sage Verify.py -lvl=1 -test -i=98
##########################################
# Testing level 1 SQIsignHD verification #
# p = 5*2**248 - 1                       #
##########################################

Verifying response number 98

Starting verification:
Commitment and public key basis generation time 0.027514934539794922 s
Challenge and challenge basis computation time 0.0976099967956543 s
Response image recovery time 0.06902694702148438 s
4-dimensional isogeny computation time 1.5125019550323486 s
Do F1 and F2_dual have the same codomain?
True
Codomain matching verification time 0.00016307830810546875 s
Is the point evaluation correct ?
True
Time image verification 0.12107706069946289 s

End verification. Total verification time 1.8280739784240723 s


All tests passed.

Average verification time 1.828143835067749 s
```

```
sage Verify.py -lvl=1 -bench -n=10
###############################################
# Benchmarking level 1 SQIsignHD verification #
# p = 5*2**248 - 1                            #
###############################################

Verifying response number 0
Verifying response number 1
Verifying response number 2
Verifying response number 3
Verifying response number 4
Verifying response number 5
Verifying response number 6
Verifying response number 7
Verifying response number 8
Verifying response number 9

All tests passed.

Average verification time 1.7706276178359985 s
```

### Implementation restriction

The implemented `sagemath` verification algorithm does not recompute the challenge as a hash of the public key, commitment and input message as it should. Indeed, for unexplained reasons, the `SHAKE256` hash function from common python libraries (e.g. `Cryptodome`) do not match with the one from the `Signature/src/common/fips202.c` file used in the signature. For that reason, the challenge is appended to the signature to make the verification work. This should be solved when a full `C` implementation of the verification is available. 

## License 

SQIsignHD is licensed under Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE). Third party code is used in each sublibrary.

The `Signature` sublibrary is deeply based on the SQISign2D-West implementation <https://github.com/SQISign/the-sqisign>:

Apache-2.0: Copyright (c) 2024, The Isogeny Gringos  
(Andrea Basso, Pierrick Dartois, Luca De Feo, Antonin Leroux,
Luciano Maino, Giacomo Pope, Damien Robert and Benjamin Wesolowski)

Within the `Signature` sublibrary, third party code is used in some test and common code files:

- `Signature/src/common/aes_c.c`; MIT: "Copyright (c) 2016 Thomas Pornin <pornin@bolet.org>"
- `Signature/src/common/fips202.c`: Public Domain
- `Signature/src/common/randombytes_system.c`: MIT: Copyright (c) 2017 Daan Sprenkels <hello@dsprenkels.com>

The `Verification` sublibrary imports functions from the `Theta_dim4` submodule <https://github.com/Pierrick-Dartois/Theta_dim4>:

Apache-2.0: Copyright (c) 2024, Pierrick Dartois

Within the `Theta_dim4` submodule, third party code is also used

- `montgomery_isogenies`; MIT: "Copyright (c) 2023 Giacomo Pope"
- `utilities`; "Copyright (c) 2023 Andrea Basso, Luciano Maino and Giacomo Pope"
- `isogenies_dim2`; MIT: "Copyright (c) 2023 Pierrick Dartois, Luciano Maino, Giacomo Pope and Damien Robert"

## References

[1] Pierrick Dartois. Fast computation of higher dimensional isogenies for cryptographic applications. PhD thesis, University of Bordeaux, 2025. In preparation.

[2] Pierrick Dartois, Antonin Leroux, Damien Robert and Benjamin Wesolowski, SQISignHD: New Dimensions in Cryptography, In Advances in Cryptology – EUROCRYPT 2024. https://eprint.iacr.org/2023/436

[3] Andrea Basso, Pierrick Dartois, Luca De Feo, Antonin Leroux, Luciano Maino, Giacomo Pope, Damien Robert and Benjamin Wesolowski, SQIsign2D-West: The Fast, the Small, and the Safer, In Advances in Cryptology – ASIACRYPT 2024. https://eprint.iacr.org/2024/760

[4] Marius A. Aardal, Gora Adj, Diego F. Aranha, Andrea Basso, Isaac Andrés Canales Martínez, Jorge Chávez-Saab, Maria Corte-Real Santos, Pierrick Dartois, Luca De Feo, Max Duparc, Jonathan Komada Eriksen, Tako Boris Fouotsa, Décio Luiz Gazzoni Filho, Basil Hess, David Kohel, Antonin Leroux, Patrick Longa, Luciano Maino, Michael Meyer, Kohei Nakagawa, Hiroshi Onuki, Lorenz Panny, Sikhar Patranabis, Christophe Petit, Giacomo Pope, Krijn Reijnders, Damien Robert, Francisco Rodríguez Henríquez, Sina Schaeffler and Benjamin Wesolowski. SQIsign: Algorithm specifications and supporting documentation, Version 2.0. https://sqisign.org/spec/sqisign-20250205.pdf

[5] Pierrick Dartois, Fast computation of 2-isogenies in dimension 4 and cryptographic applications, IACR Cryptology ePrint Archive, Paper 2024/1180, 2024. https://eprint.iacr.org/2024/1180



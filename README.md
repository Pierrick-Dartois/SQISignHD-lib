# SQISignHD library

This library contains an implementation of FastSQISignHD (Fast Version of the Short Quaternion Isogeny Signature in Higher Dimension), as described in the SQISignHD paper [1]. This library contains two independent sublibraries:
- `Signature`: a `C` implementation of the signature based on the SQISign implementation <https://github.com/SQISign/the-sqisign>.
- `Verification-sage`: an experimental `sagemath` implementation of dimension 4 $2$-isogeny chains in the Theta model of level $2$ applied to the verification of SQISignHD.

For any further information on those libraries and how to run them, one may refer to the `README.md` files of each sublibrary.  

## License 

SQIsignHD is licensed under Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE). Third party code is used in each sublibrary. Credits for third party code are given in the `License` section of the `README.md` files of each sublibrary.

## References

[1] P. Dartois, A. Leroux, D. Robert, and B. Wesolowski. SQIsignHD: New Dimensions in Cryptography. Cryptology ePrint Archive, Paper 2023/436. 2023. url: https://eprint.iacr.org/2023/436.
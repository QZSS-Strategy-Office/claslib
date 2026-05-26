# CLASLIB (CLAS test library)-P/PS

## Modifications by yoronneko

- adds ssr2osr and ssr2obs from [CLAS Test Library](https://qzss.go.jp/en/technical/dod/clas/clas_test-library.html)
- [rearranges direcrory](https://s-taka.org/en/first-trial-of-claslib/#investigation-of-source-code-derived), [organizes compile options](https://s-taka.org/en/first-trial-of-claslib/#organize-compile-options)

## Overview

CLAS test library is a reference implementation of CLAS, and it contains
conversion utilities or a post-processing PPP-RTK positioning utility depending
on its version. Please note it doesn't assure anything in terms of its 
effectiveness or reliability.

CLAS is the Centimeter-Level Augmentation Service of QZSS.
The service is provided as CLAS messages (Ref. Note 1), the format is
defined in RTCM STANDARD 10403.2(MT4073) and called 'Compact SSR'(Ref. Note 2).

The library is derived from RTKLIB (version 2.4.2p13, Ref. Note 3) and
GSILIB (version.1.0.3, Ref. Note 4).

## Document files

### Top level

- [readme_claslib-0.7.3-PS.txt](readme_claslib-0.7.3-PS.txt)
- [readme_rtklib.txt](readme_rtklib.txt)

### doc directory
- [claslib-sourcecode-licence.pdf](docs/claslib-sourcecode-licence.pdf)
- [manual-claslib.pdf](docs/manual-claslib.pdf)
- [relnotes/](docs/relnotes/)


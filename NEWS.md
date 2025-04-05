

# rBiasCorrection NEWS

## Unreleased (2025-04-05)

#### New features

-   working on implementing better guess for nls gauss-newton method
    ([522bce3](https://github.com/kapsner/rBiasCorrection/tree/522bce3b4d47b649929c91d4902d9b7b1117f8fc))
-   added implementation of minpack.lm nls solving with good guess of
    start values
    ([2082832](https://github.com/kapsner/rBiasCorrection/tree/2082832f025ee71311aba105e4fcd9351243c0e9))
-   preparing implementation of alternative algorithm using minpack.lm
    ([426f0b9](https://github.com/kapsner/rBiasCorrection/tree/426f0b95881e24bdabdcbd080aaabc0192df076f))

#### Bug fixes

-   fix issue with deprecated aes_string from ggplot2
    ([1a5d4f0](https://github.com/kapsner/rBiasCorrection/tree/1a5d4f01b17f99fb40b2e68c35d4ad16b9f393ae))
-   removing snapshot-value tests due to cran issues
    ([94084af](https://github.com/kapsner/rBiasCorrection/tree/94084afc0d189326cd8a408d767be8f08833be6e))
-   added seed before call to minpack.lm
    ([a7a9319](https://github.com/kapsner/rBiasCorrection/tree/a7a93193e08e6948c11e047d9342135443dfbf40))
-   fixed implementation of nls_solver
    ([d35ab8b](https://github.com/kapsner/rBiasCorrection/tree/d35ab8bcb3b82c097dc9abb528e0a15299e7246a))

#### Refactorings

-   enhanced code, removed redundancies
    ([c3466e4](https://github.com/kapsner/rBiasCorrection/tree/c3466e4e1c150c22cfac31dd55c0eaae686621ab))
-   cubic-minmax also using nls_solver
    ([b3a5a93](https://github.com/kapsner/rBiasCorrection/tree/b3a5a93e8a82d1c32ebc7ee6c32b63bc0a3f9385))
-   added zzz.r and option to select nls algo
    ([45dd674](https://github.com/kapsner/rBiasCorrection/tree/45dd674ee2064b3e694aaadda01c1d751f57d634))
-   common nls2-code to nls_solver fun
    ([2b2f862](https://github.com/kapsner/rBiasCorrection/tree/2b2f8620d4047a3658fd6303f98841485df6871b))
-   stopifnot with understandable msgs
    ([4fa38fe](https://github.com/kapsner/rBiasCorrection/tree/4fa38fe3e9ad38b9bb2a4ae8111f59525a929efb))

#### Tests

-   added snapshots again
    ([d00ec13](https://github.com/kapsner/rBiasCorrection/tree/d00ec13d6ab0d3ba83505d2050bfa3faaf75eaf0))
-   enhanced unit tests
    ([ab66b88](https://github.com/kapsner/rBiasCorrection/tree/ab66b884e14120d1e6833612b7a244d2e77e39c6))
-   added unit tests for new computation options
    ([c482deb](https://github.com/kapsner/rBiasCorrection/tree/c482deb40047330e1a641180ad982e0d3f5b9c00))

#### CI

-   updated gha
    ([1e2a366](https://github.com/kapsner/rBiasCorrection/tree/1e2a366a2ec564ebe92303656bd501ad19ff550d))
-   updated lint-stage
    ([7740e79](https://github.com/kapsner/rBiasCorrection/tree/7740e79b75167dfb4f9d2d81e3113eaaada516a9))
-   updated gha
    ([1504faa](https://github.com/kapsner/rBiasCorrection/tree/1504faa78692a6ff3450748dcbe78d070c90cde3))

#### Docs

-   beginning with documentation of new nls options
    ([cc6f7e5](https://github.com/kapsner/rBiasCorrection/tree/cc6f7e5584f5ef95f9cddea7f500338c5a924079))

#### Other changes

-   automated gen of readme
    ([01b1dac](https://github.com/kapsner/rBiasCorrection/tree/01b1dac2517bb118f32b31bd8fd55c9696a06f0c))
-   fixed cran errors url
    ([3873d8f](https://github.com/kapsner/rBiasCorrection/tree/3873d8f363859990468eee41b94a192d28fda92f))
-   updated cran-comments
    ([3674986](https://github.com/kapsner/rBiasCorrection/tree/3674986988fb2d24976a7a90fdb96f10b65f3d63))

Full set of changes:
[`v0.3.4...d00ec13`](https://github.com/kapsner/rBiasCorrection/compare/v0.3.4...d00ec13)

## v0.3.4 (2022-06-20)

#### Refactorings

-   removing deprecated future-multiprocess
    ([c8e8ea1](https://github.com/kapsner/rBiasCorrection/tree/c8e8ea1df3d25c254e90972a0e9664d0b505ad84))

#### Tests

-   fixed test errors (linting errors)
    ([c4fe62f](https://github.com/kapsner/rBiasCorrection/tree/c4fe62f83f4e4a172e0a04654622bd50c3f2e925))

#### Docs

-   updated vignette
    ([0a77950](https://github.com/kapsner/rBiasCorrection/tree/0a779507cc7efebb97cf24ac3d651171ca5e09a4))
-   updated url to cran hosted vignette
    ([2d4c2f3](https://github.com/kapsner/rBiasCorrection/tree/2d4c2f3034a35220bfabc397127cf26d5f17b7f9))

#### Other changes

-   updated news.md
    ([08f4f48](https://github.com/kapsner/rBiasCorrection/tree/08f4f48cf1e58d2a2c56cb9e3a12d699d0ba58ef))
-   added dependencies badge
    ([faa1eb7](https://github.com/kapsner/rBiasCorrection/tree/faa1eb75ebf0227c89383753ad96c15967c4e2fb))
-   updated lintr
    ([5fde604](https://github.com/kapsner/rBiasCorrection/tree/5fde604021ae9dc5f083e6133672398f8b8bae91))
-   added news.md
    ([d843819](https://github.com/kapsner/rBiasCorrection/tree/d84381935bd9e06c9d6f74827d047523c4777d57))
-   updating package code, badges, gha
    ([cfac06c](https://github.com/kapsner/rBiasCorrection/tree/cfac06c04e58ff91c09f81066dc4f02aaf288015))

Full set of changes:
[`v0.3.3...v0.3.4`](https://github.com/kapsner/rBiasCorrection/compare/v0.3.3...v0.3.4)

## v0.3.3 (2022-02-16)

Full set of changes:
[`v0.3.2...v0.3.3`](https://github.com/kapsner/rBiasCorrection/compare/v0.3.2...v0.3.3)

## v0.3.2 (2021-08-03)

#### Bug fixes

-   unit tests; update to 0.3.2
    ([715a46d](https://github.com/kapsner/rBiasCorrection/tree/715a46d9f6517a1ca465fad1aa4b2a52bb1fef9d))

Full set of changes:
[`v0.3.1...v0.3.2`](https://github.com/kapsner/rBiasCorrection/compare/v0.3.1...v0.3.2)

## v0.3.1 (2021-06-21)

Full set of changes:
[`v0.3.0...v0.3.1`](https://github.com/kapsner/rBiasCorrection/compare/v0.3.0...v0.3.1)

## v0.3.0 (2021-05-17)

## v0.2.9 (2021-03-29)

Full set of changes:
[`v0.2.8...v0.2.9`](https://github.com/kapsner/rBiasCorrection/compare/v0.2.8...v0.2.9)

## v0.2.8 (2021-03-27)

Full set of changes:
[`v0.2.7...v0.2.8`](https://github.com/kapsner/rBiasCorrection/compare/v0.2.7...v0.2.8)

## v0.2.7 (2021-03-01)

Full set of changes:
[`v0.2.6...v0.2.7`](https://github.com/kapsner/rBiasCorrection/compare/v0.2.6...v0.2.7)

## v0.2.6 (2021-02-13)

Full set of changes:
[`v0.2.5...v0.2.6`](https://github.com/kapsner/rBiasCorrection/compare/v0.2.5...v0.2.6)

## v0.2.5 (2020-12-17)

Full set of changes:
[`v0.2.4...v0.2.5`](https://github.com/kapsner/rBiasCorrection/compare/v0.2.4...v0.2.5)

## v0.2.4 (2020-11-16)

Full set of changes:
[`v0.2.3...v0.2.4`](https://github.com/kapsner/rBiasCorrection/compare/v0.2.3...v0.2.4)

## v0.2.3 (2020-09-22)

Full set of changes:
[`v0.2.2...v0.2.3`](https://github.com/kapsner/rBiasCorrection/compare/v0.2.2...v0.2.3)

## v0.2.2 (2020-09-13)

Full set of changes:
[`v0.2.1...v0.2.2`](https://github.com/kapsner/rBiasCorrection/compare/v0.2.1...v0.2.2)

## v0.2.1 (2020-07-22)

Full set of changes:
[`v0.2.0...v0.2.1`](https://github.com/kapsner/rBiasCorrection/compare/v0.2.0...v0.2.1)

## v0.2.0 (2020-07-17)

Full set of changes:
[`v0.1.7...v0.2.0`](https://github.com/kapsner/rBiasCorrection/compare/v0.1.7...v0.2.0)

## v0.1.7 (2020-06-17)

Full set of changes:
[`v0.1.6...v0.1.7`](https://github.com/kapsner/rBiasCorrection/compare/v0.1.6...v0.1.7)

## v0.1.6 (2020-01-18)

Full set of changes:
[`v0.1.5...v0.1.6`](https://github.com/kapsner/rBiasCorrection/compare/v0.1.5...v0.1.6)

## v0.1.5 (2019-12-16)

Full set of changes:
[`v0.1.4...v0.1.5`](https://github.com/kapsner/rBiasCorrection/compare/v0.1.4...v0.1.5)

## v0.1.4 (2019-11-27)

Full set of changes:
[`v0.1.3...v0.1.4`](https://github.com/kapsner/rBiasCorrection/compare/v0.1.3...v0.1.4)

## v0.1.3 (2019-11-09)

Full set of changes:
[`v0.1.0.9001...v0.1.3`](https://github.com/kapsner/rBiasCorrection/compare/v0.1.0.9001...v0.1.3)

## v0.1.0.9001 (2019-07-28)

Full set of changes:
[`v0.0.1.9002...v0.1.0.9001`](https://github.com/kapsner/rBiasCorrection/compare/v0.0.1.9002...v0.1.0.9001)

## v0.0.1.9002 (2019-06-24)

Full set of changes:
[`f80b0c1...v0.0.1.9002`](https://github.com/kapsner/rBiasCorrection/compare/f80b0c1...v0.0.1.9002)

# Coding Theory

## [BCH code](https://en.wikipedia.org/wiki/BCH_code)
```
.
├── README.md
├── bch.py
├── finitefield.py
├── table_probability.py
├── awgn.py
└── report
    ├── tasks
    └── reports
```

This repository contains software implementation of the Bose-Chaudhuri-Hocquenghem (BCH) code systematic encoder with Berlekamp-Massey decoding algorithm over the GF(2ᵐ).

Program can be run on both Python3 and Python2 without any third-party libraries, only standard ones.

- `finitefield.py` – implementation of basic operations and basic polynomial operations in GF(2ᵐ).
- `bch.py` – encoding and decoding. Verbose output can be enabled using global variable `VERBOSE` from this file.
- `table_probability.py` – test that demonstrates the error probability table with `t = (d-1)/2` and `t + 1` errors in received codewords for different BCH-codes. Carried out at 10,000 iterations.
- `awgn.py` – simulation of data transmission over an additive white Gaussian noise channel (AWGN). In `main` there is a sample how to execute it with multiprocessing: one process on each Eb/N0 value performs 10,000 simulations. Flag `PRE_DEFINED` means to use precomputed result. Flag `PLOT` is option for showing the result as plot.

### Recomendations

It's highly recommended to run huge tests and AWGN simulations with `pypy`, but not with common `cpython`, for the sake of speed increasing.
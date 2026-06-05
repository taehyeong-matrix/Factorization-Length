# Structured Matrix Factorization Length
<p align="center"><a href=""><img src='https://img.shields.io/badge/arXiv-Paper-red?logo=arxiv&logoColor=white' alt='arXiv'></a></p>

This repository contains the code of the paper "Structured Matrix Factorization Length".

## Overview

We study structured matrix factorization length and related factorization varieties for several classical matrix structures. The repository includes Python code for the numerical Toeplitz alternating-minimization examples and Macaulay2 code for the algebraic computations appearing in the paper.

## Files

- `python/`
  - `make_toeplitz_als_figures.py`: generates the Toeplitz alternating-minimization figures
- `output/`
  - Generated figures of the paper
- `macaulay2/`
  - `Example 3.14.m2`: verifies that `diag(1,2,3)` is not a product of two Toeplitz matrices
  - `Example 3.16.m2`: verifies that the upper 3-diagonal matrix `U_3` is not a product of two upper bidiagonal matrices
  - `Example 3.17.m2`: verifies that the 5-diagonal matrix `S` is not a product of two tridiagonal matrices
  - `Example 4.14.m2`: computes numerical evidence for `deg mu_2(T_4) >= 74`
  - `Example 4.15.m2`: computes defining equations for `mu_2(Lambda_3)` and `mu_2(Lambda_4)` by elimination

## Usage

### Python figures

1. Install the required Python packages:

```bash
pip install numpy scipy matplotlib
```

2. Run the figure-generation script:

```bash
python python/make_toeplitz_als_figures.py
```

3. The generated figures are written to `output/`.

To choose a different output directory:

```bash
python python/make_toeplitz_als_figures.py --output-dir path/to/output
```

### Macaulay2 computations

Run each Macaulay2 script from the `macaulay2/` directory. For example:

```bash
cd macaulay2
M2 --script "Example 3.14.m2"
M2 --script "Example 4.14.(1).m2"
```

The scripts include comments describing the expected outputs.

## Generated Figures

- `output/toeplitz_3x3_residual.png`
  - Residual curve for the $3\times 3$ Toeplitz example from Ye--Lim.
- `output/toeplitz_5x5_r15_residual.png`
  - Residual curve for the $5\times 5$ Ye--Lim Example 2 target with `r=15`.
- `output/toeplitz_5x5_r15_factors.png`
  - Heatmap of the computed 15 Toeplitz factors on a common color scale.

## Default Settings

- `3x3` example:
  - `r=2`
  - `als_sweeps=16`
  - `max_nfev=6000`
  - `tol=1e-14`
- `5x5` example:
  - `r=15`
  - `als_sweeps=12`
  - `max_nfev=12000`
  - `tol=1e-14`

## Environment

- Python: `3.13.13`
  - NumPy: `2.4.6`
  - SciPy: `1.17.1`
  - Matplotlib: `3.10.9`
- Macaulay2: `1.25.11`
  - NumericalAlgebraicGeometry: `1.24`
    
## Sample Results

The table below summarizes one run of `python/make_toeplitz_als_figures.py`.

| Example | Length | Output | Final residual |
|---|---:|---|---:|
| Ye--Lim $3\times 3$ Toeplitz example | 2 | `toeplitz_3x3_residual.png` | `2.088e-14` |
| Ye--Lim Example 2 $5\times 5$ target | 15 | `toeplitz_5x5_r15_residual.png`, `toeplitz_5x5_r15_factors.png` | `6.906e-14` |

For faster checks, reduce `als_sweeps` or `max_nfev` inside `python/make_toeplitz_als_figures.py`.

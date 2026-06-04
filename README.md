# Structured Matrix Factorization Length: Toeplitz ALS Figures
<p align="center"><a href=""><img src='https://img.shields.io/badge/arXiv-Paper-red?logo=arxiv&logoColor=white' alt='arXiv'></a></p>

This repository contains the Python code for generating the Toeplitz alternating-minimization figures in the manuscript "Structured Matrix Factorization Length".

## Overview

We run an alternating least-squares method for products of Toeplitz matrices, refine the ALS output by nonlinear least squares, and generate the PNG figures used in the Toeplitz examples of the manuscript.

The script is a cleaned-up standalone version of the local experimental code used during drafting.

## Files

- `make_toeplitz_als_figures.py`
  - Implements the Toeplitz parameterization `vec(T)=Sb`
  - Runs blockwise ALS updates for Toeplitz factors
  - Applies nonlinear least-squares refinement
  - Writes residual plots and a factor heatmap as PNG files

## Usage

1. Install the required Python packages:

```bash
pip install numpy scipy matplotlib
```

2. Run the figure-generation script:

```bash
python make_toeplitz_als_figures.py
```

3. The generated figures are written to `output/`.

To choose a different output directory:

```bash
python make_toeplitz_als_figures.py --output-dir path/to/output
```

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

## Sample Results

The table below summarizes one run of `make_toeplitz_als_figures.py`.

| Example | Length | Output | Final residual |
|---|---:|---|---:|
| Ye--Lim $3\times 3$ Toeplitz example | 2 | `toeplitz_3x3_residual.png` | `2.088e-14` |
| Ye--Lim Example 2 $5\times 5$ target | 15 | `toeplitz_5x5_r15_residual.png`, `toeplitz_5x5_r15_factors.png` | `6.906e-14` |

For faster checks, reduce `als_sweeps` or `max_nfev` inside `make_toeplitz_als_figures.py`.

# CORAL-MOFS

This repository provides the core MATLAB implementation of **CORAL**, a rank-memory search framework for multi-objective feature selection.

## Contents

- `CORAL.m`: main MATLAB implementation of the proposed CORAL algorithm.
- `parameter_settings.txt`: parameter settings used in the manuscript.
- `README.md`: basic usage instructions.

## Requirements

- MATLAB
- PlatEMO platform

## How to use

1. Install PlatEMO.
2. Copy `CORAL.m` into the PlatEMO algorithm directory or the corresponding algorithm folder.
3. Run CORAL with the multi-objective feature selection problem setting used in the manuscript.

## Main parameters

The main CORAL parameters used in the manuscript are:

- `B = 5`: number of sparsity layers in stratified environmental selection.
- `Nm = 3`: number of elite individuals used for local refinement.
- `rho = 0.4`: rank-memory update coefficient.
- `gamma = 0.4`: rank-memory guidance strength.

## Fixed implementation constants

The fixed constants used in the released `CORAL.m` include:

- `F = 0.5`
- `CR = 0.9`
- `sigma = 0.1`
- `Fk = 0.3`
- `Kmin = 1`
- `KmaxFrac = 0.2`

Please see `parameter_settings.txt` for a brief explanation.

## Notes on datasets

The benchmark datasets used in the manuscript are public datasets. Due to file size and redistribution considerations, the datasets are not included in this repository. Users may plug in their own datasets following the PlatEMO problem interface.

## Availability

The code is provided for academic research and reproducibility.
# Wang–Landau Simulation of the Potts Model

This project implements the Wang–Landau algorithm in C++ to estimate the density of states of the q-state Potts model. The resulting density of states is used to compute thermodynamic quantities at a given temperature.

The project was developed for the SI2530 Computational Physics course at KTH Royal Institute of Technology.

## Features

- Simulation of the q-state Potts model on a lattice of arbitrary dimension.
- Wang–Landau sampling to estimate the logarithmic density of states.
- Histogram flatness criterion and configurable convergence tolerance.
- Serial and OpenMP-based parallel simulation methods.
- Computation of thermodynamic quantities:
  - Internal energy
  - Heat capacity
  - Helmholtz free energy
  - Entropy
  - Partition function
  - Energy probability distribution
- Export of histograms and logarithmic densities of states to CSV files.

## Requirements

- C++20-compatible compiler
- CMake >= 4.0
- OpenMP

## Build and Run

From the project directory:

```bash
cmake -S . -B build
cmake --build build
./build/SI2530PROJECT

# SCORE in Lattigo

This repository contains proof-of-concept Go implementations
of bootstrapping algorithms for encrypted real vectors in the CKKS scheme,
based on the [Lattigo library](https://github.com/tuneinsight/lattigo).
It is built around a variant of the SlotToCoeff operation called SCORE.
It focuses on the conventional CKKS bootstrapping [1] and the SPRU bootstrapping [2],
and accompanies the paper [3].

## Table of Contents

- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [References](#references)

## Installation

Simply ensure that you have [Go 1.23](https://go.dev/dl/) or higher installed.

## Repository Structure

Implementations and benchmarks:

- **`boot_tests/main.go`**:
Contains implementations and benchmarks
comparing the conventional bootstrapping procedure [1]
with the SCORE-based variant [3].

- **`spru_tests/main.go`**:
Contains implementations and benchmarks
comparing the SPRU bootstrapping procedure [2]
with the SCORE-based variant [3].

Supporting code:

- **`lattigo-main/`**:
A copy of [Lattigo](https://github.com/tuneinsight/lattigo)
with modifications enabling SCORE and SPRU functionality.

- **`CHANGES.md`**:
A summary of all local modifications made to the Lattigo library.

- **`go.mod` and `go.sum`**:
Standard Go module files managing dependencies and build configuration.

## Usage

The experiments in `boot_tests/main.go` and `spru_tests/main.go`
compare the performance of the conventional CKKS bootstrapping
and the SPRU bootstrapping against their SCORE variants.
All experiments use the parameter sets defined in [3].
The results are printed directly in the terminal as benchmark tables
summarizing timing and precision.

To execute the experiments, run the following from the repository root:
```bash
go run ./boot_tests
```
Or:
```bash
go run ./spru_tests
```

## References

1. Jung Hee Cheon, Kyoohyung Han, Andrey Kim, Miran Kim, and Yongsoo Song.  
   *Bootstrapping for Approximate Homomorphic Encryption*.  
   In *Advances in Cryptology – EUROCRYPT 2018.* Springer, 2018.

2. Jean-Sébastien Coron and Robin Köstler.  
   *Low-Latency Bootstrapping for CKKS using Roots of Unity*, 2025.  
   https://eprint.iacr.org/2025/651.

3. Author(s) anonymous for review.  
   *SCORE: A SlotToCoeff Optimization for Real-Vector Encryption in CKKS*, 2025.
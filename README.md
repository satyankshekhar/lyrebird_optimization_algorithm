# Parallel Lyrebird Optimization Algorithm (LOA)

Sequential, OpenMP, and CUDA implementations of the **Lyrebird Optimization Algorithm** — a bio-inspired metaheuristic that models the escape and hiding behavior of lyrebirds — benchmarked across 10 standard test functions.

## Overview

LOA is a population-based metaheuristic with two search phases:
- **Exploration (escape):** a candidate moves toward a randomly chosen *better* solution in the population, simulating a lyrebird fleeing toward safer ground.
- **Exploitation (hiding):** a candidate makes a small local perturbation around its current position, simulating a lyrebird freezing in place and hiding.

Each individual in the population picks one of the two phases per iteration (roughly 50/50, controlled by a random draw), and keeps the new position only if it improves fitness.

The reference (sequential) version evaluates and updates the population one solution at a time, which scales poorly as population size, dimensionality, or iteration count grows. This repo implements two parallel versions of the same logic:

| Version | Approach |
|---|---|
| **Sequential** | Single-threaded CPU, implemented directly from the original LOA paper |
| **OpenMP** | Multi-core CPU — `#pragma omp parallel for` distributes the per-solution update step across threads |
| **CUDA C** | GPU — each candidate solution is mapped to one CUDA thread, so all N solutions are updated concurrently |

All three versions share the same update rules; only the execution model differs, so results are directly comparable.

## Repository structure

Each notebook covers one benchmark function and contains all three implementations (sequential → OpenMP → CUDA) back to back, compiled and run via Colab's `%%writefile` + `!g++`/`!nvcc` shell magic — the same pattern shown below for Sphere.

```
.
├── Sphere.ipynb
├── Ackley.ipynb
├── Beale.ipynb
├── Booth.ipynb
├── Griewank.ipynb
├── Levy.ipynb
├── Rastrigin.ipynb
├── Rosenbrock.ipynb
├── Schwefel.ipynb
└── Zakharov.ipynb
```

> Rename the list above to match your actual filenames if they differ.

Inside each notebook you'll find three self-contained blocks:
1. `*_serial.cpp` — sequential baseline (`g++`)
2. `*_omp.cpp` — OpenMP version (`g++ -fopenmp`)
3. `*_cuda.cu` — CUDA version (`nvcc`)

Each writes its own source file via `%%writefile`, compiles it, and runs it, printing progress every 1000 iterations, the final best solution, the best fitness value, and total execution time. The CUDA version additionally logs the best fitness and best position at *every* iteration to a CSV file (e.g. `sphere_log.csv`), useful if you want to plot convergence curves.

## Benchmark functions

| Function | Formula | Global Minimum | Domain |
|---|---|---|---|
| Sphere | f(x,y) = x² + y² | f(0,0) = 0 | x,y ∈ [-5.12, 5.12] |
| Rastrigin | 20 + x² + y² − 10(cos 2πx + cos 2πy) | f(0,...,0) = 0 | x,y ∈ [-5.12, 5.12] |
| Rosenbrock | (1−x)² + 100(y−x²)² | f(1,1) = 0 | unbounded |
| Schwefel | 418.9829·2 − [x·sin(√\|x\|) + y·sin(√\|y\|)] | f(420.96, 420.96) = 0 | x,y ∈ [-500, 500] |
| Zakharov | x² + y² + (0.5x+1.5y)² + (0.5x+1.5y)⁴ | f(0,0) = 0 | x,y ∈ [-5, 10] |
| Booth | (x+2y−7)² + (2x+y−5)² | f(1,3) = 0 | x,y ∈ [-10, 10] |
| Ackley | −20·exp(−0.2√(0.5(x²+y²))) − exp(0.5(cos 2πx + cos 2πy)) + e + 20 | f(0,0) = 0 | x,y ∈ [-5, 5] |
| Griewank | 1 + (1/4000)Σx² − Πcos(xᵢ/√i) | f(0,...,0) = 0 | unbounded |
| Lévi | sin²(3πx) + (x−1)²(1+sin²(3πy)) + (y−1)²(1+sin²(2πy)) | f(1,1) = 0 | x,y ∈ [-10, 10] |
| Beale | (1.5−x+xy)² + (2.25−x+xy²)² + (2.625−x+xy³)² | f(3, 0.5) = 0 | x,y ∈ [-4.5, 4.5] |

In the actual notebooks the search runs in 5 dimensions (`DIM = 5`) with bounds `[-100, 100]`, not the 2D textbook domains above — the table just identifies which functions are being tested.

## Requirements

- A C++17-capable compiler (`g++`)
- OpenMP support (bundled with most modern GCC builds)
- NVIDIA CUDA Toolkit (developed/tested on **12.2**) and an NVIDIA GPU
- Tested on **Google Colab** with a Tesla T4 GPU runtime — no local setup needed if you run it there

## Running the notebooks

**On Google Colab (recommended, zero setup):**
1. Open a notebook, e.g. `Sphere.ipynb`.
2. `Runtime → Change runtime type → GPU` (needed for the CUDA cells).
3. Run all cells top to bottom. Each implementation writes its source file, compiles, and executes in sequence.

**Locally, if you have `nvcc`/`g++` installed:**
```bash
# Sequential
g++ sphere_serial.cpp -o sphere_serial
./sphere_serial

# OpenMP
g++ -fopenmp sphere_omp.cpp -o sphere_omp
./sphere_omp

# CUDA (adjust -arch for your GPU's compute capability; sm_75 = Turing/T4)
nvcc -arch=sm_75 sphere_cuda.cu -o sphere_cuda
./sphere_cuda
```

Default parameters across all three versions: population size `256`, `8000` iterations, dimension `5`.

## Results

Execution time (seconds), population = 256, iterations = 8000, on a Tesla T4 (Colab):

| Function | Sequential | OpenMP | CUDA |
|---|---|---|---|
| Ackley | 6.830 | 5.407 | 0.852 |
| Beale | 7.583 | 7.514 | 0.754 |
| Booth | 5.804 | 6.376 | 0.697 |
| Griewank | 9.793 | 8.900 | 0.863 |
| Levy | 16.850 | 13.239 | 1.272 |
| Rastrigin | 6.206 | 6.031 | 0.963 |
| Rosenbrock | 17.266 | 13.461 | 1.201 |
| Schwefel | 15.055 | 12.967 | 1.473 |
| Sphere | 7.340 | 6.488 | 0.748 |
| Zakharov | 7.697 | 6.766 | 0.715 |

**Takeaways:**
- CUDA delivers the largest speedup across the board — roughly 8–13x over the sequential baseline — since all 256 candidates are evaluated and updated concurrently on the GPU.
- OpenMP gives a smaller, more modest gain. The per-solution work here is light (low dimension, simple fitness functions), so thread-management overhead eats into a chunk of the theoretical multi-core benefit.
- Solution accuracy (best fitness found) is comparable across all three versions — parallelizing the execution doesn't change the search logic or the quality of the result, only how fast it gets there.

## Reference

Dehghani, M., Bektemyssova, G., Montazeri, Z., Shaikemelev, G., Malik, O. P., & Dhiman, G. (2023). Lyrebird Optimization Algorithm: A New Bio-Inspired Metaheuristic Algorithm for Solving Optimization Problems. *Biomimetics*, 8(6), 507. https://doi.org/10.3390/biomimetics8060507

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

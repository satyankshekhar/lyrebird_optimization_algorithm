# Parallel Implementation of Lyrebird Optimization Algorithm (LOA)

## DECLARATION
This is to certify that the project work entitled "Parallel Implementation of the Lyrebird Optimization Algorithm using OpenMP and CUDA", submitted in partial fulfillment of the requirements for the award of the Degree of Master of Technology (M.Tech) in Computer Science and Engineering, Department of Computer Science and Engineering, Sardar Vallabhbhai National Institute of Technology (SVNIT), Surat, is an authentic and original work carried out under Prof. Anugrah Jain supervision and guidance. To the best of my knowledge, the content of this Project does not form a basis for the award of any previous Degree to anyone else.

## CERTIFICATE OF APPROVAL
This is to certify that the project work entitled "Parallel Implementation of the Lyrebird Optimization Algorithm using OpenMP and CUDA", submitted by Team, in partial fulfillment of the requirements for the award of the Master of Technology (M.Tech) degree in Computer Science and Engineering, Department of Computer Science and Engineering, Sardar Vallabhbhai National Institute of Technology (SVNIT), Surat, is an authentic and original piece of work carried out by the team under my supervision and guidance. It is understood that by this approval, the undersigned do not necessarily endorse any conclusion drawn or opinion expressed therein, but approve the project for the purpose for which it has been submitted.

---

## ABSTRACT
In this work, we present a parallel implementation of the Lyrebird Optimization Algorithm (LOA), a recently developed bio-inspired metaheuristic that models the escape and hiding behavior of lyrebirds in nature. While the original LOA has shown strong performance in solving complex optimization problems, its sequential design limits its efficiency for large-scale tasks. To address this, we develop a parallel version of LOA using OpenMP for multi-core CPUs and CUDA C for GPU acceleration. The exploration (escape) and exploitation (hiding) phases of LOA were redesigned to support data-level and thread-level parallelism without affecting the core search logic.

Experimental results show that the parallel LOA significantly reduces computation time while maintaining and in many cases improving solution accuracy compared to the original sequential LOA. GPU-based parallelization using CUDA provides the highest performance gain, especially for higher-dimensional problems.

**Keywords:** Parallel LOA, metaheuristic optimization, OpenMP, CUDA C, GPU acceleration, benchmark functions, exploration-exploitation.

---

## INTRODUCTION
Optimization problems appear in many areas, including engineering design, data analysis, industry automation, and scientific computation. Traditional deterministic optimization methods struggle when the problem becomes nonlinear, high-dimensional, or complex. Because of these limitations, researchers increasingly rely on metaheuristic algorithms, which use randomization and nature-inspired strategies to search for high-quality solutions.

The Lyrebird Optimization Algorithm (LOA) is inspired by the escape and hiding strategies of lyrebirds in the wild. These behaviors naturally translate into exploration and exploitation. However, its execution time increases significantly when the number of dimensions or benchmark functions grows, making the standard version less suitable for large-scale applications.

---

## DESCRIPTION OF AVAILABLE SOLUTION (SEQUENTIAL LOA)
The currently available solution for the Lyrebird Optimization Algorithm (LOA) is a sequential implementation. In this version, all computational operations such as generating the initial population, evaluating solutions, and updating positions are executed one after another on a single CPU core.

### Computational Bottleneck
The major performance bottleneck appears because fitness evaluation and solution updates are done solution by solution, not in parallel. This results in:
* Long execution time, especially for large populations.
* Slow convergence, since each iteration takes longer.
* Poor scalability, because it cannot use available CPU/GPU hardware.
* Limited real-world usability for high-dimensional tasks.

### Pseudocode of Sequential LOA
1. Input problem information: variables, objective function, and constraints.
2. Set LOA population size (N) and iterations (T).
3. Generate the initial population matrix at random.
4. Evaluate the objective function and determine the best candidate solution.
5. For t = 1 to T:
    * For i = 1 to N:
        * Determine defense strategy (Phase 1 or Phase 2).
        * Update positions based on Exploration or Exploitation logic.
6. Output the best quasi-optimal solution.

7. ---

## PROPOSED PARALLEL SOLUTION
To address the limitations of sequential processing, this work proposes a parallel implementation of LOA using OpenMP for multi-core CPU parallelization and CUDA C for GPU-based acceleration. The core idea is to retain the logic of the original LOA while distributing computational tasks across multiple processing units. Instead of evaluating one solution at a time, many solutions are evaluated simultaneously.

### OpenMP Architecture
OpenMP is a standardized parallel programming model for shared-memory systems. It uses compiler directives (such as `#pragma omp parallel`) to execute tasks in parallel across multi-core processors. It is especially helpful for accelerating algorithms that contain repeated operations, such as optimization methods.

### CUDA C Architecture
CUDA C is a parallel computing platform developed by NVIDIA for executing general-purpose computations on GPUs. Unlike CPUs, GPUs contain thousands of lightweight cores, making them ideal for computation-heavy tasks and scientific simulations.
* **Key Features of CUDA C:**
    * Designed for GPU-based parallel computing.
    * Enables thousands of threads to run in parallel.
    * Provides hierarchical thread organization (grids, blocks, threads).
    * Supports shared, global, and constant memory for optimized performance.

---

## DESCRIPTION OF IMPLEMENTED SOLUTION
Two parallel execution models were implemented to overcome the sequential bottleneck:

1. **OpenMP-Based Parallel LOA:** Parallelizes major computational parts using multi-core CPU threads. Population fitness evaluation and position updates are distributed using `#pragma omp parallel for`.
2. **CUDA C-Based Parallel LOA:** Offloads expensive operations to the GPU. Each candidate solution is mapped to thousands of CUDA threads inside a kernel function. Memory transfers between CPU and GPU are optimized to reduce overhead.

### Parallel Pseudocode
1. Read problem definition (variables, bounds, objective function).
2. Initialize and generate initial population X randomly.
3. Evaluate fitness of all N solutions.
4. Identify global best solution.
5. For t = 1 to T:
    * # CPU/GPU parallel region depending on mode.
    * Parallel For (i = 1 to N): 
        * Generate random number $r_p$.
        * If $r_p \le 0.5$: Exploration Phase (Identify safe areas $SA_i$, compute new position $XP_{1,i}$).
        * Else: Exploitation Phase (Compute hidden move $XP_{2,i}$).
        * Update $X_i$ if fitness is improved.
    * End Parallel For.
6. Update global best solution and end.

---

## BENCHMARK FUNCTIONS
Benchmark functions are standard mathematical test problems used to evaluate optimization algorithms. They test specific capabilities: unimodal functions measure local search, while multimodal functions test global search ability.

| Function | Formula | Global Minimum | Search Domain |
| :--- | :--- | :--- | :--- |
| **Rastrigin** | $20+x^2+y^2-10(\cos(2\pi x)+\cos(2\pi y))$ | $f(0,...,0)=0$ | $x,y \in [-5.12, 5.12]$  |
| **Rosenbrock** | $f(x,y)=(1-x)^2+100(y-x^2)^2$ | $f(1,1)=0$ | $\infty \le x_i \le \infty$  |
| **Schwefel** | $418.9829 \times 2 - [x \cdot \sin(\sqrt{|x|}) + y \cdot \sin(\sqrt{|y|})]$ | $f(420.96, 420.96)=0$ | $x,y \in [-500, 500]$  |
| **Zakharov** | $f(x,y)=x^2+y^2+(0.5x+1.5y)^2+(0.5x+1.5y)^4$ | $f(0,0)=0$ | $x,y \in [-5, 10]$  |
| **Sphere** | $f(x,y)=x^2+y^2$ | $f(0,0)=0$ | $x,y \in [-5.12, 5.12]$  |
| **Booth** | $f(x,y)=(x+2y-7)^2+(2x+y-5)^2$ | $f(1,3)=0$ | $x,y \in [-10, 10]$  |
| **Ackley** | $-20\exp(-0.2\sqrt{0.5(x^2+y^2)}) - \exp(0.5(\cos(2\pi x)+\cos(2\pi y))) + e + 20$ | $f(0,0)=0$ | $x,y \in [-5, 5]$  |
| **Griewank** | $f(x)=1+(1/4000)\sum(x^2)-\prod \cos(x_i/\sqrt{i})$ | $f(0,...,0)=0$ | $-\infty \le x_i \le \infty$  |
| **Lévi** | $f(x,y)=\sin^2(3\pi x)+(x-1)^2(1+\sin^2(3\pi y))+(y-1)^2(1+\sin^2(2\pi y))$ | $f(1,1)=0$ | $x,y \in [-10, 10]$  |


---

## RESULTS AND DISCUSSION
The CUDA and OpenMP implementations of the Parallel Lyrebird Optimization Algorithm (LOA) were tested using different population sizes, dimensions, and benchmark functions to evaluate performance in terms of execution time and solution accuracy. The results were compared with the original sequential LOA to measure speedup and efficiency. These experiments validate how effectively the parallel approach handles high-dimensional optimization tasks.

### Hardware and Software Specifications
* **Platform:** Google Colaboratory 
* **Runtime:** Python 3 
* **GPU:** NVIDIA Tesla T4 
* **GPU Memory:** 16 GB (T4) 
* **CUDA Version:** 12.2 
* **Parallel Libraries:** OpenMP (enabled through GCC) and CUDA C kernels 

### Performance Comparison Table
The table below shows the execution time (in seconds) for a fixed population size of 256 and 8000 iterations across various benchmark functions.

| Benchmark Function | Sequential LOA (sec) | OpenMP LOA (sec) | CUDA C LOA (sec) |
| :--- | :--- | :--- | :--- |
| **Ackley** | 6.83021 | 5.40686 | 0.851905 |
| **Beale** | 7.58307 | 7.51366 | 0.754362 |
| **Booth** | 5.80362 | 6.37571 | 0.697420 |
| **Griewank** | 9.79295 | 8.89989 | 0.862677 |
| **Levy** | 16.8502 | 13.2389 | 1.272140 |
| **Rastrigin** | 6.20622 | 6.03135 | 0.963288 |
| **Rosenbrock** | 17.2656 | 13.4611 | 1.201019 |
| **Schwefel** | 15.0554 | 12.9668 | 1.472910 |
| **Sphere** | 7.34009 | 6.48816 | 0.747855 |
| **Zakharov** | 7.69702 | 6.76557 | 0.715493 |

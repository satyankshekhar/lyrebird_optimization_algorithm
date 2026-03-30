# Parallel Implementation of Lyrebird Optimization Algorithm (LOA)

## DECLARATION
[cite_start]This is to certify that the project work entitled "Parallel Implementation of the Lyrebird Optimization Algorithm using OpenMP and CUDA", submitted in partial fulfillment of the requirements for the award of the Degree of Master of Technology (M.Tech) in Computer Science and Engineering, Department of Computer Science and Engineering, Sardar Vallabhbhai National Institute of Technology (SVNIT), Surat, is an authentic and original work carried out under Prof. Anugrah Jain supervision and guidance[cite: 23]. [cite_start]To the best of my knowledge, the content of this Project does not form a basis for the award of any previous Degree to anyone else[cite: 24].

## CERTIFICATE OF APPROVAL
[cite_start]This is to certify that the project work entitled "Parallel Implementation of the Lyrebird Optimization Algorithm using OpenMP and CUDA", submitted by Team, in partial fulfillment of the requirements for the award of the Master of Technology (M.Tech) degree in Computer Science and Engineering, Department of Computer Science and Engineering, Sardar Vallabhbhai National Institute of Technology (SVNIT), Surat, is an authentic and original piece of work carried out by the team under my supervision and guidance[cite: 33]. [cite_start]It is understood that by this approval, the undersigned do not necessarily endorse any conclusion drawn or opinion expressed therein, but approve the project for the purpose for which it has been submitted[cite: 34].

---

## ABSTRACT
[cite_start]In this work, we present a parallel implementation of the Lyrebird Optimization Algorithm (LOA), a recently developed bio-inspired metaheuristic that models the escape and hiding behavior of lyrebirds in nature[cite: 52]. [cite_start]While the original LOA has shown strong performance in solving complex optimization problems, its sequential design limits its efficiency for large-scale tasks[cite: 53]. [cite_start]To address this, we develop a parallel version of LOA using OpenMP for multi-core CPUs and CUDA C for GPU acceleration[cite: 54]. [cite_start]The exploration (escape) and exploitation (hiding) phases of LOA were redesigned to support data-level and thread-level parallelism without affecting the core search logic[cite: 55].

[cite_start]Experimental results show that the parallel LOA significantly reduces computation time while maintaining and in many cases improving solution accuracy compared to the original sequential LOA[cite: 57]. [cite_start]GPU-based parallelization using CUDA provides the highest performance gain, especially for higher-dimensional problems[cite: 58].

[cite_start]**Keywords:** Parallel LOA, metaheuristic optimization, OpenMP, CUDA C, GPU acceleration, benchmark functions, exploration-exploitation[cite: 60].

---

## INTRODUCTION
[cite_start]Optimization problems appear in many areas, including engineering design, data analysis, industry automation, and scientific computation[cite: 62]. [cite_start]Traditional deterministic optimization methods struggle when the problem becomes nonlinear, high-dimensional, or complex[cite: 64]. [cite_start]Because of these limitations, researchers increasingly rely on metaheuristic algorithms, which use randomization and nature-inspired strategies to search for high-quality solutions[cite: 65].

[cite_start]The Lyrebird Optimization Algorithm (LOA) is inspired by the escape and hiding strategies of lyrebirds in the wild[cite: 69]. [cite_start]These behaviors naturally translate into exploration and exploitation[cite: 70]. [cite_start]However, its execution time increases significantly when the number of dimensions or benchmark functions grows, making the standard version less suitable for large-scale applications[cite: 71, 72].

---

## DESCRIPTION OF AVAILABLE SOLUTION (SEQUENTIAL LOA)
[cite_start]The currently available solution for the Lyrebird Optimization Algorithm (LOA) is a sequential implementation[cite: 80]. [cite_start]In this version, all computational operations such as generating the initial population, evaluating solutions, and updating positions are executed one after another on a single CPU core[cite: 81].

### Computational Bottleneck
[cite_start]The major performance bottleneck appears because fitness evaluation and solution updates are done solution by solution, not in parallel[cite: 128]. This results in:
* [cite_start]Long execution time, especially for large populations[cite: 132].
* [cite_start]Slow convergence, since each iteration takes longer[cite: 133].
* [cite_start]Poor scalability, because it cannot use available CPU/GPU hardware[cite: 134].
* [cite_start]Limited real-world usability for high-dimensional tasks[cite: 135].

### Pseudocode of Sequential LOA
1. [cite_start]Input problem information: variables, objective function, and constraints[cite: 93].
2. [cite_start]Set LOA population size (N) and iterations (T)[cite: 93].
3. [cite_start]Generate the initial population matrix at random[cite: 93].
4. [cite_start]Evaluate the objective function and determine the best candidate solution[cite: 93].
5. For t = 1 to T:
    * For i = 1 to N:
        * [cite_start]Determine defense strategy (Phase 1 or Phase 2)[cite: 93].
        * [cite_start]Update positions based on Exploration or Exploitation logic[cite: 93].
6. [cite_start]Output the best quasi-optimal solution[cite: 93].

7. ---

## PROPOSED PARALLEL SOLUTION
[cite_start]To address the limitations of sequential processing, this work proposes a parallel implementation of LOA using OpenMP for multi-core CPU parallelization and CUDA C for GPU-based acceleration[cite: 144, 145, 146]. [cite_start]The core idea is to retain the logic of the original LOA while distributing computational tasks across multiple processing units[cite: 147]. [cite_start]Instead of evaluating one solution at a time, many solutions are evaluated simultaneously[cite: 148].

### OpenMP Architecture
[cite_start]OpenMP is a standardized parallel programming model for shared-memory systems[cite: 154]. [cite_start]It uses compiler directives (such as `#pragma omp parallel`) to execute tasks in parallel across multi-core processors[cite: 155, 156]. [cite_start]It is especially helpful for accelerating algorithms that contain repeated operations, such as optimization methods[cite: 157].

### CUDA C Architecture
[cite_start]CUDA C is a parallel computing platform developed by NVIDIA for executing general-purpose computations on GPUs[cite: 173]. [cite_start]Unlike CPUs, GPUs contain thousands of lightweight cores, making them ideal for computation-heavy tasks and scientific simulations[cite: 174, 175].
* **Key Features of CUDA C:**
    * [cite_start]Designed for GPU-based parallel computing[cite: 178].
    * [cite_start]Enables thousands of threads to run in parallel[cite: 179].
    * [cite_start]Provides hierarchical thread organization (grids, blocks, threads)[cite: 180].
    * [cite_start]Supports shared, global, and constant memory for optimized performance[cite: 181].

---

## DESCRIPTION OF IMPLEMENTED SOLUTION
[cite_start]Two parallel execution models were implemented to overcome the sequential bottleneck[cite: 207]:

1. [cite_start]**OpenMP-Based Parallel LOA:** Parallelizes major computational parts using multi-core CPU threads[cite: 209]. [cite_start]Population fitness evaluation and position updates are distributed using `#pragma omp parallel for`[cite: 210].
2. [cite_start]**CUDA C-Based Parallel LOA:** Offloads expensive operations to the GPU[cite: 216]. [cite_start]Each candidate solution is mapped to thousands of CUDA threads inside a kernel function[cite: 217]. [cite_start]Memory transfers between CPU and GPU are optimized to reduce overhead[cite: 219].

### Parallel Pseudocode
1. [cite_start]Read problem definition (variables, bounds, objective function)[cite: 236].
2. [cite_start]Initialize and generate initial population X randomly[cite: 237, 238].
3. [cite_start]Evaluate fitness of all N solutions[cite: 240].
4. [cite_start]Identify global best solution[cite: 243].
5. For t = 1 to T:
    * [cite_start]# CPU/GPU parallel region depending on mode[cite: 245].
    * [cite_start]Parallel For (i = 1 to N): [cite: 246]
        * [cite_start]Generate random number $r_p$[cite: 247].
        * [cite_start]If $r_p \le 0.5$: Exploration Phase (Identify safe areas $SA_i$, compute new position $XP_{1,i}$)[cite: 248, 249, 250].
        * [cite_start]Else: Exploitation Phase (Compute hidden move $XP_{2,i}$)[cite: 252, 253].
        * [cite_start]Update $X_i$ if fitness is improved[cite: 254, 255].
    * [cite_start]End Parallel For[cite: 256].
6. [cite_start]Update global best solution and end[cite: 257, 258].

---

## BENCHMARK FUNCTIONS
[cite_start]Benchmark functions are standard mathematical test problems used to evaluate optimization algorithms[cite: 283]. [cite_start]They test specific capabilities: unimodal functions measure local search, while multimodal functions test global search ability[cite: 284, 285].

| Function | Formula | Global Minimum | Search Domain |
| :--- | :--- | :--- | :--- |
| **Rastrigin** | $20+x^2+y^2-10(\cos(2\pi x)+\cos(2\pi y))$ | $f(0,...,0)=0$ | [cite_start]$x,y \in [-5.12, 5.12]$ [cite: 294] |
| **Rosenbrock** | $f(x,y)=(1-x)^2+100(y-x^2)^2$ | $f(1,1)=0$ | [cite_start]$\infty \le x_i \le \infty$ [cite: 294] |
| **Schwefel** | $418.9829 \times 2 - [x \cdot \sin(\sqrt{|x|}) + y \cdot \sin(\sqrt{|y|})]$ | $f(420.96, 420.96)=0$ | [cite_start]$x,y \in [-500, 500]$ [cite: 294] |
| **Zakharov** | $f(x,y)=x^2+y^2+(0.5x+1.5y)^2+(0.5x+1.5y)^4$ | $f(0,0)=0$ | [cite_start]$x,y \in [-5, 10]$ [cite: 294] |
| **Sphere** | $f(x,y)=x^2+y^2$ | $f(0,0)=0$ | [cite_start]$x,y \in [-5.12, 5.12]$ [cite: 294] |
| **Booth** | $f(x,y)=(x+2y-7)^2+(2x+y-5)^2$ | $f(1,3)=0$ | [cite_start]$x,y \in [-10, 10]$ [cite: 297] |
| **Ackley** | $-20\exp(-0.2\sqrt{0.5(x^2+y^2)}) - \exp(0.5(\cos(2\pi x)+\cos(2\pi y))) + e + 20$ | $f(0,0)=0$ | [cite_start]$x,y \in [-5, 5]$ [cite: 297] |
| **Griewank** | $f(x)=1+(1/4000)\sum(x^2)-\prod \cos(x_i/\sqrt{i})$ | $f(0,...,0)=0$ | [cite_start]$-\infty \le x_i \le \infty$ [cite: 297] |
| **Lévi** | $f(x,y)=\sin^2(3\pi x)+(x-1)^2(1+\sin^2(3\pi y))+(y-1)^2(1+\sin^2(2\pi y))$ | $f(1,1)=0$ | [cite_start]$x,y \in [-10, 10]$ [cite: 297] |


---

## RESULTS AND DISCUSSION
[cite_start]The CUDA and OpenMP implementations of the Parallel Lyrebird Optimization Algorithm (LOA) were tested using different population sizes, dimensions, and benchmark functions to evaluate performance in terms of execution time and solution accuracy[cite: 300]. [cite_start]The results were compared with the original sequential LOA to measure speedup and efficiency[cite: 301]. [cite_start]These experiments validate how effectively the parallel approach handles high-dimensional optimization tasks[cite: 302].

### Hardware and Software Specifications
* [cite_start]**Platform:** Google Colaboratory [cite: 304]
* [cite_start]**Runtime:** Python 3 [cite: 305]
* [cite_start]**GPU:** NVIDIA Tesla T4 [cite: 306]
* [cite_start]**GPU Memory:** 16 GB (T4) [cite: 307]
* [cite_start]**CUDA Version:** 12.2 [cite: 308]
* [cite_start]**Parallel Libraries:** OpenMP (enabled through GCC) and CUDA C kernels [cite: 309]

### Performance Comparison Table
[cite_start]The table below shows the execution time (in seconds) for a fixed population size of 256 and 8000 iterations across various benchmark functions[cite: 222, 316].

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

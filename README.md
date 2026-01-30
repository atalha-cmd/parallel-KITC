# parallel-KITC

**Parallel Kernel-Independent Treecode (KITC)** for fast and scalable evaluation of long-range particle interactions.

---

## Overview
This repository presents a **parallel implementation of a kernel-independent barycentric Lagrange treecode (KITC)** designed to efficiently accelerate pairwise interaction computations. The method replaces the direct **O(N²)** summation with a hierarchical, tree-based approximation of **O(N log N)** complexity, enabling scalable simulations for large particle systems.

Unlike classical multipole-based methods, KITC is **kernel-independent** and employs interpolation-based approximations, making it broadly applicable to a wide class of interaction kernels without requiring analytic expansions.

---

## Parallelization
The implementation supports scalable parallel execution through:
- **OpenMP** for shared-memory parallelism  
- **MPI** for distributed-memory parallelism  
- **GPU acceleration via Kokkos** *(in progress)*  

---

## Reference
Wang, L. (2020). *“A Kernel-Independent Treecode Based on Barycentric Lagrange Interpolation.”*  Communications in Computational Physics, **28**, 1415–1436.

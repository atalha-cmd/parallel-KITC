# parallel-KITC

**Parallel Kernel-Independent Treecode (KITC)** for fast and scalable evaluation of long-range particle interactions.

---

## Overview
This repository presents a **Parallel implementation of the Kernel-independent Barycentric Lagrange treecode (KITC)** designed to efficiently accelerate pairwise interaction by replacing the direct **O(N²)** summation with a hierarchical, tree-based approximation of **O(N log N)** complexity.

Unlike classical multipole-based methods, KITC is **kernel-independent** and relies on interpolation-based approximations rather than analytic expansions. This design makes the method flexible and broadly applicable to a wide class of interaction kernels.

---

## Parallelization
The implementation supports scalable parallel execution through:
- **OpenMP** for shared-memory parallelism  
- **MPI** for distributed-memory parallelism  
- **GPU acceleration via Kokkos** *(in progress)*  

---

## Reference
Wang, L. (2020). *“A Kernel-Independent Treecode Based on Barycentric Lagrange Interpolation.”* Communications in Computational Physics, **28**, 1415–1436.

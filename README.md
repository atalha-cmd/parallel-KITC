# parallel-KITC

**Parallel Implementation of Kernel-Independent Treecode (KITC)** for fast and scalable evaluation of long-range particle interactions.

---

## Overview
This repository implements a **Parallel version of a Kernel-Independent Barycentric Lagrange Treecode Algorithm** that accelerates pairwise interaction computations by replacing the direct O(N^2) summation with a hierarchical tree-based approximation of O(N log N) complexity.

Unlike classical multipole methods, KITC is **kernel-independent** and relies on interpolation-based approximations, making it flexible and applicable to a wide range of interaction kernels.

## Parallelization
The code supports scalable parallel execution using:
- **OpenMP** for shared-memory parallelism
- **MPI** for distributed-memory parallelism
- **GPU acceleration via Kokkos** (in progress)

---
## Reference
Wang, Lei. (2020). ”A Kernel-Independent Treecode Based on Barycentric Lagrange Interpolation.” Communications in Comp.Physics. 28. 1415-1436.

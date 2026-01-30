# parallel-KITC

**Parallel Implementation of Kernel-Independent Treecode (KITC)** for fast and scalable evaluation of long-range particle interactions.

---

## Overview
This repository implements a **parallel kernel-independent treecode (KITC)** that accelerates
pairwise interaction computations by replacing the direct \(\mathcal{O}(N^2)\) summation with a hierarchical tree-based approximation of \(\mathcal{O}(N \log N)\) complexity.

Unlike classical multipole methods, KITC is **kernel-independent** and relies on interpolation-based approximations, making it flexible and applicable to a wide range of interaction kernels.

---

## Key Features
- Hierarchical spatial tree construction
- Kernel-independent interpolation-based approximation
- Near-field direct interactions and far-field tree approximation
- Accuracy validation against direct summation
- Designed for large-scale particle systems

---

## Parallelization
The code supports scalable parallel execution using:
- **OpenMP** for shared-memory parallelism
- **MPI** for distributed-memory parallelism
- **GPU acceleration via Kokkos** (in progress)

---


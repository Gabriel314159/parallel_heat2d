# Parallel Two-Dimensional Heat Diffusion with MPI and OpenMP

Status: report outline. Replace prompts with the group's implementation details and measured results. No performance conclusions are available yet.

## 1. Problem and numerical method

Explain the heat equation, assumptions, initial/boundary conditions, grid, five-point stencil and explicit time-step stability condition. Define whether dimensions count interior points.

## 2. Serial implementation

Describe contiguous storage, double buffering, the update loop, parameter validation and timing scope. Give the actual build/run instructions once the source exists.

## 3. Parallel design

### 3.1 Domain decomposition

Explain row and 2D block layouts, process-grid selection, local/global coordinates and remainder handling.

### 3.2 Halo communication

Explain message partners, buffer layout, column handling and physical boundaries. Compare the implemented blocking and overlap designs.

### 3.3 OpenMP computation

Explain thread work distribution, one-writer updates, MPI thread support and the rank/thread layouts tested.

## 4. Correctness validation

Report the actual serial/parallel errors, constant-field checks, physical boundary checks, non-divisible cases and small-block cases. Separately discuss analytical discretization error and refinement.

## 5. Experimental setup

Identify Atlas or Vanda, CPU/node details, compiler and library versions, build flags, rank/thread placement, problem sizes, warmups/repetitions and timing definitions. State the commit used.

## 6. Runtime, scaling and efficiency

Insert measured figures and discuss strong scaling, weak scaling and fixed-core strategy comparisons. State the baseline for every speedup or efficiency plot. Keep raw-data references and uncertainty visible.

## 7. Discussion and limitations

Explain observed bottlenecks using evidence. Discuss communication, local block sizes, memory bandwidth, threading, problem size and the limits of the experiments. Do not claim overlap is faster before measuring it.

## 8. Conclusion

Summarize the demonstrated correctness and measured performance findings. State which strategies worked under the tested conditions and which questions remain open.

## 9. AI assistance

Summarize the factual record in [docs/ai_usage.md](../docs/ai_usage.md), including subsequent assistance and checks.

## References

Include the technical sources actually used. The study uses synthetic conditions; no real-world dataset conclusion is currently applicable.

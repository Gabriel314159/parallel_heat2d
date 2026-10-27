# AI assistance record

The assignment requires an explanation of AI use. Keep this record factual and update it as the project develops.

## Assistance already used

On 20 October 2026, OpenAI ChatGPT assisted with comparing project topics, explaining two-dimensional heat diffusion and parallelization, proposing the repository architecture, and drafting this initial documentation/configuration scaffold.

The generated scaffold contains README files, configuration files, a CMake build skeleton, experiment plans, an HPC shell template and a report outline. It contains no numerical solver implementation, benchmark runner or measured performance results.

At that stage, the team still needed to implement the solvers and run validation. Generated plans were not evidence that experiments had been performed.

## Append future entries

For each substantial use, record:

- Date and tool/model, if known.
- Task or question.
- Files or implementation areas affected.
- What was accepted, edited or rejected.
- Human review and actual checks performed.

Do not invent a model version, completed review or validation result. The final report should summarize the assistance actually used and how the team checked it.

## Implementation assistance on 27 October 2026

OpenAI ChatGPT generated the initial C++ implementation (shared modules, serial driver, MPI decomposition, blocking/nonblocking halo exchange and OpenMP kernels), Python validation/benchmark/plotting tools, VS Code configuration and updated running instructions. It also generated the integration checks and executed local builds, numerical comparisons, analytical refinement and small workflow tests. The observed checks and limitations are recorded in [verification.md](verification.md).

The group should review and understand this code, reproduce the checks on its own WSL environment, and perform the required HPC runs before drawing performance conclusions. Local generated example fields and completed test records are distinct from the pending Atlas/Vanda measurements.

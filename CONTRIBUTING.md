# Working on this project

Keep numerical behavior consistent across implementations. Make one focused change per branch or pull request, and explain its purpose and validation.

## Shared conventions

- C++17, four-space indentation and LF line endings.
- Use `double` for temperature values and numerical parameters.
- Store grid data contiguously and document local/global indexing.
- Keep the serial target independent of MPI and OpenMP.
- Keep MPI calls on the initializing thread in the initial hybrid implementation.
- Avoid changing boundary conditions, time steps or stopping rules while comparing execution strategies.
- Record changes to the planned CLI in both the code and documentation.

## Suggested branch names

- `feature/serial-solver`
- `feature/mpi-rows`
- `feature/openmp-kernel`
- `feature/cartesian-decomposition`
- `feature/overlap-exchange`
- `experiments/hpc-scaling`

## Before merging numerical changes

Run the implemented correctness checks relevant to the change. Record the command, environment and observed result. A matching serial and parallel output does not independently validate the shared numerical formula; retain analytical and invariant-based checks as well.

No solver checks exist in the initial scaffold. Do not mark planned checks as passed until they have actually run.

## Results and documentation

Commit compact raw benchmark CSV files and the selected figures needed to reproduce the report. Keep temporary outputs, large field dumps and build artifacts out of Git. Record exact configurations and commits for reported measurements.

Update [docs/ai_usage.md](docs/ai_usage.md) when AI assistance is used. The team remains responsible for understanding and validating the implementation.

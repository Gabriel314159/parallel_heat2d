# AI assistance record

The assignment requires an explanation of AI use. Keep this record factual and update it as the project develops.

## Assistance already used

On 20 October 2026, OpenAI ChatGPT assisted with comparing project topics, explaining two-dimensional heat diffusion and parallelization, proposing the repository architecture, and drafting this initial documentation/configuration scaffold.

The generated scaffold contains README files, configuration files, a CMake build skeleton, experiment plans, an HPC shell template and a report outline. It contains no numerical solver implementation, benchmark runner or measured performance results.

The team still needs to review the numerical conventions, confirm platform-specific instructions, implement the solvers and run the validation and HPC experiments. Generated plans are not evidence that experiments have been performed.

On 26–27 October 2026, Anthropic Claude (Claude Code, cloud session) explained the existing problem formulation, proposed extending the model from heat diffusion to advection–diffusion with a prescribed divergence-free flow, and drafted the revised `docs/design.md` (cell-centred finite-volume grid, streamfunction face velocities, upwind flux scheme, Dirichlet/no-flux walls, combined time-step rule). A throwaway NumPy prototype was run to check discrete divergence, mass conservation, constant-field invariance and the absence of new extrema. Accepted by the team: _[fill in]_. Human review performed: _[fill in]_.

On 27 October 2026, Anthropic Claude (Claude Code, local session; Claude Sonnet 5.5, later Claude Opus 5.5) assisted as follows:

- Documentation: converted the notation in `docs/design.md` to LaTeX; at the team's direction, restricted the model to the time-periodic double gyre with no-flux walls only (Dirichlet walls removed; the analytic reference changed from sine to cosine decay) and updated the configuration, experiment files and documentation to match.
- C++ tutorial: a team member wrote `src/common/velocity.cpp` and parts of `src/common/config.cpp` with step-by-step guidance and tests written by Claude.
- Implementation by Claude: the remaining shared components (`grid.hpp`, `physics`, `stencil.hpp`, `io`, the `FlowField` fast velocity fill), the serial solver (`src/serial/`), the MPI + OpenMP solver (`src/parallel/`: decomposition, halo exchange, blocking and overlapping modes), the CMake build, the unit tests in `tests/`, `scripts/check_correctness.py`, `scripts/plot_field.py`, `scripts/animate_field.py` and `hpc/big_run.sh`.
- Checks run during the session, on an Apple M-series laptop: all `ctest` suites and the serial correctness checks passed. Relative mass change at most 2e-13 after up to 42,500 steps; constant field kept to 4e-14; cosine decay matched the exact discrete solution to 2e-15; the error against the continuous solution fell by a factor of about 3.95 per grid refinement (second order). Animation frames were inspected visually. With Open MPI 5.0 and OpenMP on the same laptop, every correctness case (row and 2D decompositions, blocking and overlapping exchange, up to 6 ranks × 2 threads, uneven and 3 × 3 grids) reproduced the serial field exactly (maximum difference 0), including coarsened snapshots.
- Not done: the benchmark runner, scaling plots and any HPC measurements. Informal laptop timings are not report results.

Accepted by the team: _[fill in]_. Human review performed: _[fill in]_.

On 28 October 2026, Anthropic Claude (Claude Code, local session; Claude Opus 5.5) continued the work:

- Design: drafted the computed wind-driven flow for `docs/design.md` (vorticity–streamfunction form, Arakawa Jacobian, third-order Adams–Bashforth, free-slip walls, corner storage and the corner-aware halo exchange) and restructured the document into the team's writing style. At the team's direction, the double gyre was dropped as a model, bottom friction was left out, and the β term was removed; the β-free description of the flow was revised accordingly, with a placeholder for the final parameters.
- Implementation by Claude: the computed flow in both solvers (`src/common/flow.*`, `corners.hpp`, `stability.*`), the halo exchange with corners, the fixed time step with a run-time speed check, snapshots of ψ and the vorticity, the updated plotting scripts, configuration, experiments and documentation. The Poisson equation for ψ was first solved directly with a sine transform; at the team's request this was replaced by a matrix-free, warm-started conjugate gradient solver (`src/common/poisson.hpp`), shared by both solvers, with the MPI sums and halo exchanges supplied by the parallel solver.
- Checks run during the session (Apple M-series laptop, Open MPI 5.0): all `ctest` suites, including new ones for the CG solver (exact eigenmode, true residual within the tolerance, warm start, zero right-hand side, second-order convergence, a stop with a message when the iteration limit is reached), the flow (Arakawa energy and enstrophy conservation, rest stays at rest, linear steady state against its series solution, a manufactured solution converging at second order) and the halo exchange on 1, 4 and 6 ranks; `scripts/check_correctness.py` with serial/parallel differences of at most about 5e-14. Measured CG cost: about 180 iterations per step on 128 × 256 cells and 350 on 256 × 512, with the Poisson solve taking over 99% of the run time. Parameter scans on a 128 × 256 grid showed that, without β, the steady wind spins up two smooth gyres over a long viscous time, with no eddies within 100 time units up to Re of about 4000.

Accepted by the team: _[fill in]_. Human review performed: _[fill in]_.

## Append future entries

For each substantial use, record:

- Date and tool/model, if known.
- Task or question.
- Files or implementation areas affected.
- What was accepted, edited or rejected.
- Human review and actual checks performed.

Do not invent a model version, completed review or validation result. The final report should summarize the assistance actually used and how the team checked it.

# Initial scaffold checks

Historical record for the 20 October documentation-only package. The implemented version and its later tests are described in [verification.md](verification.md).

Checked on 20 October 2026 in a Linux development environment. These checks concern repository files and build configuration, not the unimplemented heat solver.

Completed:

- Parsed all five JSON files and resolved 35 planned parallel resource/case combinations.
- Checked the planned numerical parameters against the explicit-method stability condition.
- Checked weak-scaling invariants: fixed spatial spacing, time step and interior cells per CPU core.
- Checked explicit parallel backend selection in the three performance manifests.
- Checked relative Markdown links, UTF-8 text, LF line endings and final newlines.
- Checked that the CSV template contains a header only.
- Checked HPC template syntax with `bash -n`.
- Configured the default CMake scaffold and ran its empty build successfully using CMake 4.4.4.
- Confirmed that enabling solver builds stops with the intended missing-source message.
- Confirmed that scaffold mode produces no solver executables.

Not performed: numerical simulation, solver compilation, MPI/OpenMP runtime checks, benchmark-tool execution, HPC submission or performance measurement. The solver sources and experiment utilities still need to be implemented.

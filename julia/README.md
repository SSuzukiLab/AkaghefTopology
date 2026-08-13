# Topology.jl migration workbench

This directory is an in-progress MATLAB-to-Julia migration. It is not yet a
complete or fidelity-verified port.

The authoritative scope and current evidence are in
[`LIVE_SCRIPT_SCOPE.md`](LIVE_SCRIPT_SCOPE.md) and
[`LIVESCRIPT_PORTS.md`](LIVESCRIPT_PORTS.md). The original five generated SVGs
are rejected prototypes: they omit the MATLAB/Sage o-graph layout and crossing
semantics and must not be treated as valid output.

Current smoke tests can be run with:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Passing these tests only validates the currently implemented subset. Migration
acceptance additionally requires all 17 topology/dependency Live Scripts and all
44 o-graph plot calls to match the MATLAB structural and visual baselines.
The full-resolution manual review of the 46 dynamic Julia plot instances is
recorded in [`VISUAL_VERIFICATION.md`](VISUAL_VERIFICATION.md).

Canonical ownership and the personal/iPad return path are defined in
[`CANONICAL_BOUNDARY.md`](CANONICAL_BOUNDARY.md). In particular, Mac replay
evidence is not Windows Math Harness verification, and personal-surface edits
cannot overwrite either this library or canonical Harness state without an
explicit import, rerun, and review on Windows.

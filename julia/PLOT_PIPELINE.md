# O-graph plot pipeline, parity, and the Windows constraint

Audit dates: 2026-08-20 (initial), 2026-08-20 (revised after two retractions
and the Windows finding).
Evidence status: `mac-replayed` for every measurement below, except section 7,
which is read from Windows-side logs supplied on 2026-08-20.

Environment used for the Mac measurements:

- Julia 1.6.3, `/Applications/Julia-1.6.app`, project `julia/`
- SageMath 10.6, `/Applications/SageMath-10-6.app`
- repository at commit `e941b76`

Harness: `julia/scripts/sage/sage_layout_parity.py` and
`julia/scripts/sage/julia_layout_parity.jl`.

---

## 1. Summary

Three things are established.

**The Julia port runs.** The `VLExample.mlx` equivalent executes all 19 plot
calls without exception and writes 19 SVGs. All 19 satisfy P1 (every segment
axis-parallel) and P5 (crossing positions distinct). Three fail P3: one edge
per diagram ends on the wrong crossing (section 5, IS7).

**The layout algorithm was ported faithfully for diagrams without self-loops.**
Given identical bending numbers, Julia reproduces the fork's corner sequence
exactly on the hopf link and on both C251020 figures. It does **not** agree on
diagrams containing self-loops (section 4).

**The MATLAB plot path cannot run on Windows.** This is architectural, not a
configuration fault (section 7). The consequence that matters most for this
migration is that the canonical Windows `REdgeTable.Position` fixture, which
several documents still list as the pending acceptance gate, cannot be produced
at all under the current design. Those gates are void.

Acceptance criterion, decided 2026-08-20: the target is **not** reproducing
MATLAB's figures. It is that the algorithm delivers the figure quality it
specifies. Mismatches caused by the minimum-bend MILP's degeneracy are
accepted. MATLAB's output is a prior implementation, not an oracle.

---

## 2. Two retractions

Recorded because both claims were circulated before being checked, and the
later sections were rewritten around them.

**R1 — "Given the same bending numbers the Julia layout reproduces Sage
exactly" was over-generalised.** It was measured on hopf, borromean and
whitehead, none of which contain a self-loop. On `s2xs1_calculation`, which
does, the two disagree. Section 4 states the corrected scope.

**R2 — the Sage-side P6 results were produced by an invalid method, and the
"Sage emits non-orthogonal segments" claim is withdrawn.** `allocate_pos.py`
does not emit polylines. It extracts **bezier control points**:

```python
bezier_coords = [im[0][0][0], im[0][0][1], im[-1][0][0], im[-1][0][1]]  # cubic
bezier_coords = [p, im[c+1][0][0], q]                                   # quadratic
```

Control points other than the endpoints do not lie on the curve, so running a
segment-intersection test over control polygons measures nothing about the
drawn figure. Upstream Sage renders these with `bezier_path`; orthogonality of
self-loops was never part of its contract. No upstream Sage defect has been
demonstrated by this work.

What survives R2: **endpoints of a control polygon do lie on the curve**, so an
endpoint check (P3) is valid against the fork's output. That check has not been
run on the Sage side; it is the cheapest remaining comparison (section 8).

---

## 3. As-is debt split

### 3.1 MATLAB path

| Concern | Owner | Location |
|---|---|---|
| Data model, o-data, `REdgeTable`, weights, moves, invariants | MATLAB | `Manifold/VirtualLink/VirtualLink.m` |
| Process bridge | MATLAB | `SW/SageWrapper.m` — `pyenv` loads the Sage venv interpreter **into the MATLAB process**, imports `sage.all`, and runs code with `py.builtins.exec` against `sage.all.__dict__`. Diagram state lives in that namespace under `sageName`, tracked by `sageLinked`. This in-process, stateful design is what fails on Windows (section 7) |
| Layout geometry | Python | `Manifold/VirtualLink/allocate_pos.py` — a **modified fork** of Sage's `Link.plot` internals, injected as a string template by `VirtualLink.calcPos` |
| Topological data, MILP solver | Sage | `Link.pd_code`, `.orientation`, `.regions`, `._isolated_components`, `MixedIntegerLinearProgram` |
| Curve reconstruction | MATLAB | `VirtualLink.calcPos` → `connectPolylines`, which **concatenates the returned coordinate lists as polylines**. `bezier.m` and `bernstein.m` exist in the same folder but are not called from `VirtualLink.m` |

The last row is the orthogonalisation: the fork's output is bezier control
data, and the MATLAB side consumes it as a polyline. For 2- and 3-point entries
that is faithful (a quadratic control triple and a polyline corner have the
same corner sequence). For the 4-point self-loop entry it is not — read as a
polyline it becomes a triangle with diagonal edges.

### 3.2 What the fork changed relative to upstream Sage 10.6

Measured by diffing `allocate_pos.py`'s `allocate_positions` against
`sage/knots/link.py` lines 3630–3860.

| Change | Direction |
|---|---|
| `bendingNumbers` override parameter | added by the fork |
| `solver=solver` on both MILPs | **removed by the fork.** Upstream `plot(self, gap=0.1, component_gap=0.5, solver=None, ...)` (`link.py:3414`) threads a backend through both `MixedIntegerLinearProgram` calls (`link.py:3639`, `:3756`) |
| `MLP.get_values(v, convert=ZZ, tolerance=1e-3)` (`link.py:3681`) | **removed by the fork**, leaving the raw `get_values(v)` |
| `bezier_path(...)` drawing | replaced by control-point extraction |
| `gap` / `tailshort` / `headshort` crossing-gap logic | removed; MATLAB applies its own gap |
| `_isolated_components` driver, component x-offset, unknot circles, abalone reversal | added by the fork, at module level |

Consequence: **IS2 and IS3 below are fork-introduced regressions, not upstream
defects.** Upstream lets the caller pin the backend and converts the solution
to integers. Restoring both on the MATLAB side is a few lines.

### 3.3 Julia path

No Sage, no Python. Dependencies are `AbstractAlgebra`, `GLPK`, `JSON`,
`JuMP`, `LinearAlgebra` — all pure Julia or JLL-packaged, with no system
prerequisite.

| Stage | Location | Parity / status |
|---|---|---|
| Virtual planarisation | `src/virtual_planarization.jl` — `planarize_virtual_gauss` | not compared |
| PD code | `src/virtual_link.jl` — `pd_code` | matches Sage on the measured cases |
| Regions | `src/virtual_link.jl` — `regions` | matches `Link.regions()` on hopf |
| MILP 1: minimum bend | `src/orthogonal_layout.jl` — `minimal_bending_numbers` | degenerate, see IS1 |
| MILP 2: segment lengths | `src/orthogonal_layout.jl` — `_piece_lengths` | matches on self-loop-free diagrams |
| Routing | `src/orthogonal_layout.jl` — `_route_edges` | matches on self-loop-free diagrams; **misroutes one edge on three diagrams, see IS7** |
| Arc aggregation | `src/svg_plot.jl` — `_real_arc_layouts` | not compared |
| SVG emission | `src/svg_plot.jl` — `plot_svg` | not compared |

Not ported: `_isolated_components`, `component_gap`, unknot circles, abalone
reversal. `_route_edges` raises `ArgumentError` on a disconnected PD instead.

---

## 4. Measured parity with the fork

### Part A — MILP 1 is degenerate and backend-dependent

Every row is an optimal solution; the objectives tie.

| Case | Sage / Coin (default) | Sage / GLPK | Sage / PPL | Julia / GLPK | objective |
|---|---|---|---|---|---|
| hopf | `[0,4,2,2]` | `[2,2,4,0]` | `[0,4,2,2]` | `[2,2,4,0]` | 8 |
| trefoil † | `[0,2,0,2,1,3]` | `[1,3,0,2,0,2]` | `[1,3,0,2,0,2]` | `[1,3,0,2,0,2]` | 8 |
| borromean | `[3,0,1,0,3,0,0,1,1,1,0,2]` | `[1,1,0,2,3,0,1,0,3,0,0,1]` | `[2,0,0,2,1,1,1,1,4,0,0,0]` | `[1,1,1,1,2,2,0,0,4,0,0,0]` | 12 |
| whitehead | `[-1,-1,-1,-1,0,-3,0,0,3,0]` | `[-1,-1,0,-2,0,-3,0,1,2,0]` | `[0,-1,-1,-2,-1,-2,0,0,3,0]` | `[-1,-1,-1,-1,0,-3,0,0,3,0]` | 10 |

Backend parity alone does not make the sides agree: Julia/GLPK matches
Sage/GLPK on hopf and trefoil but not on borromean. Under the accepted
criterion this is tolerated, not fixed.

† Measured on the PD matrix `[1 2 4 3; 3 4 6 5; 5 6 2 1]` directly. That matrix
does not survive a round trip through `set_data!` (IS4), so trefoil is excluded
from part B.

### Part B — agreement holds only for self-loop-free diagrams

With the bending numbers pinned to the same vector, the corner sequences agree
exactly, after normalising away the representation difference (the fork emits
one control entry per bend with repeated endpoints; `_route_edges` emits
corners only).

| Case | Self-loops | Crossings | Polylines |
|---|---|---|---|
| hopf | no | identical | identical, 4/4 edges |
| borromean | no | identical | identical, 12/12 edges |
| whitehead | no | identical | identical, 10/10 edges |
| `s2xs1_calculation` | **yes** | identical | **9/12 identical; edges 1, 7, 8 differ** |

On `s2xs1_calculation` with the pinned vector:

- **edges 1 and 7** (the self-loops) differ by representation only. The fork
  returns a 4-point cubic control set, which read as a polyline is a triangle;
  Julia emits a 4-corner rectilinear loop. Julia's form is what this project
  wants, and the fork's is not evidence of a Sage defect (R2).
- **edge 8 differs materially, and Julia is wrong.** The fork routes it to
  `(-5,-1),(-5,0)`, ending on crossing 6. Julia routes it to `(-3,-1),(-3,0)`,
  ending on crossing 3's position. This is IS7, and for this case it is a
  Julia-specific defect: the fork gets it right.

---

## 5. Divergences

**IS1 — MILP 1 has no unique optimum.** `minimal_bending_numbers` and the first
`MixedIntegerLinearProgram` in `allocate_pos.py` minimise a sum with no
secondary criterion; distinct optima are drawn as visibly different diagrams.
Evidence: section 4 part A. **Accepted, not scheduled for repair.** Tests must
therefore compare the objective value, never the solution vector.

**IS2 — the fork is bound to an integer-returning backend.** It dropped
upstream's `convert=ZZ, tolerance=1e-3`. Under GLPK the raw `get_values`
returns floats and the following `range(abs(s[...]) + 1)` raises `TypeError`
(reproduced). Fork-introduced; upstream is fine.

**IS3 — the fork's geometry is not machine-independent.** It also dropped
upstream's `solver=` parameter, so the Sage default decides. Fork-introduced;
upstream exposes the choice. Largely moot now that the MATLAB path is not the
oracle and cannot run on Windows.

**IS4 — `set_data!(v; pd=...)` does not preserve PD labelling.**
`[1 2 4 3; 3 4 6 5; 5 6 2 1]` round-trips to `[5 2 6 3; 3 6 4 1; 1 4 2 5]`
because the code is rebuilt through `_pd_to_gauss` and re-derived by `pd_code`.
Per-edge data carried alongside a PD matrix — bending numbers, replay
positions — must be indexed to what `pd_code` returns.

**IS5 — the Julia port normalises a Sage bucketing quirk.** In the MILP 2
region walk, upstream applies the modulo to the first branch only
(`link.py:3767-3773`, copied verbatim into the fork), so a segment reached with
`direction` in `{5,6,7}` is dropped from all four constraint sums.
`_piece_lengths` applies `% 4` uniformly. Instrumenting the walk shows
`direction` peaks at 4 on hopf, trefoil, borromean and whitehead, and 4 is
handled identically by both. There is also an argument that the bad-region
elimination leaves exactly four left turns per region, which would make 5..7
unreachable; that has not been proved. Latent, not currently firing.

**IS6 — `planarize_virtual_gauss` runs twice per plot.** `orthogonal_layout`
calls it for `real_only` links and `_real_arc_layouts` calls it again. Results
agree; the work is duplicated and the call sites can drift.

**IS7 — `_route_edges` misroutes the edge that closes a cycle.** The routing
walk assigns a crossing position on first discovery:

```julia
elseif target in available_crossings
    crossing_positions[target] = last(points)
```

For an edge returning to an already-positioned crossing, nothing verifies that
the polyline's last point equals that position. It is assumed to close, and
when it does not, the failure is silent. Measured over the 19 `VLExample` plot
calls, exactly one edge per affected diagram fails, always the head:

```
koda_virtual        edge 14 head: got (-1,-1)  want (2,2)
s2xs1_calculation   edge  8 head: got (-3,0)   want (-5,0)
koda_fig13          edge  5 head: got (1,0)    want (-1,0)
```

Verified by hand on `koda_fig13` (4 real crossings, 6 after planarisation):
edge 5 runs from crossing 5 at `(0,1)` and should reach crossing 2 at `(-1,0)`,
but ends at `(1,0)`; the other 11 edges are correct.

The failure is bending-vector dependent. On `s2xs1_calculation` the pinned
vector triggers it and the LP-solved vector does not, so **deleting the
override in `examples/VLExample.jl` removes one of the three**. Whether the same
defect exists in the fork has not been tested; part B shows the fork routing
edge 8 correctly on the one case compared, so IS7 is at least partly
Julia-specific.

**IS8 — the MATLAB plot path cannot run on Windows.** Section 7.

---

## 6. VLExample status

`examples/VLExample.jl`, the equivalent of `topology/Manifold/VLExample.mlx`.

| Check | Result |
|---|---|
| Runs without exception | 19 / 19 |
| SVG written | 19 / 19 (315–6129 bytes, none empty) |
| P1 — every segment axis-parallel | 19 / 19 |
| P5 — crossing positions distinct | 19 / 19 |
| P3 — each edge ends on its own crossing | **16 / 19** |
| P6 — no meeting except at a crossing | 16 / 19 (the same three) |

Failing: `koda_virtual`, `s2xs1_calculation`, `koda_fig13`. All three are IS7;
P6 failures are a consequence of the misrouted endpoint, not an independent
defect. Deleting the `s2xs1_calculation` override takes this to 17 / 19.

`VISUAL_VERIFICATION.md` records all 19 as manually inspected and passing.
The manual pass missed all three. A five-line endpoint assertion did not.

---

## 7. The Windows constraint

Source: Windows-side logs supplied 2026-08-20
(`sage-windows-matlab-logs-20260820.zip`) — a ChatGPT session from
2026-02-21/22 and 2026-03-19, and a Codex session from 2026-03-31.

**`pyenv` cannot reach a WSL Sage.** MATLAB's `pyenv` embeds a *Windows native*
CPython in the MATLAB process. WSL's Sage is a Linux ELF CPython. The logs
state the combination is "structurally incompatible". `SageWrapper.m` depends
on exactly this embedding, plus a persistent `sage.all.__dict__` namespace
holding diagram state under `sageName`, so the bridge has no Windows form.

**The external-process fallback was attempted and did not land.** Calling WSL
Sage through `wsl.exe` hit `PermissionError: [Errno 13]` on
`/home/akaghef/dev/sageio`, `bash: syntax error near unexpected token '('` from
nested quoting in `sage -c`, a `NameError` from broken filename quoting, and a
Windows/WSL path mismatch for `job.sage`. These are individually solvable, but
the result is a stateless RPC, which means rewriting `SageWrapper`.

**Consequences.**

1. The canonical Windows `REdgeTable.Position` fixture **cannot be produced**
   under the current design. Every gate that waits on it is void, not pending.
2. `CANONICAL_BOUNDARY.md` requires the Windows Harness to rerun and compare
   structural, numerical **and visual** evidence. The visual clause is
   unexecutable on Windows for the MATLAB path.
3. Julia is therefore the only path that can produce plots on the canonical
   platform. The dependency set is free of Python, Sage and system libraries,
   so this is expected to work — but **it has not been verified on Windows**.
   The one thing to check is whether
   `julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'` succeeds
   there. If `GLPK_jll` fails to resolve, `HiGHS` or `Cbc` are drop-in
   alternatives.

Two escape routes exist and are not recommended: a stateless `wsl.exe` RPC
(rewrites the bridge, invests in a path being retired), or running MATLAB
itself inside WSL (the logs call this architecturally cleanest, but note
licence-activation friction, and it moves the Harness off the Windows-native
location the ADR names).

---

## 8. Not verified

- **P3 against the fork's output.** Valid to run, because control-polygon
  endpoints lie on the curve. Would establish whether IS7 has a counterpart on
  the MATLAB side. Cheapest remaining comparison.
- **`_real_arc_layouts` and the SVG conventions.** Never compared with
  `VirtualLink.calcPos`. The under-strand gap of `0.20` in pre-projection
  units, the `0.3 * scale` corner radius, and the arrow placed at the midpoint
  *index* rather than the arc-length midpoint are all unexamined.
- **`planarize_virtual_gauss`.** All part-A/B cases were supplied as ordinary
  Gauss codes, so the planariser was a pass-through. 22 of the 45 plot
  instances in the acceptance corpus take that path.
- **The Whitehead visual rejection.** Its bending vector equals the one the
  fork's default backend returns, and it passes P1/P3/P5/P6, so neither IS1 nor
  IS7 explains the 2026-08-14 rejection.
- **Julia on Windows.** Section 7.

---

## 9. To-be

Sage leaves the plot path; Julia owns every stage.

| Layer | as-is | to-be |
|---|---|---|
| `pd_code`, `regions` | Sage | Julia (verified) |
| MILP 1, MILP 2 | Sage MILP, fork default backend | Julia JuMP, degeneracy accepted |
| Routing | `allocate_pos.py` | Julia, after IS7 is fixed |
| Isolated components, unknot circles, abalone reversal | `allocate_pos.py` | Julia (**unported**) |
| Drawing | MATLAB | Julia SVG |
| Sage | required, Mac-only | not required |

**PR1 — assert P3 in `orthogonal_layout`.** Five lines: each edge's polyline
must start at its tail crossing and end at its head crossing. Failing loudly
beats returning silently broken coordinates. This is the check that found IS7.

**PR2 — delete the `bending_numbers` overrides in `examples/`.** They were
hand-found to dodge the fork's default solution and are harmful here: on
`s2xs1_calculation` the override causes IS7 and the LP-solved vector is clean.

**PR3 — fix IS7.** `koda_fig13` is the smallest reproducer at 4 real crossings.
Determine whether `_piece_lengths` gives the closing edge the wrong length or
`_route_edges` gives it the wrong turn/direction.

**PR4 — verify Julia on Windows.** Section 7.

**PR5 — port the missing layout features.** `_isolated_components`,
`component_gap`, unknot circles, abalone reversal.

**PR6 — decide IS5 rather than inherit it.** Either match the fork's bucketing
or keep the normalisation, and record which, with the reason, at the call site.
Settling the reachability argument in IS5 may close this for free.

Suggested order: PR1, PR2, PR3, PR4 in parallel with them, then PR5, PR6.

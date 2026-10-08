# Session notes: half-section (0.5 x 0.5) debugging, quad mesh, matched loads

Companion to `README.md` (not a replacement). Records what was found, what was
changed, and what is still open, from the debugging session on the Abaqus 1D
(B31) vs Abaqus 3D (C3D8R) vs FEniCS warping-solver comparison.

## 1. Code changes made

| File | Change |
|---|---|
| `loads_io.py` (new, in `batch/`) | `export_accepted_loads(summary_csv, out_csv)` writes accepted `case,Vx,Vy,Mx`; `load_fixed_loads(csv)` reads them back as the `draws` DataFrame. |
| `batch_driver.py` | Imports `loads_io`; new CLI flag `--loads-csv PATH`; the `draws` block uses `load_fixed_loads()` when the flag is given, otherwise the original `rng.uniform(...)` draw; after `summary.csv` is saved, `accepted_loads.csv` is written next to it. |
| `warping_core.py` (`FEniCSfiles/`) | (a) node displacement evaluation uses `x_sol.compute_vertex_values(mesh)` instead of `x_sol(Point(...))` per node; (b) the `x1=0` line query takes mesh vertices on that line (falls back to `Point()` sampling if none); (c) VTU export picks the meshio cell type (`triangle`/`quad`) from the mesh instead of hardcoding `"triangle"`. |

Added then fully reverted (per request): Abaqus licence-retry logic in
`abaqus_io.py` / `batch_driver.py`; the `MeshEditor(order=False)` change to
`load_mesh()` in `warping_core.py`.

Edits made by the user (not by Claude): `*Transverse Shear Stiffness` in
`1d_fx_fy_mx_half_ar.inp` changed `0.4167 -> 0.1042`; duplicate `*Static` card
removed from `3d_fx_fy_mx_half_ar.inp`; `config.py` repointed to the
`_half_ar` / `_finer` templates and `slice_ar_50_quad.xdmf`.

**Status:** the three `warping_core.py` edits are syntax-checked but have not
yet been confirmed by a full end-to-end case run. Re-test with
`--cases 2` before launching a full batch.

## 2. Bugs and findings

- **Stale shear stiffness (real bug, fixed):** the 0.5 section kept the 1x1
  value 0.4167; correct rescale is /4 = 0.1042. It was the only
  section-halving constant found wrong in a full audit.
- **Duplicate `*Static` keyword** in the 3D template made Abaqus reject the
  input (fixed).
- **False alarms, retracted:** apparent 3D/FEniCS coordinate mismatch (the
  `*Instance` translation `0.25, -0.25, 0` had not been applied) and apparent
  off-centroid coupling point (`m_Set-1` is a separate assembly node at
  `(0,0,100)`).
- **DOLFIN 2019.1.0 and quad meshes:** `Function.__call__`/`Point()` lookup is
  simplex-only ("Intersection is only implemented for simplex meshes"), and
  `mesh.order()` can fail with "Cell is not orderable" for externally generated
  quad meshes. A dolfin-native `UnitSquareMesh.create(n, n,
  CellType.Type.quadrilateral)` (then `scale` + `translate` to centre on 0)
  avoids the ordering problem; the post-processing edits above remove the
  point-lookup dependency. Arbitrary-point queries on quads remain unsupported.
- **Slenderness compensation:** with plain `K11,K22` (no `SCF` label) Abaqus
  still applies an internal factor, so `SF/SE = f_p * K`, observed `= K/2` for
  direction 3. With `w = 12 E I / (K L^2)`: old (1x1, L=2) `w = 0.600`, new
  (0.5, L=1) `w = 0.600`, so the factor is unchanged by halving the section
  and the element length together. The exact constant mapping `w -> f_p` was
  not confirmed from documentation.
- **Open anomaly:** `SF2/SE2 = 0.25000` in every case in both
  `cases_half_ar_seed1` (K=0.1042) and `cases_half_ar_seed12_I_fine`
  (K=0.4167), independent of K, while `SF3/SE3 = K/2`. Not explained.
- **Recommended K for 0.5 x 0.5:** `K = (5/6) G A = 0.0801` (not applied). The
  original 0.4167 looks like `(5/12) E A` (E used where G belongs), implying
  kappa = 1.083 > 1.
- **3D runtime:** the `_finer` 3D template halves the length-wise spacing
  (dz 1.0 -> 0.5), doubling elements (10,000 -> 20,000) and improving the
  aspect ratio from 20:1 to 10:1. Convergence is clean (about 18 increments);
  slowness is problem size times case count at `MAX_PARALLEL=3`.

## 3. How to reuse loads between runs

```bash
python3 batch_driver.py --out-dir cases_new \
    --loads-csv /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_half_ar_123/accepted_loads.csv \
    --cases 1000
```

- Use an absolute path for `--loads-csv` (relative paths resolve from your
  shell's current directory, not `config.PROJECT`).
- Replay small section -> large section: a load accepted on the smaller
  section should also pass the larger section's strain gates.
- Same seed reproduces the same candidate loads only if `config.py` scalars
  and `--cases` match. The pre-run gate is closed-form; the post-run gates use
  solved Abaqus output and can differ if the `.inp`/mesh differ.
- "Finer 1D" run needs `INP_1D = ..._finer.inp` and `INP_3D =
  3d_fx_fy_mx_half_ar.inp` in `config.py`; "finer 3D" is the reverse.

## 4. Folder state noted

- `cases_half_ar_seed1`: 46 case dirs, all 1D/3D jobs complete, but
  `summary.csv` has 6 rows (written once at the end of `main()`; a later run
  was interrupted).
- `cases_half_ar_seed12_I_fine`: still uses the old `K = 0.4167`.
- `cases_half_ar_123`: 1000 accepted loads in `accepted_loads.csv`.

## 5. Open items

1. Re-test one case end to end after the `warping_core.py` edits.
2. Decide whether to apply `K = 0.0801` to the 0.5 templates.
3. Resolve the direction-2 (`SF2/SE2 = 0.25`) anomaly, e.g. a test case with
   `*Preprint, model=YES` to print Abaqus's computed section stiffness.
4. Regenerate `summary.csv` for `cases_half_ar_seed1` from its 46 completed
   cases, or rerun.
5. `compare.py` reads the exported `.vtu` files back; confirm it handles quad
   cells.

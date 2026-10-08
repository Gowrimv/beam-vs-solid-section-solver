# Beam vs solid vs warping batch pipeline

Compares three models of the same clamped cantilever, at a mid-span section (Z_MID = 50):

| Model | Tool | Provides |
|---|---|---|
| **1D beam** | Abaqus B31 | Section strains SE/SK, forces SF/SM and rotations UR at 11 stations |
| **3D solid** | Abaqus continuum | Stress S and strain LE at Z_MID; NFORC cut resultants at every station; ELSE energy |
| **Warping section** | FEniCS (`warping_core.py`, in the parent `FEniCSfiles` folder) | Full 3D stress/strain field on the 2D section, driven by the 1D beam's strains and forces |

Loads are random tip loads (Vx, Vy, Mx) drawn from the ranges in `config.py`. Section: square of side `config.SECTION_SIDE` (currently **0.5**), L = 100, E = 1, nu = 0.3.

---

## 1. Pipeline, per case

```
draw (Vx,Vy,Mx) ──► pre-run strain screen ──► write 1D + 3D .inp ──► Abaqus 1D + 3D (in parallel)
                                                                        │
            extract_odb.py (abaqus python) ◄────────────────────────────┘
                 │  beam_whole.csv, solid_mid.csv, solid_resultants.csv, solid_energy.csv,
                 │  beam_section_points.csv
                 ▼
        1D-vs-3D force/moment check (all stations)      ── force_moment_1d3d_error.csv
                 ▼
        ADMISSION CHECK (1D beam only) ──fail──► REJECTED_admission.json, case dropped
                 ▼ pass
        FEniCS warping solve at Z_MID                   ── *_ref.vtu, *_resultants.csv, *_section_X0.csv
                 ▼
        FEniCS vs Abaqus-3D comparison                  ── errors.csv, force_error.csv, *.png
                 ▼
        post-run gates in main() ──► summary.csv row (or dropped, and a replacement case is drawn)
```

The batch keeps drawing new cases until `--cases` of them have been accepted.

---

## 2. Files

### Driver and modules

| File | Role |
|---|---|
| `batch_driver.py` | Main driver. Handles resume, the per-case pipeline (`run_case`), draw screening and top-up (`main`) |
| `config.py` | Paths, section/material constants, load ranges, thresholds, non-dimensional scales. Stdlib only |
| `strain_checks.py` | All checks and rotation helpers: `strain_screen`, `bending_dominance`, `strain_admit`, `rotated_beam_forces`, 1D-vs-3D resultant check, stress invariants |
| `abaqus_io.py` | Writes the .inp files, submits jobs, WSL/Windows paths |
| `extract_odb.py` | ODB extraction (run with `abaqus python`) |
| `fenics_solve.py` | Builds the FEniCS problem once and solves it per case |
| `compare.py` | FEniCS vs Abaqus-3D errors and plots |
| `plotting.py` | End-of-batch plots |
| `reprocess_checks.py` | Old import name for `strain_checks`; kept so old imports still work |

### Standalone scripts

| Script | Needs | Use |
|---|---|---|
| `batch_driver.py` | Abaqus + FEniCS + pyvista | `python3 batch_driver.py --cases 1000 --seed 42 --out-dir DIR [--workers 3] [--max-post-strain 5e-3] [--min-bending-dominance-post 1] [--min-bending-dominance X] [--redo-fenics] [--redo-compare] [--dry-run]` |
| `reextract.py` | Abaqus python (reads ODBs only, no licence) | Re-run extraction on existing cases: `--cases-dir --limit --workers --dry-run --force` |
| `repost_process.py` | numpy/pandas | Rewrite `force_moment_1d3d_error.csv`, `beam_1d_mid.csv` |
| `recompare.py` | pyvista, matplotlib | Rewrite `errors.csv`, `force_error.csv` and plots, without re-solving |
| `rebuild_summary.py` | numpy/pandas | Rebuild `summary.csv` from case folders; loads come from each case's own `.inp`. `--admitted-only` re-applies the checks |
| `plot_results.py` | matplotlib | Redraw the end-of-batch plots from `summary.csv` |
| `freebody_cut.py` | Abaqus python | Checks one case's 3D section resultants two ways (NFORC free-body vs equilibrium transport). `--inp --area --stations --csv` |
| `run_freebody.sh` | WSL | Wrapper: `bash run_freebody.sh <case folder or id> [stations]` |
| `sweep_freebody.py` | WSL + Abaqus | `freebody_cut.py` over a sample: `python3 sweep_freebody.py CASES_DIR --sample 40` (folder required) |
| `recover_section.py` | FEniCS | Single-case FEniCS solve from `beam_mid.csv` |
| `monte_Carlo.py` | numpy | Screen pass rate and bending dominance for the current `config.py` section and ranges, plus the exact range safety check: `python3 monte_Carlo.py [N]` |

**Out of date: don't rely on these as they stand**

- **`recover_section.py`:** applies the `R.T @` rotation to the strains (see section 6).
- **`extract_freeBody.py`:** its Abaqus calls don't work; `freebody_cut.py` replaces it.

**Re-running FEniCS on existing folders:** `batch_driver.py` draws the loads again from `--seed`. Since the `strain_screen` sign fix, those draws no longer line up with existing `case_XXXX` folders. So use `recompare.py`, `rebuild_summary.py` and `repost_process.py` on old folders, not `--redo-fenics` or `--redo-compare`.

---

## 3. Checks and admission

| # | When | Check | Criterion | Rejects? |
|---|---|---|---|---|
| 0 | Batch start | `range_safety_check` | Exact worst case over the whole load box, at the root and the tip, is below 5e-3 | Warning only |
| 1 | Pre-run | `strain_screen` | Analytical root principal strain, using the bound \|Mx\| + \|Vy\|·L, <= `STRAIN_THRESHOLD` (5e-3) | Yes |
| 2 | Pre-run | `bending_dominance` (theory, Z_MID) | >= `--min-bending-dominance` | Only if that flag is given |
| 3 | Run | Abaqus / extraction / FEniCS failure | — | Yes, each with its own reason (below) |
| 4 | After Abaqus, before FEniCS | `strain_admit`: `bending_ok` | σ_max/τ_max at Z_MID >= 1, from 1D SF/SM: σ_max = \|SF1\|/A + (\|SM1\|+\|SM2\|)·c/I, τ_max = 1.5·\|V\|/A + \|SM3\|/(0.208·a³). The old strain ratio is logged as `bending_dominance_strain_mid_abaqus` | Yes |
| 4 | After Abaqus, before FEniCS | `strain_admit`: `small_ok` | Max **1D** principal strain over all 11 stations < `--max-post-strain` (5e-3) | Yes |
| 5 | In `main()` | `root_strain_under_threshold` | 1D principal strain at Z=0 <= 5e-3. Re-checks 4 for cases reloaded from disk | Yes |

**Principal strain** is ε/2 + √((ε/2)² + (γ/2)²).
- **1D:** ε = \|SE1\| + (\|SK1\|+\|SK2\|)·c, γ = \|SE2\|+\|SE3\|. Twist SK3 is not included.
- **3D:** the corner and centreline strains (`eps_3d_corner_post`, `eps_3d_centreline_post`) are **reference only**, not gated.

**Discard reasons** returned by `run_case()`, counted separately and listed at the end of the batch:

- `abaqus_failed`
- `extract_failed`
- `fenics_not_built`
- `fenics_diverged`
- `fenics_error`
- `rejected_admission`
- `rejected_previously`

**`REJECTED_admission.json`** stores the thresholds and `ADMISSION_RULE_VERSION` (currently 2). A rerun skips a rejected case only if both match. If you change the admission logic, increase the version number so old rejections are checked again.

---

## 4. What is compared

**FEniCS vs Abaqus 3D, on the Z_MID section**, nodes paired directly (the meshes coincide):
- **`errors.csv`:**
  - S11..S23, Mises, I1 and J3, each with L2, Linf, relL2, relL2_resultant and L2_over_scale.
  - Energy: FEniCS Φ against the 3D ELSE per unit length.
- **Plots:** `contours.png`, `centerline_overlay.png` (stress along X=0), `strain_overlay.png` (LE along X=0).

**FEniCS vs the 1D beam:** `N_resultant`/`M_resultant` rows of `errors.csv`, and `force_error.csv`.
- These agree by construction: `warping_core` forces FEniCS n0/m0 to equal the 1D forces through Lagrange multipliers (α, β). So they test the solver tolerance, not the model.

**1D vs 3D:** `force_moment_1d3d_error.csv`, all stations.
- N, Vx, Vy, Mx, My, Mz: the beam's SF/SM rotated to global, against the 3D NFORC cut resultants.
- Error measures: `pct_diff`, `pct_of_resultant` (use this one), `*_norm`.
- For information only; never rejects a case.

**Analytical vs 1D/3D:** printed at the end of the batch.
- `eps_mid_analytical` against `eps_1d_mid_post` and `eps_3d_corner_post`.
- Bending dominance: formula vs Abaqus.

### Error measures

| Name | Definition |
|---|---|
| L2 | RMS over section nodes (not area-weighted) |
| Linf | Max absolute nodal difference |
| relL2 | L2 / RMS of that component in 3D. Unreliable for near-zero components (S13, S23) |
| relL2_resultant | L2 / RMS of 3D Mises (Mises³ for J3). The same scale for every component |
| L2_over_scale | L2 / σ₀ (σ₀³ for J3) |

### Normalisation scales (`config.py`)

| Scale | Formula | 1x1 | 0.5x0.5 (current) | Used for |
|---|---|---|---|---|
| F₀ | EI/L² | 8.33e-6 | 5.21e-7 | Forces |
| M₀ | EI/L | 8.33e-4 | 5.21e-5 | Moments |
| σ₀ | E·ε_char (5e-3) | 5e-3 | 5e-3 | Stresses, Mises, I1 |
| σ₀³ | — | 1.25e-7 | 1.25e-7 | J3 |
| U₀ | F₀ | 8.33e-6 | 5.21e-7 | Energy per unit length |

With the 1/8 load scaling, loads span about ±1.2 F₀ and ±1.2 M₀ for the 0.5 section (±0.6 for 1x1).

Shear judgement: shear is much smaller than axial stress, so its relative errors are naturally larger. Judge it by `relL2_resultant` or `L2_over_scale`.

---

## 5. Frames and rotations

`R_local = rotation_matrix_from_ur(UR1, UR2, UR3)` is the Rodrigues matrix of the beam node's rotation vector, taken at that station.

| Quantity | Frame / transform |
|---|---|
| Abaqus beam SE, SK, SF, SM | Beam-local (rotating) axes |
| Beam SF/SM to global | `R @ [SF3,SF2,SF1]`, `R @ [SM1,SM2,SM3]`. Verified: `beam_mid.csv`, 21.3° rotation, axial 2.5e-8 vs 2.8e-6 with the transpose |
| 3D solid S, LE, NFORC | Global axes; node positions undeformed (used only for pairing nodes) |
| 3D NFORC moments | Deformed lever arms `(x+u) − c_def`, full r × F, about the deformed cut centroid (`moment_config = deformed`) |
| FEniCS input | Strains `R.T @ [SE3,SE2,SE1]`, `R.T @ [SK1,SK2,SK3]`; force/moment targets raw SF/SM |
| FEniCS solve | Local material frame. `postprocess` rotates the output to global: stresses and strains as `R σ Rᵀ`, forces as `R n` |
| FEniCS vs Abaqus 3D | Both in global axes, so they are consistent. Global S13/S23/S22 contain rotated bending stress; Mises, I1 and J3 don't depend on the axes |

---

## 6. Known issues / open items

1. **Strain frame (by design).** FEniCS receives the strains as `sv_mat = R_local.T @ sv_sp` and `sk_mat = R_local.T @ sk_sp` (the pull-back in `warping_core`'s FRAMES note). The force/moment targets are passed unrotated. `rotated_beam_components()` uses the same `R_local.T`.
2. **Both strains and forces are enforced.** `warping_core` enforces n0 = n0_t and m0 = m0_t through the Lagrange multipliers α and β, on top of the prescribed strains.
   - **α = β = 0:** the prescribed strains and forces agree with the FEniCS section stiffness.
   - **Otherwise:** the multipliers adjust the field so both hold.
   - α and β are not currently logged.
   - Because the forces are enforced, the FEniCS-vs-1D force comparison is ~0 by construction.
3. **Bending-dominance check:** it is a strain ratio (the stress ratio is 2.6x larger), uses average shear, and ignores torsion (SM3/SK3). In `beam_mid.csv`, torsional shear is about 7x the transverse shear. A stress-based version including torsion has been proposed.
4. **Older case folders:** they keep old outputs until re-processed. The section-point shear used to be doubled in `beam_section_points.csv`, and old `solid_resultants.csv` files lack `moment_config`. Run `reextract.py`, then `repost_process.py`, `recompare.py` and `rebuild_summary.py`.

Backups of pre-edit files: `_backup_20260925/` (latest) and `_pre_fix_backup/`.

---

## 7. Changing the section size

`config.py` is set up so the section is **one number**: `SECTION_SIDE` (currently 0.5). `A`, `I`, `c`, `G`, the load ranges and all the scales are derived from it.

### What changes physically (1 x 1 to 0.5 x 0.5)

| Quantity | 1 x 1 | 0.5 x 0.5 | Factor |
|---|---|---|---|
| A | 1 | 0.25 | 1/4 |
| I | 0.0833 | 0.00521 | 1/16 |
| c | 0.5 | 0.25 | 1/2 |
| J (0.1406·a⁴) | 0.1406 | 0.00879 | 1/16 |
| c/I (bending strain per unit moment) | 6 | 48 | **8x** |
| 1/(GA) (shear strain per unit shear) | 2.6 | 10.4 | 4x |
| F₀, M₀ | 8.33e-6, 8.33e-4 | 5.21e-7, 5.21e-5 | 1/16 |

**Load ranges.** They are defined for 1 x 1 (V ±5e-6, Mx ±5e-4) and scaled by `SECTION_SIDE³`, which keeps the bending strain distribution identical. Checked on 400,000 draws:

| Section, ranges | Pass pre-screen | Median strain (passed) | Median bending dominance |
|---|---|---|---|
| 1x1, V ±5e-6, Mx ±5e-4 | 62.4% | 3.73e-3 | 179 |
| 0.5, same ranges (not scaled) | 0.14% | — | — |
| **0.5, V ±6.25e-7, Mx ±6.25e-5 (used)** | **62.4%** | **3.73e-3** | **358** |
| 0.5, scaled 1/16 instead | 100% | 2.25e-3 | 362 |

- **Same strains, different behaviour:** at the same strain the smaller section has half the shear strain, so bending dominance doubles.
- **Rotations double:** curvature κ = ε/c doubles, so section rotations roughly double (about 21° to about 40° at mid-span for the `beam_mid.csv` case). Abaqus nlgeom may need more increments.
- **To keep rotations instead of strains:** scale the loads by `SECTION_SIDE⁴` (1/16). Strains then halve and every draw passes.

Run `python3 monte_Carlo.py` after any change to see the pass rate and the exact range safety check.

### Already done in the code

- `config.py`: `SECTION_SIDE = 0.5`; A, I, c and the ranges are derived from it.
- `abaqus_io.write_inp`: divides the 3D TRVEC value by `config.A`, because TRVEC is a traction (force per area). *This was only correct before because A = 1.*
- `freebody_cut.py`: `--area` defaults to `config.A` (read from `config.py` in the same folder), so `run_freebody.sh` and `sweep_freebody.py` need no extra argument.
- `monte_Carlo.py`: uses `config` and `strain_checks` directly; no private constants.
- `batch_driver.py`: at start-up it reads the `*Beam Section, SECTION=RECT` dimensions from `INP_1D` and **stops** if they don't equal `SECTION_SIDE`. It skips the check with a message if the file can't be read or parsed, and never stops a `--dry-run`.
- Automatic from `config` or the mesh: `strain_screen`, `bending_dominance`, `range_safety_check`, the post-run strain checks, `extract_odb.py`, `compare.py`, `plotting.py`, `rebuild_summary.py`, `recompare.py`, `repost_process.py`, and `warping_core.py` (area, centroid and I come from the mesh).

### Still to do by hand (outside this folder, in `PROJECT`)

- [ ] `1d_fx_fy_mx_half_ar.inp`: beam section dimensions, `1., 1.` to `0.5, 0.5` on the line after `*Beam Section, ..., section=RECT`. The start-up check enforces this.
- [ ] `3d_fx_fy_mx_half_ar.inp`: new part geometry and mesh, still centred on X=0, Y=0.
  - Keep an **odd** number of nodes across, so there is a node at the centre.
  - Keep a sensible element aspect ratio: if the in-plane element size halves, consider halving the axial size too.
  - Not checked automatically.
- [ ] Keep the set and load line formats the regexes in `abaqus_io.py` depend on:
  - `free, 1/2/4, value` (1D `*Cload`)
  - `free, TRVEC, value, 1., 0., 0.` and `..., 0., 1., 0.` (3D `*Dsload`)
  - `m_Set-1, 4, value` (3D moment on the coupling node)
- [ ] Re-slice the new 3D solid at Z=50 into `slice_ar_50.xmf` (`config.MESH_XMF`). Its nodes must coincide with the Abaqus nodes: `compare.py` warns if any pair is more than 1e-4 apart. Not checked at start-up.
- [ ] Pass a new `--out-dir` (now required): 1x1 and 0.5x0.5 results must not be mixed.

**Going back to 1 x 1:** set `SECTION_SIDE = 1.0` and point `INP_1D`, `INP_3D` and `MESH_XMF` back at the 1 x 1 files.

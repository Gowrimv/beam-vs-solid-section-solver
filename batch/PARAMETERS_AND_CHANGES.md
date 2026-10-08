# Beam 1D / solid 3D / FEniCS section study: parameters and change log

Last updated **8 Oct 2026**. This is a running record of what each run used, what changed between runs, and why.

Sources:
- this chat's analysis of the case folders;
- `batch/README.md` (pipeline, checks, frames);
- `batch/Claude outputs/beam_comparison_summary.tex` (an earlier chat on the shear-stiffness and section-halving bugs).

For the pipeline itself (checks, frames, error measures), see `batch/README.md`. This file covers only parameters and the reasons for changing them.

---

## 1. Current set-up (for the paired runs `cases_EI_1x1` / `cases_EI_half`, started 7–8 Oct)

The earlier `E = 1e6` set-up is obsolete (§2.7, §2.13, §2.17). Current values:

| Parameter | Value | Where | Why / see |
|---|---|---|---|
| Length L, mid-span Z_MID | 100, 50 | `config.py` | Unchanged; all errors are evaluated at Z_MID |
| Section | square, side h = 1 or 0.5, chosen with the `BATCH_SECTION_SIDE` environment variable | `config.py` | One config serves both sizes |
| E, ν | **E = 1.2e5 / h⁴** (1.2e5 at h = 1, 1.92e6 at h = 0.5), ν = 0.3 | `config.py`, templates | E·I = 1e4 for every size, so F0 = EI/L² = 1 and M0 = EI/L = 100. §2.17 |
| Loads | v_x, v_y ∈ ±0.6·F0; m_tip ∈ ±0.7·M0, drawn per the gate table (seed 123); 1000 cases + 300 spares | `loads_table.csv` (made by `make_table_loads.py`) | One file serves both sizes. §2.14 |
| 1D model | B31, 100 elements, NLGEOM; transverse shear stiffness = Abaqus default | `1d_fx_fy_mx_E1p2e5.inp` (h = 1), `1d_fx_fy_mx_half_ar_EI1e4.inp` (h = 0.5) | §2.2, §2.9 |
| 3D model | C3D8R, 10×10 across the section, 10:1 element aspect ratio, NLGEOM, tip shears as non-follower tractions, Mx through a distributing coupling | `3d_fx_fy_mx_E1p2e5.inp`, `3d_fx_fy_mx_half_ar_finer_EI1e4.inp` | §2.12 |
| FEniCS section mesh | 10×10 quads | `slice_50_quad_native.xdmf` (h = 1), `slice_ar_50_quad.xdmf` (h = 0.5) | — |
| FEniCS solve | at **E = 1**; forces divided by k = E/1 before, stresses multiplied by k after | `config.FENICS_E`, `fenics_solve.py` | The solve does not converge with a large E. §2.15 |
| **Strain input to FEniCS** | **Abaqus SE/SK as written (beam-local director basis); no Rᵀ pull-back** | `batch_driver.py` | **§2.19. Changed 8 Oct; affects every earlier FEniCS result** |
| Shear-strain correction | on: SE2 = SF2/K, SE3 = SF3/K, K = 0.8497·G·A | `CORRECT_SHEAR_STRAINS` | Harmless; left on. §2.16–2.17 |
| Force/moment constraints in the section solve | on | `warping_core.py` | With the fix their multipliers are ~1e-6, so nearly inactive. §2.19 |
| Pre-run gates (loads only) | table strain gate (c/L)·max(\|m−v_y\|+\|v_x\|, \|m\|) ≤ 0.005; shear/bending ratio ≤ 0.1; mid-span moment ≥ 0.1·M0 | `strain_checks.py`, `config.py` | §2.13–2.14 |
| Post-Abaqus gates | 1D residual ≤ 1% of applied load (`--max-1d-residual`); hourglass ALLAE/ALLIE ≤ 1% (`--max-hourglass`) | `batch_driver.py` | §2.11, §2.12 |
| Admission (after Abaqus, before FEniCS) | bending dominance and max 1D strain < 5e-3 (`strain_admit`) | `strain_checks.py` | Unchanged; does not use the rotated strains |
| Abaqus wait limits | 1D 300 s; 3D 2400 s | `abaqus_io.py` | §2.4 |
| Parallel cases | 3 | `MAX_PARALLEL` | 30 of 38 licence tokens |
| Error measures | `relL2_resultant` (÷‖Mises₃D‖), `relL2_Szz` (÷‖S33₃D‖), `L2_over_scale` (÷σ0 = E·c/L) | `compare.py` | `relL2` is diagnostic only. §2.10 |

Scales: F0 = 1 and M0 = 100 at both sizes; σ0 = 600 (h = 1) and 4800 (h = 0.5).

---

## 2. Change log: what changed and why

### 2.1 Section halving and the h³ load scaling (late Sep)
- The 0.5×0.5 section (`half_ar`) was introduced, and loads were scaled by side³ so the **bending strain** distribution matched 1×1. Pass rate was 62.4% at both sizes.
- **Side effect, noted at the time:** at equal strain, curvature and rotations scale as 1/side, so the 0.5 section rotates about twice as much. This turned out to drive the size effect (§2.5).
- `abaqus_io.write_inp` divides the 3D traction by A, because TRVEC is force per area. It was only correct before because A = 1.

### 2.2 Transverse shear stiffness (28–30 Sep, earlier chat)
- **Stale 0.5×0.5 value:** the 0.5 template kept the 1×1 value K = 0.4167. It was corrected to 0.4167/4 = 0.1042.
- **The 0.4167 itself is wrong in kind:** it implies k = K/(GA) = 1.083 > 1, so it was most likely 5/12·E·A instead of k·G·A. The correct 5/6·G·A is 0.3205 (1×1) and 0.0801 (0.5×0.5).
- **Abaqus softens explicit values:** for user-given K it applies a slenderness compensation factor f_p = 1/(1+w). This was confirmed from output: SF3/SE3 = K/2 for the 0.5 runs.
- **Unresolved anomaly:** SF2/SE2 = 0.25, a constant, regardless of K. Not explained.
- **Effect on results:** negligible here. Shear deflection over bending deflection is about 0.8·(h/L)², that is 8e-5 (1×1) and 2e-5 (0.5×0.5).

### 2.3 Matched-load replay (late Sep)
- `loads_io.py` and the `--loads-csv` option replay one run's accepted loads on another run instead of drawing new ones.
- Each batch now writes its own `accepted_loads.csv`.
- **Case numbering:** case *i* = row *i* of the loads file, counting only rows that pass the pre-run screen. Rejected rows shift the numbering.

### 2.4 3D wait limit 900 → 2400 s (6 Oct)
- **Found:** `cases_half_ar_3d_finer_123` has 758 cases instead of 1000.
  - 60 cases converged in Abaqus, but the driver stopped waiting at 900 s. Every accepted case took ≤ 850 s; every dropped one took 858–2032 s.
  - These are the hardest nonlinear cases (34 increments vs 20).
  - 176 cases were never run.
  - Cases 821–834 failed on a full disk (ODB "Extend failed").
  - Case 508 was killed with SIGTERM.
- **Changed:** `ABAQUS_TIMEOUT_3D = 2400`. Backup in `_backup_timeout_20261006/`.
- **Added `rerun_dropped.py`:** finishes cases whose Abaqus jobs completed (extraction, FEniCS and comparison only), then merges them into `summary.csv` after a timestamped backup. It reads loads from each case's own `.inp`, so it is safe whatever `config.py` currently says.

### 2.5 Load scaling h³ → h⁴ (6 Oct)
- **Found:** von Mises, I1 and S33 errors depend on the **absolute twist angle** |UR3|.
  - At equal twist, the 1×1 and 0.5×0.5 medians agree within about 0.5 percentage points.
  - At the equal size-independent twist/bending ratio, they don't.
  - The 0.5 tail is simply the cases with twist > 0.1 rad, which 1×1 never reaches.
- **Why:** strain ≈ rotation × c/L, so at equal strain the smaller section rotates twice as much.
- **About the twist:** no torque is applied. Mx is a **bending moment** about X (the 1D instance is rotated onto Z). Twist therefore comes only from large-rotation effects: the combined rotation vector (UR3 tracks UR1·UR2/2, correlation 0.8–0.9) and non-follower loads gaining a component about the deformed axis.
- **Changed:** `_LOAD_SCALE = SECTION_SIDE**4`, giving equal non-dimensional loads (±0.6 F0/M0), so rotations and twist match. Strains halve for 0.5×0.5. Backup in `_backup_h4_20261006/`.
- **Test run:** `cases_half_ar_h4_finer`, loads from `cases_full` × 0.5⁴, paired with `cases_full`.
  - **Correction:** despite the folder name, it ran on the **coarse** 3D mesh (1.24 MB inputs, 101 nodes along the length). `config.py` had been changed back before the run.
- **Result:** rotations matched case by case (ratio 1.000), and median errors matched `cases_full` (von Mises 6.74 vs 6.71, S33 6.58 vs 6.53, I1 7.52 vs 7.35). The size effect was rotation.

### 2.6 Remaining h⁴ tail and the 1D convergence problem (7 Oct)
- **Found:** the top ~1% of h⁴ cases (634, 953, 895, 346, 592, 489) were still 1.7–4.8 points worse in von Mises and up to 9 points worse in energy.
  - In these cases the 3D fields were exactly half the 1×1 values, as they should be.
  - The 1D beam, however, rotated 2–5% more.
- **Cause:** with E = 1 the 0.5×0.5 loads are about 1e-7, below Abaqus's absolute "zero force" level.
  - Abaqus then replaces the reference force with a fallback of 1e-2, so the force check always passes.
  - 1D solutions were accepted with residual forces of up to 300% of the applied shear.

  1D residual force at convergence, as a fraction of applied V:

  | Run | median | 90th pct | 99th pct | cases on the fallback |
  |---|---|---|---|---|
  | `cases_full` (1×1) | 0.02% | 0.12% | 0.5% | 0.3% |
  | `cases_half_ar_3d_finer_123` (0.5, h³) | 0.04% | 3.6% | 49% | 20% |
  | `cases_half_ar_h4_finer` (0.5, h⁴) | 0.2% | 23% | 92% | 55% |

- **Consequence:** FEniCS is driven by the 1D section forces, so any 0.5×0.5 run at E = 1 has some under-converged cases. This also inflated the earlier h³ tail.

### 2.7 E = 1 → E = 1e6 (7 Oct)
- **Why not change Abaqus's convergence controls instead:** the model is linear elastic and every output is normalised (F0/M0 or relative error), so scaling E and the loads together changes no result and lifts every force well above the zero-force level. Lowering the threshold through `*Controls, parameters=field` was rejected: it relies on obscure parameter meanings, has to go in every template, and leaves the 3D zero-force warnings.
- **Code check:** every hard-coded small threshold in the code (1e-9, 1e-20 and so on) is a "too small to divide by" guard, so larger values are safe.
- **New templates** in `batch_test/`; the originals are untouched:
  - `1d_fx_fy_mx_E1e6.inp`
  - `1d_fx_fy_mx_half_ar_E1e6.inp`
  - `3d_fx_fy_mx_E1e6.inp`
  - `3d_fx_fy_mx_half_ar_finer_E1e6.inp`
- **`config.py`** (backup in `_backup_Escale_20261007/`):
  - E = 1e6.
  - `BATCH_SECTION_SIDE` picks the section, the templates and the FEniCS mesh.
  - Load ranges are set directly in F0/M0 units.

### 2.8 Which loads to replay (7 Oct)
Three options were compared:

| Option | Cases | Verdict |
|---|---|---|
| A. h⁴ / `cases_full` loads, screened | 797 paired | 203 of the 1000 exceed the 1×1 strain screen (up to 0.0072), because `cases_full` predates the strain-screen sign fix |
| B. h⁴ loads, `STRAIN_THRESHOLD` raised to 7.5e-3 | 1000 paired | Keeps the high-twist tail but departs from the 0.5% limit |
| **C. `cases_quadfix` loads** | **1000 paired** | **Chosen.** Respects the 0.5% limit at 1×1; links to the new quad-mesh 1×1 run |

- Files: `loads_E1e6_1x1_from_quadfix.csv` and `loads_E1e6_half_from_quadfix.csv`. The non-dimensional loads are identical (max difference 1e-11).
- The option A files (`loads_E1e6_*_from_h4*.csv`) are kept for reference.
- `loads_h4_from_cases_quadfix.csv` (E = 1) is **obsolete**.
- **Low bending dominance:** in `cases_quadfix`, cases 47, 469 and 819 have strain-based bending dominance < 10 (5.2–8.0, against a median of 150). They pass the admission check (≥ 1).
  - They are zero-moment cases: the moment crosses zero near Z = 50 (mid-span moment / (V·50) = 0.05–0.08, against a median of 1.48).
  - Keep them, but report them separately: their errors are 12–16% von Mises and up to 34% energy.
  - At 0.5×0.5 their bending dominance roughly doubles (σ/τ ∝ 1/h at equal non-dimensional loads). That is expected, not an error.

### 2.9 Transverse shear stiffness → Abaqus default (7 Oct)
- **Considered:** explicit 5/6·G·A, Cowper k = 10(1+ν)/(12+11ν) = 0.850 (the most accurate for a solid square), or Abaqus's default.
- **Chosen:** Abaqus's default, as the simplest; it is the same rule at both sizes and matches the h⁴ run.
- **Effect:** about 1e-5 on rotations and section forces, invisible in the errors.
- **Caveat:** Abaqus's slenderness compensation depends on element length relative to section size. Both 1D meshes have element length 1, so the compensation differs slightly between 1×1 and 0.5×0.5. The effect is negligible.

### 2.10 Stress scale σ0 = E·c/L, and errors on a common per-case scale (7 Oct)
- **Why:** the loads are non-dimensional through F0 = EI/L² and M0 = EI/L. The matching stress scale is σ0 = M0·c/I = **E·c/L**, the bending stress produced by M = M0. Then S/σ0 = m for pure bending at any section size.
  - The old σ0 = E·ε_char (5e-3) left a factor c/L, so at equal non-dimensional loads the 0.5×0.5 stresses came out half the 1×1 ones.
  - Example, case 634: peak Szz/E = 1.99e-3 (1×1) vs 9.95e-4 (0.5×0.5); peak Szz/σ0 = 0.398 at both sizes.
- **`config.py`:** `STRESS_SCALE = MOMENT_SCALE * c / I`, giving 5000 for h = 1 and 2500 for h = 0.5 at E = 1e6. `J3_SCALE = σ0³` follows automatically.
  - For h = 1 the value happens to equal the old one, because c/L = ε_char = 5e-3. Only h = 0.5 changes.
- **`compare.py`:** new `errors.csv` column **`relL2_Szz`** = ‖ΔSij‖ / ‖S33,3D‖. It uses the same denominator for every component and invariant (‖S33‖³ for J3), so shear errors are judged against the governing bending stress, not their own near-zero norm.
  - The existing `relL2_resultant` (denominator ‖Mises_3D‖) does the same with Mises.
  - Synthetic check: a shear error worth 50% of the shear's own norm is 0.3% of ‖Szz‖.
- **`batch_driver.py` and `rebuild_summary.py`:** `summary.csv` now carries `relL2_Szz_*` and `relL2_resultant_*` for Mises, inv_1, J3, S33, S13 and S23. The reload path in `batch_driver.py` previously skipped `relL2_resultant_*`; it now copies it too.
- **Which numbers to use:**
  - **Accuracy statistics across cases and sizes:** `relL2_Szz_*` or `relL2_resultant_*`.
  - **Plots and absolute errors:** stresses divided by σ0, through `L2_over_scale_*`. This is comparable across sizes now.
  - Per-component `relL2_*`: diagnostic only.
- **Old folders:** their `errors.csv` files have no `relL2_Szz`, and their `L2_over_scale` used the old σ0. To regenerate without re-solving, run `recompare.py --force --cases-dir <dir>`, then `rebuild_summary.py`.
- **Backups:** `_backup_stressscale_20261007/` holds config.py, compare.py, batch_driver.py and rebuild_summary.py as they were before.

### 2.11 New checks: 1D convergence, mid-span moment floor, section-frame stresses, discard record (7 Oct)
- **1D convergence gate (on by default).** `abaqus_io.residual_check_1d` reads the last equilibrium iteration in `case_XXXX_1d.msg`.
  - It computes res_ratio_1d = max(residual force / F_ref, residual moment / M_ref), with F_ref = max(|V|, |Mx|/L) and M_ref = max(|Mx|, |V|·L).
  - **New flag `--max-1d-residual`** (default 0.01; negative disables). A case above the limit is discarded with the new reason `unconverged_1d`. Results go to `convergence_1d.csv` either way.
  - **Columns added to `summary.csv`:** `res_ratio_1d`, `res_force_1d`, `res_moment_1d`, `tavg_force_1d`, `tavg_fallback_1d` (True means Abaqus used its 1e-2 zero-force fallback).
  - Checked on old cases: h⁴ case 634 gives 0.58 with the fallback flagged; `cases_full` case 693 gives 6.7e-5.
  - `rerun_dropped.py` has the same flag. On E = 1 folders it will now drop unconverged cases unless you pass `--max-1d-residual -1`.
- **Mid-span moment floor (optional, off by default).** `strain_checks.mid_moment_nd` = (|Mx − Vy·(L − Z_MID)| + |Vx·(L − Z_MID)|) / M0, logged as `m_mid_nd` in every summary.
  - **New flag `--min-mid-moment-nd X`:** skips draws below X before running. It works in M0 units, so with equal non-dimensional loads the same draws are skipped at both sizes and the pairing survives.
  - In the `cases_quadfix` loads the median is 0.37; 59 cases are below 0.1 and 19 are below 0.05. Cases 47, 469 and 819 are at 0.013–0.018.
  - Left off because the zero-moment cases were to be kept and reported separately; filter on `m_mid_nd` afterwards instead.
- **Section-frame check.** `compare.py` takes the beam rotation R at Z_MID (`R=R_local`) and pulls the global stresses back as Rᵀ·σ·R.
  - It writes `section_frame_check.csv` with `transverse_ratio_3d` and `transverse_ratio_fen` = RMS(√(S11² + S22² + 2·S12²)) / RMS(S33) in the section frame. These are copied into `summary.csv`.
  - Synthetic check at 0.3 rad: the global-frame value is 0.096, the pulled-back value ≈ 0.
  - **Suggested acceptance:** < 0.05. It is reported only, not gated.
- **Discard record.** Every batch writes `discard_counts.csv`: the counts for every pre-run screen, run-time discard reason and post-run gate, plus the number accepted. This shows whether a gate is biasing the sample, for example by dropping the large-rotation cases.
- **`rebuild_summary.py`** also adds `m_mid_nd`, the section-frame ratios and `res_ratio_1d`. It re-reads the residual from `.msg` for older folders that have no `convergence_1d.csv`.
- **Not changed:**
  - **3D peak strain away from the clamp:** 3D strain is not gated (the gates use 1D strain at every station), and extraction covers only the mid-span section, so the clamp-corner peak never enters.
  - **Bending-gate redundancy:** that applied to the proposed gate table; the code's pre-run bending gate is off by default.
- **Backups:** `_backup_checks_20261007/` holds abaqus_io.py, batch_driver.py, compare.py, strain_checks.py and rebuild_summary.py as they were before.

### 2.12 Gates per the gate table; root-strain re-check removed (7 Oct)
- **Removed:** the post-run root-strain re-check in `main()`. It repeated `small_ok`, which uses the max 1D strain over all stations including Z = 0 with the same 5e-3 limit, so it never rejected anything new. `root_strain_under_threshold` is still written to `summary.csv`. Backup in `_backup_rootcheck_20261007/`.
- **Load generation** (random draws and top-ups only; replayed `--loads-csv` runs are unaffected): m_x,tip ~ U(−0.7, 0.7) (`LOAD_ND_M`), v_x, v_y ~ U(−0.6, 0.6).
- **Pre-run gates**, all written in m = M/M0 and v = V/F0 so they hold at any section size:
  - **Strain gate** (`strain_checks.table_strain_gate`): (c/L)·max(|m − v_y| + |v_x|, |m|) ≤ 5e-3. This replaces the Mohr-circle `strain_screen` as the reject rule; `strain_screen` is still logged as `eps_root`.
  - **Bending gate** (`table_bending_ratio`): √(v_x² + v_y²)·h / (4·L·(|m_x,mid| + |m_y,mid|)) ≤ 0.1, which is v/(400·m) at h = 1. Flag `--max-shear-bending` (default 0.1; negative disables).
  - **Mid-span floor:** |m − v_y/2| + |v_x/2| ≥ 0.1 in M0 units. This is exactly the table's (1/200)(…) ≥ 0.1·ε_tol at h = 1. It is kept in M0 units instead of strain so the same draws are dropped at every section size; written as strain, it would drop 59 cases at 1×1 but 174 at 0.5×0.5, breaking the pairing. **Now ON by default** (`MID_FLOOR_ND = 0.1`, flag `--min-mid-moment-nd`; negative disables). "Tune after a pilot."
  - **Effect on the replayed `cases_quadfix` loads:** strain gate 0, bending gate 0, floor 59 (the same 59 at both sizes). The zero-moment cases 47, 469 and 819 are among them.
- **Post-run hourglass check:** `extract_odb.py` now writes `solid_energy_history.csv` (ALLAE, ALLIE, ALLSE, ALLWK at the end of the step, from PRESELECT history output). `batch_driver.py` discards a case if ALLAE/ALLIE > 0.01 (`--max-hourglass`; negative disables), with the new reason `hourglass_energy`. If the file or the variables are missing, it warns and doesn't reject. `ALLAE_over_ALLIE` is copied into `summary.csv`.
  - **Not tested:** the Abaqus Python part can't be run here; check `solid_energy_history.csv` in the 5-case test.
- **Not added from the table:**
  - 3D peak LE, root reactions balance, and the transverse stress as a gate (the transverse stress is still reported).
  - The 1D post-run gates `bending_ok` and `small_ok` are kept.
- **Backups:** `_backup_tablegates_20261007/` holds extract_odb.py, config.py, strain_checks.py, batch_driver.py and rerun_dropped.py as they were before.

### 2.13 E = 1.2e5 per the gate table (7 Oct)
- **E:** 1e6 → **1.2e5**, so F0 = EI/L² = 1 and M0 = 100 for h = 1, as in the table's Case A. For h = 0.5: F0 = 0.0625, M0 = 6.25.
  - Loads are 0.04–60, far above the ~1e-7 where Abaqus's zero-force fallback appeared.
  - σ0 = E·c/L = 600 (h = 1) and 300 (h = 0.5).
  - Non-dimensional results are unchanged.
- **New templates** (copies of the `_E1e6` ones with only `*Elastic` changed): `1d_fx_fy_mx_E1p2e5.inp`, `1d_fx_fy_mx_half_ar_E1p2e5.inp`, `3d_fx_fy_mx_E1p2e5.inp`, `3d_fx_fy_mx_half_ar_finer_E1p2e5.inp`. `config.py` points to them; backup in `_backup_E1p2e5_20261007/`.
- **New loads files:** `loads_E1p2e5_1x1_from_quadfix.csv` and `loads_E1p2e5_half_from_quadfix.csv` (the `cases_quadfix` loads × 1.2e5, and × 1.2e5 × 0.5⁴). The `_E1e6` loads files and templates are obsolete.
- **Pre-run gates:** 941 of the 1000 loads pass at both sizes (the floor drops 59). Run with `--cases 941`; a larger number makes the batch top up with random, unpaired draws.

### 2.14 Loads generated exactly as the gate table says; 1000 cases with spares (7 Oct)
- **`make_table_loads.py` (new):** draws m_x,tip ~ U(−0.7, 0.7) and v_x, v_y ~ U(−0.6, 0.6) (seed 123) and applies the table's three gates for Case A (h = 1, L = 100, ε_tol = 0.005), drawing again until enough pass. The pass rate is 31.7%.
  - It writes `loads_table_nd.csv` (non-dimensional, with `load_id`), `loads_table_1x1.csv` and `loads_table_half.csv`. The two physical files are the same rows scaled with E = 1.2e5: F0 = 1, M0 = 100 for h = 1, and F0 = 0.0625, M0 = 6.25 for h = 0.5.
  - There are 1300 rows: 1000 plus 300 spares. All 1300 pass the batch's own gates at both sizes (max table strain 0.005 at h = 1, 0.0025 at h = 0.5). m_mid ranges from 0.10 to 0.84.
- **`batch_driver.py`:** with `--loads-csv`, every row is screened. Rows beyond `--cases` are kept as spares and replace cases rejected after the run, in file order, before any random top-up. Backup in `_backup_spares_20261007/`.
- **Pairing:** the case number is the same at both sizes until the two runs reject different cases. Pair afterwards by non-dimensional load (Vx/F0, Vy/F0, Mx/M0), which is unique per `load_id`.
- The `cases_quadfix` loads files are superseded for the new runs.

### 2.15 FEniCS always solved at E = 1; E applied in post-processing (7 Oct)
- **Why:** the FEniCS nonlinear section solve does not converge reliably with a large E. Abaqus still uses E = 1.2e5, for the zero-force reason in §2.6–2.7.
- **How** (`config.FENICS_E = 1.0`; `fenics_solve.py`):
  - The FEniCS problem is built with E = 1.
  - Before each solve, the 1D force and moment targets are divided by k = E/FENICS_E (1.2e5). The section strains passed to FEniCS don't depend on E, so they go in unchanged.
  - **After the solve, multiplied by k:**
    - n, m and Φ (returned to `batch_driver.py`, so `*_resultants.csv`, `Phi_FEniCS` and `force_error.csv` are in real units);
    - the stress fields S11–S23, Mises and inv_1 in `<case>_ref.vtu` and `<case>_deformed.vtu`, which `compare.py` reads;
    - the same columns in `<case>_section_X0.csv`.
  - Strains (E11–E33) and the displacement U are unchanged.
  - Each rescaled VTU gets a field-data marker, `stress_scale_applied = k`.
- **Exact:** the material is linear elastic, so stress ∝ E at a given strain.
- **Not rescaled:** FEniCS's own PNG plots (`*_Stress_centerline.png`) are drawn inside `warping_core.postprocess` before the rescale, so they show stresses at E = 1. `FENICS_PLOTS` is now **False** in `config.py` (backup in `_backup_fenicsplots_20261007/`), so they are no longer made. `contours.png` and `centerline_overlay.png` from `compare.py` use the rescaled VTU and are correct.
- **Not tested here:** pyvista and dolfin aren't available in this environment. In the pilot, check that FEniCS Mises ≈ 3D Mises in `centerline_overlay.png` (not 1.2e5 smaller).
- Backup in `_backup_fenicsE1_20261007/`.

### 2.16 Shear strains passed to FEniCS made consistent with the section (7 Oct)
- **Found in the pilot** (`pilot_table_1x1`, 5 cases): FEniCS showed in-plane stresses in the section frame (`transverse_ratio_fen`) of 0.009–0.086, growing with |V|. The 3D model showed about 0.001.
  - Mises error ran 6.4% (low shear) to 11.3% (high shear), and S13 error on the common scale 1–10%.
  - The 1D beam gives SF/SE = 3.14e4 ≈ 0.68·G·A (Abaqus's softened default). The FEniCS section's own shear stiffness is about 0.85·G·A.
  - FEniCS is given both the strains and the forces, so the ~25% conflict is absorbed by the Lagrange multipliers as spurious stress.
  - Bending (SM/SK = EI) and twist (SM3/SK3 ≈ GJ) are consistent.
  - The effect does not depend on section size: the 0.5×0.5 pilot gives the same numbers.
- **Fix** (`batch_driver.py`, `config.py`): before the FEniCS solve, the transverse shear strains are recomputed from the 1D shear forces as SE2 = SF2/K and SE3 = SF3/K, with K = `SHEAR_K_FACTOR`·G·A.
  - The default is Cowper's factor for a solid rectangle, 10(1+ν)/(12+11ν) = 0.8497.
  - The axial strain, curvatures and twist are untouched.
  - In the pilot, Abaqus's shear strains were exactly 1.25× the corrected ones.
  - Each case records `shear_strain_correction.csv` (copied into `summary.csv`). `CORRECT_SHEAR_STRAINS = False` restores the old behaviour.
- **To verify:** rerun FEniCS on the pilot (`--redo-fenics` keeps the ODBs). `transverse_ratio_fen` should fall to about the 3D level, and the S13/S23 errors should drop.
  - If it stays above about 0.01 for the high-shear cases, tune `SHEAR_K_FACTOR`. A single FEniCS solve with only a unit shear force would give the section's exact value.
- **Backup:** `_backup_shearfix_20261007/` (batch_driver.py, config.py). `strain_checks.py` only had one file name added to its diagnostics list.

### 2.17 E chosen per section so E·I is fixed (new table), and the shear-strain correction result (7 Oct)
- **E = E_REF / h⁴** with E_REF = 1.2e5, giving E·I = 1e4, **F0 = 1 and M0 = 100 at every size**:

  | Run | h | E |
  |---|---|---|
  | A | 1 | 1.2e5 |
  | B | 0.5 | 1.92e6 |
  | C | 0.25 | 3.072e7 (templates and mesh not made yet) |

  - One loads file now serves all sizes: `loads_table.csv` (Vx = v_x, Vy = v_y, Mx = 100·m; 1300 rows; `loads_table_1x1.csv` and `loads_table_half.csv` are identical copies).
  - All 1300 rows pass the gates at h = 0.5.
  - σ0 = E·c/L = 600 (h = 1) and 4800 (h = 0.5).
  - The non-dimensional problem is exactly the same as before; only the units change.
- **New templates:** `1d_fx_fy_mx_half_ar_EI1e4.inp` and `3d_fx_fy_mx_half_ar_finer_EI1e4.inp` (E = 1.92e6). h = 1 keeps the `_E1p2e5` templates.
- **`config.py`:** E is set from `BATCH_SECTION_SIDE` after the section is known, so G, F0, M0 and σ0 follow automatically. The FEniCS E = 1 rescale uses the per-run E. `make_table_loads.py` was updated to match.
- **Backups:** `_backup_EIconst_20261007/` (config.py, make_table_loads.py).
- **Re-run the 0.5×0.5 pilot** into a new folder. `pilot_table_half` was made with E = 1.2e5 and the old loads.
- **Shear-strain correction (§2.16) result:** the `--redo-fenics` pilot used the corrected strains (`SE*_abaqus_over_used` = 1.25 in every case), but `transverse_ratio_fen` and the errors were essentially unchanged.

  | Case | `transverse_ratio_fen` before → after | Mises error before → after |
  |---|---|---|
  | 0 | 0.086 → 0.086 | 9.8% → 9.8% |
  | 4 | 0.083 → 0.083 | 11.3% → 11.3% |

  - So the FEniCS stresses are set by the force targets, not the prescribed shear strains, and the 25% shear-stiffness mismatch is **not** the cause.
  - The spurious in-plane stress (equal S11 and S22, growing with |V|, absent in 3D) must come from how the section model handles the shear/flexure case. That needs a look inside `warping_core.py` (not accessible from this session).
  - The correction is left on: it is harmless and makes the inputs consistent.

### 2.18 Overlays for all components; warping_core findings (7 Oct)
- **`compare.py`** (backup in `_backup_overlays_20261007/`):
  - **`contours.png`:** all of S11, S22, S33, S12, S13, S23, Mises, I1 and J3.
  - **`centerline_overlay.png`:** all 6 components plus Mises and I1, in global axes (2×4).
  - **New `centerline_overlay_section.png`:** the 6 components pulled back to the section frame (Rᵀ·σ·R), where beam theory says S11, S22 and S12 ≈ 0.
  - **`strain_overlay.png`:** all 6 strain components. The Abaqus shear strains are halved, because LE12/13/23 are engineering strains while warping_core's E12/13/23 are tensor (Green-Lagrange) strains. The pilot confirms the factor: E23 Abaqus/FEniCS = 2.15.
  - The σ0 label now reads E·c/L.
- **S11 = S22 artifact existed before these changes:** S11 and S22 are perfectly correlated (r = 1.00) along the FEniCS centre line in `cases_full` (24 Sep), `cases_half_ar_123`, `cases_half_ar_h4_finer` and `cases_quadfix`, growing with |V|.
- **It is not negligible.** Equal in-plane stress p lowers Mises by about p and adds 2p to I1. In the pilot, I1 error rises faster than Mises with |V| (case 0: Mises 9.8%, I1 16.2%; case 4: 11.3% and 17.0%), while the low-shear cases stay at about 6.5% for both.
- **warping_core `m_con`:** the active third component ∫x1·x2 dA = 0 does not stop in-plane rigid rotation of a square section.
  - Under a rotation θ, ∫x1'x2' = ½·(I22 − I11)·sin 2θ + ∫x1x2·cos 2θ. With I11 = I22 (square, circle, any section symmetric under a 90° turn) it stays 0 for every θ.
  - It only fixes the rotation when I11 ≠ I22.
  - It does constrain in-plane shear distortion (x1 → x1 + γ·x2 changes it by γ·I11).
  - The commented-out ∫(X1·u2 − X2·u1) dA = 0 is the general rigid-rotation constraint.

### 2.19 Strain-frame fix: no Rᵀ pull-back of Abaqus SE/SK (8 Oct)
- **Root cause of the S11 = S22 artifact, the "twist" problem and most of the apparent size effect.** Abaqus B31 SE/SK (and SF/SM) are already components in the beam-local director basis dᵢ = R·eᵢ, i.e. they *are* v0 = Rᵀv and k0 = Rᵀk of Arora, Kumar & Steinmann (2019). The driver applied `R_local.T` to them again, so FEniCS got strains rotated twice while the force/moment targets (n0, m0) were not rotated. The force constraints then forced the section to absorb the mismatch, which showed up as equal in-plane stresses S11 ≈ S22.
- **Change:** `batch_driver.py` and `recover_section.py` now use `sv_mat = sv_sp`, `sk_mat = sk_sp`. `R_local` is still used on the output side (`postprocess(R_local)`, `abq_n = R_local @ n_t_sp`, `abq_m = R_local @ m_t_sp`). `strain_checks.rotated_beam_components()` also returns the unrotated director-basis values now; the column names in `beam_1d_mid.csv` (`*_rot_global`) are kept for file compatibility. This changes `eps_1d_mid_post` (the small-strain gate value) slightly. `trial_noforce.py` now defaults to `--strain-frame local`. The original `FEniCSfiles/warping_core.py` is untouched; its FRAMES note (`sv_mat = R_local.T @ sv_spatial`) is wrong for Abaqus beam output and should be read with this note.
- **Trial evidence (pilot, local strains):**
  - `transverse_ratio` (FEniCS) 0.0019–0.0057, down from up to 0.086;
  - force multipliers α ≈ 1e-6, so the constraints are nearly inactive;
  - Mises 6.3–7.4% and I1 ≈ Mises; S13 (÷Szz) 0.2–2.4%;
  - without force constraints the moments still match to 1e-4·M0.
  - Small leftover mismatches: the strain-only shear stiffness is 0.64–0.84 G·A, and N3 is 0.09–0.45 from the Green–Lagrange κ² terms. With the force constraints on these are absorbed with negligible stress effect, so the constraints stay on.
- **Earlier conclusions to revise:** the "high-twist cases are bad" and "size effect = geometric nonlinearity" readings (README, §2.18) were mostly this bug. In `cases_quadfix` and `cases_half_ar_123`, an Rᵀ curvature-mixing metric correlates with the Mises error (r = 0.82–0.97) and the I1 error (r = 0.95–0.98).
- **All earlier FEniCS results are affected.** Re-run them with `--redo-fenics` (keeps the Abaqus .odb files).
- Backup: `_backup_strainframe_20261008/`.

### 2.20 Tests for the 1D-convergence, hourglass and mid-moment gates (8 Oct)
- **New `test_gates.py`.**
  - Part A runs synthetic pass/fail inputs for each gate, including Abaqus's 1e-2 time-average fallback. All pass.
  - `--scan DIR …` applies the 1D and hourglass gates to existing case folders and writes `gate_scan.csv` in each.
- **Real positives for the 1D gate:** in `cases_half_ar_h4_finer` (E = 1), 546 of 1000 cases hit the fallback and 411 would be discarded (worst: case 578 at 1.68).
- **Hourglass:** no real positives yet. Only the new runs write `solid_energy_history.csv`; the pilot gives 0.02–0.26%.
- **End-to-end driver test:** `batch_test/gate_test_1x1` holds copies of the pilot cases 0–4 (with their ODBs). `batch_test/loads_gate_test.csv` holds the 5 pilot loads, one row that fails only the mid-moment floor (m_mid = 0.003), and 10 spares (`loads_table` rows 5–14). Run it with thresholds that split the pilot cases (command in the chat of 8 Oct).

---

## 3. Run inventory (`batch_test/`)

| Folder | Section | Loads | E | 3D mesh | 1D shear K | Status / validity |
|---|---|---|---|---|---|---|
| `cases_full` | 1×1 | ±0.6 F0/M0 (random) | 1 | coarse (10:1) | 0.4167 | Older code (pre strain-screen sign fix, pre deformed-moment fix). 1D–3D moment comparison unreliable; stress errors fine |
| `cases_half`, `cases_3d_half`, `cases_1d_half_5` | 1×1 | ±0.6 | 1 | half mesh | 0.4167 | Older reference runs |
| `cases_half_ar_123` | 0.5 | h³ (±1.2 F0/M0) | 1 | coarse (20:1) | 0.1042 | 1D partly under-converged |
| `cases_half_ar_1d_finer_quad_123` | 0.5 | h³ | 1 | coarse (20:1) | — | 1000 cases; same caveats |
| `cases_half_ar_3d_finer_123` | 0.5 | h³ | 1 | finer (10:1) | — | 758 cases; 60 recoverable with `rerun_dropped.py`. **Do not resume with the current `config.py`** |
| `cases_half_ar_h4_finer` | 0.5 | h⁴ (±0.6), paired with `cases_full` | 1 | **coarse** (20:1), despite the name | default | Rotation test done; 1D under-converged in its tail |
| `cases_quadfix` | 1×1 | ±0.6 (random, screened at 1×1) | 1 | coarse (10:1) | 0.4167 | Quad FEniCS mesh; source of the option C loads |
| `pilot_table_1x1` | 1 | `loads_table` rows 0–4 | 1.2e5 | coarse (10:1) | default | 5 cases. Abaqus 7 Oct; **FEniCS re-solved 7 Oct 19:07 after the strain-frame fix** (results in §7) |
| `pilot_table_half` | 0.5 | `loads_table` rows 0–4 (made with the old E = 1.2e5) | 1.2e5 | finer | default | **Obsolete**: E should be 1.92e6; FEniCS not re-solved after the fix. Re-run into a new folder |
| `gate_test_1x1` | 1 | `loads_gate_test.csv` (5 pilot loads + 1 floor-failing row + 10 spares) | 1.2e5 | coarse | default | Copy of the pilot, for the gate end-to-end test (§2.20). No `discard_counts.csv` yet, so that test **has not been run** |
| `cases_EI_1x1` | 1 | `loads_table.csv`, 1000 + 300 spares | 1.2e5 | coarse | default | **Full run in progress** (started 7 Oct 19:11, about 860 case folders by 8 Oct 05:12). Uses the strain-frame fix |
| `cases_EI_half` *(to run)* | 0.5 | `loads_table.csv` | 1.92e6 | finer | default | Planned, after `cases_EI_1x1` finishes |
| `cases_quadfix_E1e6` *(superseded)* | 1×1 | `cases_quadfix` loads × 1e6 | 1e6 | coarse (10:1) | default | Planned |
| `cases_half_ar_E1e6_finer` *(superseded)* | 0.5 | `cases_quadfix` loads × 1e6 × 0.5⁴ | 1e6 | finer (10:1) | default | Planned |

---

## 4. Commands (current)

Run one after the other: each uses 30 of 38 licence tokens, and the two runs need roughly 70 GB of disk. Rejected cases are replaced from the 300 spares first, then by random top-up draws. If a run stops, repeat the same command; finished cases are reloaded.

```
cd /mnt/c/Users/macLab/Desktop/FEniCSfiles/batch
BATCH_SECTION_SIDE=1 python3 batch_driver.py --cases 1000 \
  --loads-csv /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/loads_table.csv \
  --out-dir /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_EI_1x1
BATCH_SECTION_SIDE=0.5 python3 batch_driver.py --cases 1000 \
  --loads-csv /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/loads_table.csv \
  --out-dir /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_EI_half
```

Other commands:
- **Re-solve FEniCS only** (keeps the Abaqus ODBs): add `--redo-fenics` to a run command with the same `--out-dir`. This is how old folders get the strain-frame fix.
- **Gate unit tests:** `python3 test_gates.py` (synthetic, no Abaqus), and `python3 test_gates.py --scan DIR …` to apply the 1D and hourglass gates to existing folders.
- **Gate end-to-end test** (deliberately tight limits so that some pilot cases fail):
  ```
  BATCH_SECTION_SIDE=1 python3 batch_driver.py --cases 5 --redo-fenics \
    --max-1d-residual 1.5e-4 --max-hourglass 0.0015 \
    --loads-csv .../batch_test/loads_gate_test.csv --out-dir .../batch_test/gate_test_1x1
  ```
  Expected in `discard_counts.csv`: `pre_mid_moment_floor` = 1, `unconverged_1d` ≥ 2 (cases 0 and 2), `hourglass_energy` ≥ 1 (case 4), and 5 accepted cases at the end. About three spares need Abaqus. Do not use these limits for real runs.

**After the full runs, check that:**
- `discard_counts.csv` shows no 1D-fallback cases (`tavg_fallback_1d`) among the accepted ones;
- `transverse_ratio_fen` is below about 0.01 for essentially all cases;
- Mises and I1 errors are close to each other, and the 0.5 and 1 errors agree case by case (pair by non-dimensional load);
- the tail of the error histogram: if it is still long, look at the largest cases.

---

## 5. Open items

- **Pilot not yet re-checked at 0.5×0.5** with E = 1.92e6 and the strain-frame fix. Compare with `pilot_table_1x1` by non-dimensional load. The full `cases_EI_half` run will show this.
- **Gate end-to-end test not run** (§2.20). The hourglass gate has no real positive case yet.
- **Re-solve old folders** (`cases_quadfix`, `cases_half_ar_123`, `cases_half_ar_h4_finer`) with `--redo-fenics` to see how much of the earlier error tail and size effect was the Rᵀ bug.
- **Shear-flexure in the section model:** the strain-only shear stiffness is 0.64–0.84·G·A (it varies by case) and N3 is 0.09–0.45 from the Green–Lagrange κ² terms. Both are absorbed by the force constraints with negligible effect on the stresses (§2.19).
- **SF2/SE2 = 0.25 anomaly:** direction-2 shear stiffness ignores K (§2.2). Unexplained; test with `*Preprint, model=YES` or `*Section Print`.
- **S13/S23 gap:** in h⁴ vs `cases_full`, S13 and S23 errors stayed 4–8% higher even in well-converged cases. Re-check after the fix: in the new pilot the S13/S23 errors are 0.2–2.8% of ‖S33‖.
- **3D first increment:** align the 3D `*Static` first increment between the 1×1 and 0.5×0.5 templates (1.0 vs 0.01). Noted only; results should be unaffected.
- **`warping_core.py` FRAMES note** (`sv_mat = R_local.T @ sv_spatial`) is wrong for Abaqus beam output. The original file is untouched; this file is the record.
- **Older README conclusions** ("twist" and "size effect = geometric nonlinearity", `README.md` and the outputs README) predate §2.19 and are partly wrong; see §7.

## 6. Backups

| Folder (in `batch/`) | Contents |
|---|---|
| `_backup_timeout_20261006/` | `abaqus_io.py` before the 2400 s limit |
| `_backup_h4_20261006/` | `config.py` before h⁴ scaling |
| `_backup_Escale_20261007/` | `config.py` before E = 1e6 / `BATCH_SECTION_SIDE` |
| `_backup_checks_20261007/` | `abaqus_io.py`, `batch_driver.py`, `compare.py`, `strain_checks.py`, `rebuild_summary.py` before the convergence / floor / section-frame checks |
| `_backup_stressscale_20261007/` | `config.py`, `compare.py`, `batch_driver.py`, `rebuild_summary.py` before σ0 = E·c/L and `relL2_Szz` |
| `_backup_strainframe_20261008/` | `batch_driver.py`, `strain_checks.py`, `recover_section.py`, `trial_noforce.py`, this file, before the strain-frame fix (§2.19); `PARAMETERS_AND_CHANGES_before_record.md` is this file before the 8 Oct rewrite of §1, §3–§5 |
| `_backup_20260925/`, `_pre_fix_backup/` | Earlier pipeline fixes (see `README.md`) |

---

## 7. Observations (what we found, in order)

1. **Tails in the old histograms** (`cases_half_ar_1d_finer_quad_123`, `cases_half_ar_3d_finer_123`): the worst cases (460, 47, 787, 35) were mostly near-zero mid-span moment (relative errors blow up), under-converged 1D runs (Abaqus 1e-2 zero-force fallback at E = 1) and high shear. Dropped 3D cases were timeouts at 900 s. Fixes: §2.3–§2.4, E scaling §2.6–2.7, floor §2.13.
2. **Size effect** (0.5 section worse than 1.0): first read as geometric nonlinearity from larger rotations (§2.5). **Revised:** most of it came from the strain-frame bug (obs. 6).
3. **Loads must be non-dimensional** (F0 = EI/L², M0 = EI/L) so both sizes see the same bending strain and rotation. E is chosen per size so that E·I is fixed (§2.17).
4. **Gates exactly as the table says** (strain, shear/bending ratio, mid-span moment ≥ 0.1·M0), plus 1D-residual and hourglass gates after Abaqus. 1300 loads were made, 31.7% pass rate (§2.14).
5. **S11 ≈ S22 artifact in FEniCS:** perfectly correlated along the centre line, growing with |V|, absent in 3D (`transverse_ratio_fen` up to 0.086 vs 0.001 in 3D). It was present in every earlier folder (`cases_full`, `cases_half_ar_123`, `cases_half_ar_h4_finer`, `cases_quadfix`). The shear-stiffness correction (§2.16) did not change it (§2.17), and neither did the angle-constraint variants of `m_con` (§2.18).
6. **Root cause (§2.19):** Abaqus B31 SE/SK are already v0/k0 (director basis) of Arora, Kumar & Steinmann (2019). The driver applied `R_local.T` again, a double rotation, while the force and moment targets were not rotated, so the section solve had to absorb the conflict as spurious stress. Check on the 2019 paper: n0, v0 etc. are the undeformed-configuration material quantities; R is only needed on output.
7. **Effect of the fix, trial and then `pilot_table_1x1` re-solved with the real driver** (7 Oct, 19:07):

   | Quantity | Before | After |
   |---|---|---|
   | `transverse_ratio_fen` | 0.009–0.086 | 0.0019–0.0057 (3D: 0.0007–0.0017) |
   | Mises error (÷‖Mises₃D‖) | 6.4–11.3% | 6.3–7.4% |
   | I1 error | up to 17% | 5.8–6.6%, close to Mises |
   | S33 error (÷‖S33₃D‖) | — | 6.3–6.9% |
   | S13, S23 errors (÷‖S33₃D‖) | 1–10% | S13 0.2–2.4%, S23 1.6–2.8% |
   | S11, S22 errors (÷‖S33₃D‖) | — | 0.25–0.74% |
   | Force multipliers α (trial) | — | about 1e-6 |

   Without force constraints the moments matched to 1e-4·M0 (trial). Keeping the constraints costs nothing, so they stay on.
8. **The Rᵀ curvature-mixing metric explained the old errors:** in `cases_quadfix` and `cases_half_ar_123` it correlates with the Mises error (r = 0.82–0.97) and the I1 error (r = 0.95–0.98).
9. **Gate checks on old data:** `test_gates.py --scan` on `cases_half_ar_h4_finer` (E = 1): 546 of 1000 cases hit the 1e-2 fallback and 411 would be discarded by the 1% residual gate (worst: case 578 at 1.68). The synthetic tests for all three gates pass. The pilot energies are 0.02–0.26% hourglass.
10. **Load file does not change with the fix.** The strain-frame fix touches only how FEniCS reads Abaqus strains. The pre-run gates (loads only), the post-run admission checks (use raw Abaqus values) and the 1D/hourglass gates do not use the rotated strains, so the same cases pass or fail. The old 1000-case folders are re-solved with `--redo-fenics`, not re-drawn.
11. **Run state (8 Oct 05:13):** `cases_EI_1x1` is running with the fix (about 860 case folders). `cases_EI_half` and the gate end-to-end test are not run. `pilot_table_half` is obsolete (old E).

---

## 8. Results of `cases_EI_1x1` (finished 8 Oct 07:03; analysed 8 Oct)

**Run:** 1000 accepted, 0 discards at any gate, except one `abaqus_failed` (case 469). Its 1D and 3D jobs both finished successfully but nothing was post-processed; recover it with `rerun_dropped.py`. Case 1000, the first spare, took its place. Largest 1D residual 0.65% of the load; hourglass ALLAE/ALLIE 0.012–0.53%; `transverse_ratio_fen` 0.0007–0.0073 (3D: 0.0002–0.0022).

**Error distribution, 1 × 1 section (percent, 1000 cases):**

| Measure | min | median | p90 | p99 | max |
|---|---|---|---|---|---|
| Mises (÷‖Mises₃D‖) | 5.3 | 6.4 | 7.1 | 7.8 | 8.4 |
| I1 | 5.4 | 6.2 | 6.6 | 7.4 | 8.2 |
| S33 (÷‖S33₃D‖) | 5.4 | 6.4 | 6.9 | 7.5 | 8.3 |
| J3 (÷‖J3₃D‖) | 2.6 | 3.5 | 3.7 | 4.2 | 4.8 |
| S13 (÷‖S33₃D‖) | 0.01 | 0.9 | 2.2 | 2.8 | 3.2 |
| S23 (÷‖S33₃D‖) | 0.1 | 1.6 | 2.6 | 3.5 | 4.0 |
| Energy | 0.01 | 1.4 | 3.0 | 4.4 | 5.9 |

There is a floor, not a tail: the earlier tails (Mises to 18–21%, S13 to 170%, energy to 27%) are gone.

**What drives the spread:**
- **Floor (about 5.5–6.5%):** the 3D reference. Along the section centre line the 3D S33 slope is 0.962 of FEniCS (median; 758 cases with strong y-bending). At the outer fibre the 3D stress is 0.905 of FEniCS. The 3D profile departs from a straight line by 6% of its peak (FEniCS: 0.4%), which is the staircase of single-point C3D8R elements. The gradient shortfall rank-correlates with the S33 error at −0.88. In 41 clean bending cases, FEniCS is at 0.988 and 3D at 0.949 of the beam-theory slope M/I. This is consistent with the 23 Sep report (about 7% bending-gradient shortfall of C3D8R, 10 elements across). **Not yet proven:** the test is a refined or C3D8I/C3D20R 3D run.
- **Spread above the floor (Mises 5.3 → 8.4%):** grows with |Mx| (Spearman 0.81) and with the twist ratio |UR3|/√(UR1²+UR2²) (0.78). Regression on |Mx|, 1/lever arm and twist gives R² = 0.78 for Mises. Top Mises/I1/J3 cases (769, 406, 141, 887, 450, 83): large tip moment, but the moment at mid-span is near the 0.1·M0 floor (moment gradient large relative to moment: lever arm M_mid/V about 20 against a median of 95). The section solver assumes strain varies slowly along the beam, so this is where that assumption is weakest.
- **S13/S23 tails (2.5–3.2% and 3.0–4.0%):** the high-shear cases (24, 27, 9, 438, 425, 245, 704, 305): |V| about 0.8, twist ratio 0.4–0.6, bending dominance about 4. Shear errors correlate with |Vx| (0.81), bending dominance (−0.92) and hourglass energy (0.92: both just track the shear share).
- **Energy tail (4–6%):** small |Mx|, small moment at mid-span (cases 804, 382, 754, 687, 575, 213, 763, 961). Energy is then small and a fixed absolute error shows up as a larger share. Correlation with |Mx| is −0.84.
- **J3:** relL2_resultant_J3 is 2.6–4.8%. It is normalised by ‖Mises₃D‖³, unlike the 23 Sep own-norm figure (21%), so the two are not comparable. The cubic-amplification explanation of 23 Sep has not been re-tested on these data.

**Open question 1 (shear stiffness anomaly): not present any more.** In all 995 cases with shear, SF2/SE2 = SF3/SE3 = 31384.6 = 0.6800·G·A in both directions. The 0.25-versus-K/2 pattern belonged to the old templates with an explicit `*Transverse Shear Stiffness` line; the current templates use Abaqus's default and show no anomaly. The cause in the old templates was never found, but it no longer matters. The shear-strain correction divides by 0.8497·G·A, so Abaqus/used = 1.2495 in every case.

**Open question 2 (S13/S23 gap): resolved for h = 1.** Errors are 0.01–4% of the axial stress. Whether a gap between sizes remains needs `cases_EI_half`.

**Still to do:** `cases_EI_half`; recover case 469; the refined-3D test; `--redo-fenics` on older folders; the gate end-to-end test.


### 8.1 Edge-node test of the 6% floor (8 Oct, no new Abaqus runs)
- **Reasoning:** C3D8R gives one stress per brick; a node at the surface (y = ±0.5) takes the value at the nearest brick centre (y = ±0.45), 90% of the true value. Measured edge ratio 3D/FEniCS = 0.905. Arithmetic for 11 node rows with the two edge rows 10% low gives a 6.7% rms error (observed median 6.4%) and a slope 4.5% low (observed about 4%). With 20 bricks across the prediction is about 2.5%.
- **Data test:** S33 error along the centre line (x1 = 0), 1000 cases: 6.62% median with all 11 points; **1.00%** without the two edge points. Full-section value in summary.csv: 6.40%.
- **Twist trend:** with all points the error rises with twist (6.1 → 7.4% from the lowest to the highest quarter); without the edge points it does not (1.6 → 0.5%). So the rise with twist comes from the edge nodes.
- **Total rotation does not matter:** sqrt(UR1² + UR2² + UR3²) has R² = 0.007 against the Mises error. The earlier "large rotations" suspicion is not supported. Twist ratio and |Mx| are the predictors.
- **Still to do (needs Abaqus):** same 10 × 10 mesh with C3D8I elements (change `*Element, type=C3D8R`), and 20 × 20 bricks across. Expected floor about 1–2% (C3D8I) and about 2.5% (20 bricks).

---

## 9. Choosing the 3D element type; element-comparison run (8 Oct)

**Changes** (backups in `_backup_element3d_20261008/`):
- `config.ELEMENT_3D` (environment variable `BATCH_ELEMENT_3D`) and `batch_driver.py --element-3d TYPE`: replace the `*Element, type=` line of the 3D template when each `.inp` is written (`abaqus_io.set_solid_element`). Default is unchanged (the template's C3D8R).
- Only 8-node bricks are accepted: `C3D8R`, `C3D8`, `C3D8I`, `C3D8H`, `C3D8RH` and the like. The template's connectivity is first-order, so quadratic bricks (C3D20R) and tetrahedra are refused with an explanation; they need a new mesh and a matching FEniCS section mesh.
- Each run appends a line to `<out-dir>/run_settings.txt` (section, E, 3D template, element type, loads file).
- New `make_subset_loads.py` (pick cases from a finished batch) and `compare_element_runs.py` (pair the replay with the original and print errors side by side).
- Notes: C3D8 and C3D8I have no hourglass modes, so ALLAE = 0 and the hourglass gate always passes. Run the new element into a **new** `--out-dir`.

**Run (C3D8I test; 30 cases: 5 named, 11 more of the worst Mises, 14 of the best):**
```
cd /mnt/c/Users/macLab/Desktop/FEniCSfiles/batch
python3 make_subset_loads.py --src /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_EI_1x1 \
  --out /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/loads_subset_elem.csv \
  --cases 769,406,24,438,804 --worst 15 --best 15
BATCH_SECTION_SIDE=1 python3 batch_driver.py --cases 30 --element-3d C3D8I \
  --loads-csv /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/loads_subset_elem.csv \
  --out-dir /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_EI_1x1_C3D8I
python3 compare_element_runs.py --base /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_EI_1x1 \
  --new /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_EI_1x1_C3D8I \
  --loads-csv /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/loads_subset_elem.csv
```
Other element: change `--element-3d` (and the out-dir name), e.g. `C3D8` (full integration) or `C3D8H`. Prediction for C3D8I: the floor falls from about 6% to 1–2% and the rise with twist disappears. If the floor stays near 6%, the explanation in §8.1 is wrong.


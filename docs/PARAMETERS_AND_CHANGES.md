# Beam 1D / solid 3D / FEniCS section study: parameters and change log

Last updated **7 Oct 2026**. This is a running record of what each run used, what changed between runs, and why.

Sources:
- this chat's analysis of the case folders;
- `batch/README.md` (pipeline, checks, frames);
- `batch/Claude outputs/beam_comparison_summary.tex` (an earlier chat on the shear-stiffness and section-halving bugs).

For the pipeline itself (checks, frames, error measures), see `batch/README.md`. This file covers only parameters and the reasons for changing them.

---

## 1. Current set-up (for the paired runs started 7 Oct)

| Parameter | Value | Where | Why |
|---|---|---|---|
| Length L | 100 | `config.py` | Unchanged since the start |
| Mid-span station Z_MID | 50 | `config.py` | Where all errors are evaluated |
| Section | square, side 1 or 0.5, chosen with the `BATCH_SECTION_SIDE` environment variable | `config.py` | One config file serves both sizes; no hand edits between runs |
| E, ν | **1e6**, 0.3 | `config.py` and every `*_E1e6.inp` | Was E = 1. See §2.7 |
| Load ranges | V: ±0.6·F0; Mx: ±0.6·M0 (F0 = EI/L², M0 = EI/L) | `config.py` (`LOAD_ND = 0.6`) | Equal non-dimensional loads give equal rotations at both sizes. See §2.5 |
| Loads actually used | Replayed from `cases_quadfix` (option C) | `batch_test/loads_E1e6_*_from_quadfix.csv` | Paired case by case; all 1000 pass the strain screen at both sizes. See §2.8 |
| 1D model | B31, 100 elements (length 1), NLGEOM, one load increment (`*Static 1., 1., 1e-05, 1.`) | `1d_fx_fy_mx_E1e6.inp` (1×1), `1d_fx_fy_mx_half_ar_E1e6.inp` (0.5×0.5) | — |
| 1D transverse shear stiffness | **Abaqus default** (no `*Transverse Shear Stiffness` line) | both 1D templates | See §2.2 and §2.9 |
| 3D model | C3D8R, 10×10 elements across the section, NLGEOM, tip shears as non-follower tractions, Mx through a distributing coupling | `3d_fx_fy_mx_E1e6.inp` (1×1, 100 elements along L) and `3d_fx_fy_mx_half_ar_finer_E1e6.inp` (0.5×0.5, 200 elements along L) | Same 10:1 element aspect ratio at both sizes |
| 3D first increment | 1.0 (1×1) vs 0.01 (0.5×0.5) | templates | Not intentional. Abaqus cuts back as needed, so results should be unaffected; noted for consistency |
| FEniCS section mesh | 10×10 quads (121 nodes) | `slice_50_quad_native.xdmf` (1×1), `slice_ar_50_quad.xdmf` (0.5×0.5) | The "quadfix" quad mesh |
| Strain screen (pre-run) | root principal strain ≤ 5e-3 | `STRAIN_THRESHOLD` | Study limit |
| Admission (after Abaqus, before FEniCS) | stress-ratio bending dominance ≥ 1; max 1D strain < 5e-3 | `strain_admit` | Unchanged |
| Abaqus wait limits | 1D: 300 s; 3D: **2400 s** | `abaqus_io.py` | Was 900 s. See §2.4 |
| Parallel cases | 3 | `MAX_PARALLEL` | 3 × 2 jobs × 5 tokens = 30 of 38 licence tokens |

Non-dimensional scales at E = 1e6:

| | 1×1 | 0.5×0.5 |
|---|---|---|
| F0 = EI/L² | 8.333 | 0.5208 |
| M0 = EI/L | 833.3 | 52.08 |
| V range | ±5.0 | ±0.3125 |
| Mx range | ±500 | ±31.25 |

---

## 2. Change log: what changed and why

### 2.1 Section halving and the h³ load scaling (late Sep)
- The 0.5×0.5 section (`half_ar`) was introduced, and loads were scaled by side³ so the **bending strain** distribution matched 1×1. Pass rate was 62.4% at both sizes.
- **Side effect, noted at the time:** at equal strain, curvature and rotations scale as 1/side, so the 0.5 section rotates about twice as much. This turned out to drive the size effect (§2.5).
- `abaqus_io.write_inp` divides the 3D traction by A, because TRVEC is force per area. It was only correct before because A = 1.

### 2.2 Transverse shear stiffness (28–30 Sep, earlier chat)
- **Stale 0.5×0.5 value:** the 0.5 template kept the 1×1 value K = 0.4167. It was corrected to 0.4167/4 = 0.1042 in test runs (`cases_half_ar_seed1`). The later 0.5 batches (`cases_half_ar_123` onward) have no explicit K, so they use Abaqus's default. All 1×1 batches up to `cases_quadfix` use 0.4167.
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

---

## 3. Run inventory (`batch_test/`)

| Folder | Section | Loads | E | 3D mesh | 1D shear K | Status / validity |
|---|---|---|---|---|---|---|
| `cases_full` | 1×1 | ±0.6 F0/M0 (random) | 1 | coarse (10:1) | 0.4167 | Older code (pre strain-screen sign fix, pre deformed-moment fix). 1D–3D moment comparison unreliable; stress errors fine |
| `cases_half`, `cases_3d_half` | 1×1 | ±0.6 | 1 | 50 elements along L (20:1) | 0.4167 | Older reference runs |
| `cases_1d_half_5` | 1×1 | ±0.6 | 1 | coarse (10:1) | 0.4167 | Older reference run |
| `cases_half_ar_123` | 0.5 | h³ (±1.2 F0/M0) | 1 | coarse (20:1) | default | 1D partly under-converged |
| `cases_half_ar_1d_finer_quad_123` | 0.5 | h³ | 1 | coarse (20:1) | default | 1000 cases; same caveats |
| `cases_half_ar_3d_finer_123` | 0.5 | h³ | 1 | finer (10:1) | default | 758 cases; 60 recoverable with `rerun_dropped.py`. **Do not resume with the current `config.py`** |
| `cases_half_ar_h4_finer` | 0.5 | h⁴ (±0.6), paired with `cases_full` | 1 | **coarse** (20:1), despite the name | default | Rotation test done; 1D under-converged in its tail |
| `cases_quadfix` | 1×1 | ±0.6 (random, screened at 1×1) | 1 | coarse (10:1) | 0.4167 | Quad FEniCS mesh; source of the option C loads |
| `cases_quadfix_E1e6` *(to run)* | 1×1 | `cases_quadfix` loads × 1e6 | 1e6 | coarse (10:1) | default | Planned |
| `cases_half_ar_E1e6_finer` *(to run)* | 0.5 | `cases_quadfix` loads × 1e6 × 0.5⁴ | 1e6 | finer (10:1) | default | Planned |

---

## 4. Commands for the paired runs

Run one after the other (each uses 30 of 38 licence tokens). About 104 GB is free and the two runs need roughly 70 GB. Try `--cases 5` with a throwaway `--out-dir` first.

```
cd /mnt/c/Users/macLab/Desktop/FEniCSfiles/batch
BATCH_SECTION_SIDE=1 python3 batch_driver.py --cases 1000 \
  --loads-csv /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/loads_E1e6_1x1_from_quadfix.csv \
  --out-dir /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_quadfix_E1e6
BATCH_SECTION_SIDE=0.5 python3 batch_driver.py --cases 1000 \
  --loads-csv /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/loads_E1e6_half_from_quadfix.csv \
  --out-dir /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_half_ar_E1e6_finer
```

Compare only cases 0–999. Rejected cases are replaced by random top-up draws numbered 1000 and up, which have no partner.

**After the runs, check that:**
- the 1D residual force is below 1% of the applied V (grep `LARGEST RESIDUAL FORCE` / `TIME AVG. FORCE` in `*_1d.msg`; no `1.000E-02` fallback);
- the 0.5×0.5 and 1×1 errors match case by case;
- the zero-moment cases (47, 469, 819) are reported separately.

---

## 5. Open items

- **SF2/SE2 = 0.25 anomaly:** direction-2 shear stiffness ignores K (§2.2). Unexplained; test with `*Preprint, model=YES` or `*Section Print`.
- **S13/S23 gap:** in h⁴ vs `cases_full`, S13 and S23 errors stay 4–8% (relative) higher even in well-converged cases. Unexplained; possibly the different shear stiffness settings. Re-check in the E1e6 pair.
- **Convergence check in the driver:** reject any case with 1D residual force > 1% of V. Proposed, not implemented.
- **3D first increment:** align the 3D `*Static` first increment between the 1×1 and 0.5×0.5 templates (1.0 vs 0.01).
- **3D-finer folder:** run `rerun_dropped.py` on `cases_half_ar_3d_finer_123` to recover the 60 dropped cases. Optional, since that folder is superseded.

## 6. Backups

| Folder (in `batch/`) | Contents |
|---|---|
| `_backup_timeout_20261006/` | `abaqus_io.py` before the 2400 s limit |
| `_backup_h4_20261006/` | `config.py` before h⁴ scaling |
| `_backup_Escale_20261007/` | `config.py` before E = 1e6 / `BATCH_SECTION_SIDE` |
| `_backup_20260925/`, `_pre_fix_backup/` | Earlier pipeline fixes (see `README.md`) |

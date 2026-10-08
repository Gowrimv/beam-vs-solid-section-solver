# 1D vs 3D beam comparison: error tails and missing cases

Analysis of two Abaqus/FEniCS batch folders (data read only, nothing in them was modified):

| Folder | What is refined | Case dirs | Accepted (`summary.csv`) |
|---|---|---|---|
| `cases_half_ar_1d_finer_quad_123` | 1D beam mesh | 1012 (0–1011) | 1000 |
| `cases_half_ar_3d_finer_123` | 3D solid mesh (ODB ~39 MB vs ~20 MB) | 835 (0–834) | 758 |

All error values are the rel-L2 errors in `summary.csv`, evaluated on the mid-span section (Z = 50).

---

## 1. What causes the long tails in the error histograms

The tails are not random. The same cases sit in the tail in both folders, so the tails come from the load cases, not from either mesh. Three load regimes explain almost all of them.

### 1a. High twist + weak bending → tails in von Mises, I1, S33 (and part of S23)

Two quantities predict these errors almost perfectly (Spearman rank correlation, 1D folder; 3D folder within ±0.03):

| Predictor | von Mises | I1 | S33 | S23 |
|---|---|---|---|---|
| `bending_dominance_mid_abaqus` | **−0.97** | −0.96 | −0.83 | −0.79 |
| twist ratio = \|UR3\| / √(UR1² + UR2²) | **+0.93** | +0.89 | +0.82 | +0.73 |

Median error vs twist ratio (1D folder):

| Twist ratio | ≤0.05 | 0.05–0.1 | 0.1–0.2 | 0.2–0.3 | 0.3–0.5 | >0.5 |
|---|---|---|---|---|---|---|
| von Mises | 6.4% | 6.9% | 8.1% | 10.4% | 13.6% | 16.5% |
| I1 | 6.5% | 8.3% | 11.6% | 17.5% | 25.4% | 28.4% |

- Of the cases in the top 10% of at least three of {von Mises, I1, S33, S23}, **all 63** (1D folder) and **52 of 54** (3D folder) have `bending_dominance_mid_abaqus` < 6 (the median over all cases is ~13).
- 56 of 63 (1D) and 44 of 54 (3D) have twist ratio > 0.3.

**No torque is applied.** The 1D beam instance is rotated onto the Z axis, so the loads are two tip shears (Vx, Vy) and a **bending moment Mx about X**. The twist UR3 therefore comes only from large-rotation (geometrically nonlinear) effects:

- **How rotations combine:** bending about X and Y at the same time gives the rotation vector a Z component, roughly ∝ UR1·UR2. UR3 tracks UR1·UR2/2 with correlation 0.88 (1×1) and 0.80 (0.5×0.5).
- **Loads that keep a fixed direction:** `follower=NO` loads and the global Mx gain a component about the deformed beam axis once the tip rotates. That is a real induced torque (the small Mz in `force_moment_1d3d_error.csv`).

Interpretation: twist is a direct measure of how nonlinear a case is. The FEniCS section model is linear (small rotation) and has no such twist, so the error grows with it. A large twist also goes with a weak bending dominance, because the induced torque adds shear stress.

### 1b. Large moment relative to shear → the S13 tail

S13 has a different driver: the moment-to-shear ratio |Mx| / √(Vx² + Vy²) (Spearman +0.93). Mx is a bending moment, so this is a lever-arm length, not a torque. Tail cases (S13 > ~120%) have a median ratio of ~130, against ~48 for the rest. A large ratio means the shear force, and hence S13, is small in absolute terms, so the relative error is large. Worst cases: 560, 769, 29, 696, 673, 424, 280.

### 1c. Zero bending moment at mid-span → strain-energy outliers (and high von Mises/S33/S23)

These 11 cases are **identical in both folders**: **12, 35, 47, 113, 413, 423, 460, 547, 568, 643, 787**.

In these cases the bending moment passes through zero almost exactly at Z = 50, which is where all errors are evaluated:

- Mid-span bending moment / (V × 50) is 0.04–0.18 for these cases, against a median of ~1.3 over all cases.
- Case 460: Mx goes from 3.0e-5 at the root to −2.7e-5 at the tip, and is only ~1e-6 at Z = 50.
- The section strain energy is ~200× smaller than in a typical case (Φ ≈ 1.7e-10 against 4e-8 for case 1).
- What energy remains is mostly from transverse shear, which is where the 1D beam model is weakest. A small absolute miss therefore becomes a 14–27% relative error.

| Case | Energy err (1D folder) | Energy err (3D folder) | von Mises err (3D folder) |
|---|---|---|---|
| 460 | 20.2% | 27.3% | 18.1% |
| 47 | 14.0% | 16.1% | 21.0% (worst in folder) |
| 787 | 2.2% | 24.7% | 15.3% |
| 35 | 8.8% | 13.1% | 6.6% |

Likely contributor (not confirmed): the energy slice is 2.0 thick in the 1D folder and 1.0 thick in the 3D folder. Near a moment zero crossing, averaging M² over the slice adds a term ∝ V²h²/12 that depends on slice thickness. This may be why case 787 differs so much between the folders.

### 1d. What is not explained

The S23 tail is mostly outside all three groups. Those cases have a *lower* moment/shear ratio (~40) than average. Their cause is not pinned down.

---

## 2. Do 1D and 3D forces/moments match in the outlier cases?

Yes, in absolute terms. Comparison is from `force_moment_1d3d_error.csv`.

- **Vx, Vy** match within 1.7% at every Z, and **My** within 0.5% at mid-span. This is the same as typical cases.
- **Mx** has a nearly constant absolute offset along the whole beam. For case 47 it is ~1.0e-7 (1D folder) and ~2.1e-7 (3D folder), which is 0.5–1% at the root. At Z = 50, where Mx ≈ 0, the same offset appears as a 9–31% difference (9–21% of the total resultant moment). The moment doesn't really disagree more here; it is divided by a much smaller number.
- **Axial force N** in case 787 (1D folder) and case 35 differs by 15–27% of the resultant at mid-span. N should be ≈ 0, and 5% of all cases show a mismatch this large, so it is not specific to the outliers.

Typical mid-span mismatch as % of resultant (median / 95th percentile):

| | N | Vx | Vy | Mx | My | Mz |
|---|---|---|---|---|---|---|
| 1D folder | 0.15 / 10.4 | 0.02 / 1.7 | 0.03 / 2.5 | 0.09 / 1.0 | 0.12 / 0.42 | 0.14 / 0.41 |
| 3D folder | 0.24 / 4.8 | 0.02 / 0.54 | 0.03 / 1.5 | 0.30 / 2.1 | 0.35 / 0.89 | 0.16 / 0.52 |

---

## 3. Why the 3D-finer folder has 758 cases instead of 1000

**These cases did not fail to converge.** Abaqus finished every case it actually ran to completion. The 242 missing cases break down as follows:

| Reason | Cases | Count |
|---|---|---|
| Never run: the batch stopped at case 834 (the 1D folder ran to 1011) | 835–1011 | 176 |
| Abaqus completed, but the run was too slow and the case was dropped (see below) | 30, 56, 58, 73, 86, 91, 93, 123, 127, 151, 156, 158, 159, 161, 169, 180, 188, 189, 192, 195, 218, 224, 225, 268, 311, 331, 390, 403, 408, 420, 447, 448, 450, 462, 476, 484, 486, 497, 498, 509, 510, 512, 519, 523, 526, 527, 528, 531, 543, 562, 567, 570, 572, 573, 577, 581, 585, 605, 607, 782 | 60 |
| Batch interrupted: case 508 1D job got SIGTERM; 507 and 511 not post-processed | 507, 508, 511 | 3 |
| End of batch: ODB write failure ("Extend failed", typically disk full) at 822; 823–834 aborted or never started | 821–834 | 14 |

(The 1000 − 758 count does not reconcile exactly because the folders also differ in which cases were accepted: 7 cases (22, 730, 815–819) are accepted only in the 3D folder.)

### The 60 "too slow" cases

There is a clean split on Abaqus wall-clock time for the 3D job (from the `.dat` files):

| | n | Increments (median) | Wall-clock (s) |
|---|---|---|---|
| Accepted | 758 | 20 | 64–**850** (median 378) |
| Completed but dropped | 60 | 34 | **858**–2032 (median 965) |

No accepted case took more than 850 s, and every dropped case took more than 858 s. **This looks like a per-case timeout in the batch driver** (~15 min including pre/post overhead). The driver gives up waiting and skips post-processing, even though Abaqus later finishes successfully. The `.sta` files show `THE ANALYSIS HAS COMPLETED SUCCESSFULLY` and `.msg` has no errors. Only the post-processing outputs (CSVs, VTUs, PNGs) are missing, which is why these folders hold 13 files instead of ~35.

Why these cases passed in the 1D folder (coarser 3D mesh) but not here:

- They are the **hardest nonlinear cases**. In the 1D folder they also needed more increments (median 28 vs 16.5) and longer runs (median 387 s vs 180 s, max 718 s). They also have larger strains (eps_root_abaqus median 0.0035 vs 0.0026) and higher bending dominance.
- The finer 3D mesh roughly doubles the model size and adds more increments (34 vs 20). That pushes these cases from ~400 s to ~900–2000 s, past the timeout. Under the coarser mesh they stayed below it.
- 56 of the 60 were accepted in the 1D folder.

**Fix:** raise or remove the per-case timeout (≥ 2100 s would cover all 60). Alternatively, re-run only post-processing on these 60 folders, since their `_3d.odb` files are complete. Also free disk space before resuming from case 821.

---

## 4. Output files

| File | Contents |
|---|---|
| `histograms_grouped_1d.png`, `histograms_grouped_3d.png` | Original error histograms, stacked by group: zero mid-span moment (orange, also ▼ on the top edge), high twist + weak bending (blue), moment/shear > 100 (green), other (grey). Dashed line = 90th percentile. |
| `histograms_by_twist_1d.png`, `histograms_by_twist_3d.png` | Same histograms coloured by twist angle \|UR3\| (light → dark blue). The green outline (S13 panel only) marks moment/shear > 100; ▼ marks the zero-moment cases. |
| `size_check_error_vs_twist.png` | Median error vs twist angle (top) and vs twist ratio (bottom) for the 1×1 and 0.5×0.5 runs (§5). |
| `tail_cases.csv` | Every case in a group or in the top 10% of any error metric, both folders. Columns: group, `n_metrics_in_tail`, twist ratio, bending dominance, `moment_over_shear`, loads, all errors, per-metric p90 flags. |

Group definitions (applied in this order):

1. **Zero mid-span moment**: the 11 cases listed in §1c (mid-span moment / (V·50) < 0.2).
2. **High twist + weak bending**: twist ratio > 0.3 and `bending_dominance_mid_abaqus` < 6.
3. **High moment/shear**: |Mx| / |V| > 100 (Mx is the bending moment about X, so this ratio is a length).
4. **Other**: everything else.

---

## 5. Why the 0.5×0.5 section shows larger errors than 1×1 (same L = 100)

**Load scaling.** The 0.5 runs scaled the loads by h³ (`_LOAD_SCALE = SECTION_SIDE**3`), which keeps the strain equal across sizes. In the non-dimensional units F0 = EI/L² and M0 = EI/L, that is twice the load:

| | V range | Mx range |
|---|---|---|
| 1×1 (`cases_full`) | ±0.6 F0 | ±0.6 M0 |
| 0.5×0.5, h³ scaling | ±1.2 F0 | ±1.2 M0 |

Rotation ≈ M/M0 and strain ≈ (M/M0)·(c/L). At a fixed L, equal strain therefore means twice the rotation for the half-size section. Strain and rotation can only both match if L/h is also kept fixed (full geometric similarity).

**Evidence that the extra error comes from the larger rotations.** Median von Mises error by twist angle |UR3|:

| Twist (rad) | ≤0.01 | 0.01–0.02 | 0.02–0.04 | 0.04–0.06 | 0.06–0.1 | 0.1–0.15 | 0.15–0.3 |
|---|---|---|---|---|---|---|---|
| 1×1 `cases_full` | 6.1 | 6.5 | 7.0 | 7.9 | 8.9 | — | — |
| 1×1 `cases_half` | 5.8 | 6.3 | 6.9 | 7.7 | 8.9 | — | — |
| 0.5 `cases_half_ar_123` | 6.1 | 6.5 | 7.1 | 7.9 | 9.4 | 11.8 | 14.4 |
| 0.5 3D-finer | 6.1 | 6.4 | 7.1 | 7.9 | 9.2 | 11.5 | 14.0 |

- **Von Mises and S33:** the two sizes agree at the same twist angle. They do not agree at the same size-independent twist ratio. The 0.5 tail is simply the cases with twist > 0.1 rad, which the 1×1 runs never reach.
- **I1:** only partly lines up (about 10–12% vs 11–15% at 0.05–0.08 rad), so something else also contributes there. Candidates are the coarse 0.5 mesh (elements 20:1 instead of 10:1) and the shear-stiffness settings in the older 1×1 folders.
- **Caveat:** this compares two different sets of cases, not a controlled test.

**Controlled test (set up 6 Oct, to be run):**

- `config.py`: `_LOAD_SCALE = SECTION_SIDE**4`, giving V and Mx ranges of ±0.6 F0 and ±0.6 M0, the same as the 1×1 runs.
- `config.py`: `INP_3D` set to `3d_fx_fy_mx_half_ar_finer.inp`, which has 10:1 elements, the same shape as `cases_full`.
- Loads replayed from `batch_test/loads_h4_from_cases_full.csv`, which is the `cases_full` loads × 0.5⁴. New case *i* is the same non-dimensional problem as row *i* of `cases_full/summary.csv`. Checked: the maximum difference in V/F0 and M/M0 is 0.
- Bending rotations (and therefore the induced twist) match by construction. Strains are halved (∝ h/L), which is harmless for a linear-elastic material.

```
python3 batch_driver.py --cases 1000 \
  --loads-csv /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/loads_h4_from_cases_full.csv \
  --out-dir /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_half_ar_h4_finer
```

How to read the result:
- **Errors match `cases_full`:** the size effect comes from rotation (geometric nonlinearity).
- **A gap remains:** look at the shear-stiffness settings or the section mesh.

## 6. Changes made to the batch code (6 Oct)

| File | Change | Backup |
|---|---|---|
| `abaqus_io.py` | `ABAQUS_TIMEOUT_3D` 900 → 2400 s | `_backup_timeout_20261006/` |
| `rerun_dropped.py` (new) | Finishes cases whose Abaqus jobs completed but were never post-processed (extraction + FEniCS + compare only; no Abaqus). Merges accepted cases into `summary.csv` after writing a timestamped backup. | — |
| `config.py` | `_LOAD_SCALE` h³ → h⁴; `INP_3D` → finer mesh | `_backup_h4_20261006/` |

Do **not** resume `cases_half_ar_3d_finer_123` with the current `config.py`: it would add h⁴-scaled loads to an h³ folder. `rerun_dropped.py` is safe to use there, because it reads each case's loads from that case's own `.inp`.

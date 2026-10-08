# Beam vs 3D Solid vs Section Solver: Detailed Report

Oct 8, 2026 · @user

## 1. Purpose and setup

The goal is to check how well a one-dimensional beam model, and a cross-section warping solver that feeds it, reproduce the stresses of a full three-dimensional solid. The test problem is a clamped square cantilever of length L = 100, with stresses compared at mid-span (Z = 50).

Three models are compared:

- **1D beam (Abaqus B31).** A Cosserat-rod beam. Its section forces and strains are used to recover stresses.
- **3D solid (Abaqus C3D8R).** Eight-node bricks with one integration point each, on a 10 × 10 section mesh. This is the reference.
- **FEniCS section solver.** A cross-section warping solver (Arora–Kumar–Steinmann 2019) that recovers full 3D stress from beam section resultants. It is solved at E = 1; forces are divided by E before the solve and stresses multiplied by E afterwards.

**Loads and scaling.** Loads are non-dimensional: F0 = EI/L² and M0 = EI/L. Young's modulus is set to E = 1.2e5/h⁴, where h is the section side, so EI = 1e4 for both section sizes (1 × 1 and 0.5 × 0.5), giving F0 = 1 and M0 = 100. Shear loads range over ±0.6 F0 and the tip moment over ±0.7 M0.

**Cases.** 1000 accepted cases plus 300 spares are drawn from `loads_table.csv` (seed 123).

## 2. Method: pipeline, gates and error measures

Each case runs through Abaqus (1D beam and 3D solid) and the FEniCS section solver, and the three stress fields are compared at mid-span.

**Gates.** Cases are filtered automatically so that bad runs cannot hide in the statistics:

- Before running: the table strain gate (≤ 0.005), shear-to-bending ratio (≤ 0.1), and mid-span moment (≥ 0.1 · M0).
- After the Abaqus runs: 1D residual (≤ 1%), hourglass energy ALLAE/ALLIE (≤ 1%), the `strain_admit` check, and FEniCS convergence.

**Error measures.**

- `relL2_resultant`: L2 error divided by the norm of the 3D Mises stress.
- `relL2_Szz`: L2 error of the axial stress, divided by the norm of the 3D S33.
- `L2_over_scale`: error divided by a fixed scale.
- `relL2` (own-norm): diagnostic only.
- `transverse_ratio_fen` / `transverse_ratio_3d`: transverse (shear) stress relative to the axial stress, for FEniCS and the 3D solid.

**Key term: the floor.** The "floor" is the error level that stays about the same in every case, however the loads change. Here it is about 6%: the errors cluster there and never drop much below it.

## 3. The bug: strains rotated twice

The old code gave 11–17% errors because the Abaqus beam strains were rotated a second time. Abaqus B31 outputs (SE, SK, SF, SM) are already in the beam's own axes (the director basis, as in Cosserat-rod strains v0 = Rᵀv and k0 = Rᵀk). The old code applied an extra R\_localᵀ pull-back on SE and SK, a double rotation, which produced a false S11 ≈ S22 stress.

**Fix.** In `batch_driver.py` the strains are now used as they come (`sv_mat = sv_sp`, `sk_mat = sk_sp`). The rotation R is applied only on output, when stresses are reported in the global frame. The same change went into `strain_checks.py`, `recover_section.py` and `trial_noforce.py`.

**Check.** In the pilot re-run the transverse ratio of FEniCS fell to 0.0019–0.0057 and the Mises error to 6.3–7.4%. Several earlier conclusions about twist and size effects were mostly this bug.

**Proof the bug was real.** Before the fix, S11 and S22 along the FEniCS centre line were perfectly correlated (r = 1.00) in `cases_full`, `cases_half_ar_123`, `cases_half_ar_h4_finer` and `cases_quadfix`, and grew with the shear force |V|. The transverse-to-axial stress ratio of FEniCS reached 0.086, against 0.001 for the 3D solid, where a beam theory predicts almost zero. After the fix, the same ratio in the full run is 0.0007–0.0073 for FEniCS (3D: 0.0002–0.0022). Changing the shear-stiffness correction (section 2.16 of the record) and the angle-constraint variants did not change the artefact; only removing the second rotation did.

## 4. Results of the 1000-case run (1 × 1 section)

After the fix, the beam and section-solver stresses match the 3D solid to 5.3–8.4% in Mises stress, with a median of 6.4% (folder `cases_EI_1x1`, 1000 accepted cases).

| Quantity | Before the fix | After the fix |
| --- | --- | --- |
| Mises error vs 3D solid | 11–17%, long tail to 18–21% | 5.3–8.4%, no tail |
| Shear stresses vs axial stress | false in-plane stress present | at most 4% |
| Shear-stiffness anomaly (shear force / shear strain, direction 2) | 0.25 of expected | gone |
| S13 / S23 error gap | 4–8% higher even when converged | 0.2–2.8% in the pilot; settled for h = 1 |

**Spread of the error (1000 cases, `relL2_resultant_Mises`).** Minimum 5.32%, quartiles 6.04% and 6.69%, median 6.36%, maximum 8.39%, standard deviation 0.49 percentage points. The error varies by only about 1 point of its 6% level across wildly different loads, which is why it is called a floor.

The highest errors (7.5–8.4%) occur in cases 769, 406, 141, 887 and 450. They are not outliers; they sit at the top of the same narrow band.

Case 469 finished in Abaqus but was not post-processed; it can be recovered with `rerun_dropped.py`.

## 5. The \~6% floor: probably the 3D reference

The error is steady at about 6% in every case, and the best explanation is that the 3D solid itself reads low at its outer nodes.

**Why.** Bending stress is zero at the centre of the section and largest at the surface. The C3D8R brick calculates stress at one point, in the middle of the brick. In the 10 × 10 mesh the outer bricks are 0.1 wide, so that point sits at 0.45 while the surface is at 0.5, and reads 0.45/0.5 = 90% of the surface stress. The arithmetic predicts an error of about 6.7%; with 20 bricks across it would be about 2.5%.

**Evidence.**

- Along the centre line, the error falls from 6.6% to 1.0% when the outer two nodes are left out.
- FEniCS agrees with beam theory to about 1%, and the 3D solid is about 5% short of both.
- The section-mesh contribution has not been separated out, so this is consistent with the brick explanation but does not exclude other causes.

**Numbers behind the floor explanation** (1000 cases, mid-span, centre line x1 = 0 unless stated):

| Check | Result |
| --- | --- |
| S33 error, all 11 centre-line points | 6.62% median |
| S33 error, 9 inner points only | 1.00% median |
| Full-section Mises error (summary.csv) | 6.40% median |
| Surface stress, 3D relative to FEniCS | 0.905 |
| Centre-line S33 slope, 3D relative to FEniCS | 0.962 median (758 strongly y-bending cases) |
| Departure of the profile from a straight line, as share of its peak | 3D 6%; FEniCS 0.4% |
| Rank correlation of slope shortfall with S33 error | −0.88 |
| Slope relative to beam theory M/I (41 clean bending cases) | FEniCS 0.988; 3D 0.949 |
| Predicted error from the 90% edge reading (11 node rows, two edge rows 10% low) | 6.7% rms, slope 4.5% low |
| Observed | 6.4% median, slope about 4% low |
| Predicted with 20 bricks across | about 2.5% |

The prediction from the brick geometry lands within 0.3 points of the observed value. A 23 September report had also found about a 7% bending-gradient shortfall for C3D8R with 10 elements across.

**What this does not prove.** It is strong evidence, not proof. The "leaving out" figure is a recomputation on saved data (all 11 centre-line points against the 9 inner ones); it removes nothing from the runs, and it is partly circular, since the points were chosen because the theory says they read low. It explains the level of the error, but not why the error rises with tip moment and hourglass energy (section 6). The C3D8I test in section 7 is the direct check, with every node kept.

**C3D8I test result (8 Oct 2026).** The same 30 loads (`loads_subset_elem.csv`: the worst cases 769, 406, 438, 804 and 141, 887, 273, 450 and others, plus 14 of the best) were rerun with `--element-3d C3D8I` on the same 10 × 10 mesh (folder `cases_EI_1x1_C3D8I`). One Abaqus job failed or timed out and a spare took its place, so 29 cases pair with the original run. The loads are identical (mismatch 0.0), the 1D residual is at most 0.6%, and no case failed the hourglass gate.

| Group (Mises error) | C3D8R before | C3D8I after |
| --- | --- | --- |
| Best cases (14) | 5.55% | 1.83% |
| Worst cases (11) | 7.77% | 2.03% |
| Hand-picked worst (4) | 8.13% | 2.64% |
| All 29 (median) | 5.64% | 1.93% |

Every one of the 29 cases improved: the paired drop was 2.7 to 6.2 points (median 4.2), and the new error is about 30% of the old one. The new errors range from 1.4% to 3.1%, against 5.3% to 8.4% before. S33 fell from 5.7% to 1.9% (median), and S13 and S23 fell from about 0.9% to 0.6%. The strain-energy error barely moved (3.9% to 3.2%).

Within these 29 cases the correlation of the error with tip moment fell from 0.86 to 0.16, and with twist from 0.59 to −0.08. Two limits apply: the cases were chosen as best and worst, so these correlations are not comparable with the 1000-case values (0.69 and 0.65), and with n = 29 a correlation below about 0.37 cannot be told from zero. C3D8I has no hourglass modes, so its hourglass ratio is exactly zero; the test therefore cannot separate hourglass from other reduced-integration effects of C3D8R. About 2% error remains unexplained. It may be the coarse 10 × 10 mesh (the arithmetic gave about 2.5% for 20 bricks with C3D8R) or differences between FEniCS and the 3D solid that are not about the brick; a finer mesh is needed to tell them apart.

**Conclusion.** Most of the 6% floor comes from the C3D8R element: switching element, with identical loads and mesh, removed about 70% of the error in every case tested. This confirms the prediction (1–2%) in the median, not in every case.

## 6. Open questions

**Error grows with tip moment and twist.** In `cases_EI_1x1` (1000 cases) the Mises error (`relL2_resultant_Mises`) rises from a median of 5.9% in the lowest fifth of tip moments to 6.7% in the highest fifth. Extrapolating the straight-line fit to zero moment gives 5.9%, so the floor exists even with no moment, and the moment adds about 0.017 percentage points per unit of |Mx|. The table below shows how well each candidate cause tracks the error across the 1000 cases (r = Pearson correlation, R² = share of the variation explained, ρ = Spearman rank correlation).

| Candidate | r | R² | ρ | Reading |
| --- | --- | --- | --- | --- |
| Hourglass energy ratio ALLAE/ALLIE | 0.70 | 0.49 | 0.60 | Tracks the error; range 0.01–0.5%, inside the 1% gate |
| Tip moment \|Mx\| | 0.69 | 0.48 | 0.81 | Tracks the error; strongest rank correlation |
| Twist \|UR3\| | 0.65 | 0.42 | 0.67 | Tracks the error, but see below |
| Total rotation | −0.08 | 0.007 | 0.00 | No simple link on its own |
| Shear force \|V\| | 0.17 | 0.03 | 0.05 | Weak |
| Abaqus shear stiffness (SE2/SE3 ratio, shear factor) | – | – | – | Identical in all 1000 cases (spread 2e-7 and 0), so it cannot explain a difference between cases |

The candidates are not independent, so I separated them: twist and the hourglass ratio are correlated at r = 0.91, and tip moment and hourglass at r = 0.41. After removing the hourglass ratio, the partial correlation of twist with the error falls to 0.05, so twist matters only through the hourglass energy. Tip moment and hourglass each keep a partial correlation of about 0.62 when the other is removed. Together they explain R² = 0.69 of the variation; adding twist gives 0.76.

**Correction to the earlier statement about rotation.** Total rotation has R² = 0.007 on its own, but it is correlated with tip moment (r = 0.41). With the tip moment removed, its partial correlation with the error is −0.55, and adding it to the fit raises R² from 0.76 to 0.83. So "total rotation has no link" was too strong: it has no simple link, and it may act together with the moment. The sign (more rotation, lower error at the same moment) is not understood.

The best-supported reading is that the part of the rise carried by twist comes from hourglass energy growing in the C3D8R bricks, and the rest comes from the tip moment itself. This is a statistical association, not a demonstration of cause. The C3D8I run (section 5) removed most of the dependence on tip moment and twist in 29 hand-picked cases, but because C3D8I has no hourglass modes it cannot separate hourglass from other C3D8R effects.

**Shear-stiffness anomaly.** In the old runs, shear force divided by shear strain in direction 2 was 0.25 of the expected value whatever stiffness was set. The cause was the double rotation; it no longer appears. In the new run the shear factor is 0.8497 (Cowper) in every case, and the Abaqus-to-used shear-strain ratios SE2 and SE3 both equal 1.2495 in every case (spread 2.5e-7). Abaqus's default transverse shear stiffness is 0.68·G·A, against 0.8497·G·A for Cowper. Because this value is the same in all 1000 cases, it cannot cause the rise of the error with moment or twist; it could only shift the overall level, and it does not affect axial stress at mid-span, which comes from the moment.

**Not yet covered.** The 0.5 × 0.5 section (`cases_EI_half`) has not been run, and the older result folders still use the buggy FEniCS input and need `--redo-fenics`.

## 7. Next steps and how to test

1. **C3D8I test (the main check).** `batch_driver.py` now takes `--element-3d C3D8I`, or the environment variable `BATCH_ELEMENT_3D`, with the same 10 × 10 mesh. `make_subset_loads.py` makes a 30-case subset (the worst cases 769, 406, 24, 438, 804 plus some of the best), and `compare_element_runs.py` compares the two runs. Predictions: the floor drops from about 6% to 1–2%, the hourglass energy ratio falls, and the rise with tip moment and twist weakens. Compare the same correlations as in section 6 (error against |Mx|, |UR3| and ALLAE/ALLIE) on the 30 cases. Result (section 5): done on 8 Oct; the median error fell from 5.6% to 1.9%.
2. Run the 0.5 × 0.5 section (`cases_EI_half`).
3. Re-solve older folders with `--redo-fenics`.
4. Run the gate end-to-end test (`gate_test_1x1`) and recover case 469.

**Caveats.** Only the 1 × 1 section has been run in full. The floor explanation is supported by the C3D8I test (median 5.6% to 1.9%), but only on 29 hand-picked cases; a random sample and a finer mesh are still needed, and about 2% is unexplained. Only C3D8\* eight-node bricks can be swapped in; quadratic or tetrahedral elements need a new mesh template. Each code change left a dated backup folder, and the record is in `PARAMETERS_AND_CHANGES.md`.

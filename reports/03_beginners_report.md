# Beam vs Solid vs Section-Solver Study: A Beginner's Report

Oct 8, 2026 · @user

## Summary

After fixing one bug, the FEniCS section solver matches the 3D Abaqus solid to about 6–7% in stress on the five pilot cases, down from 6–11%. The bug was a rotation applied twice to the beam strains before they reached FEniCS.

- **What the study does.** It compares three models of the same clamped beam under random end loads: a 1D Abaqus beam, a full 3D Abaqus solid, and a FEniCS solver that rebuilds the stress inside one cross-section from the 1D beam's results.
- **Why it matters.** If the cheap 1D route reproduces the expensive 3D stresses, a beam model can be trusted in that regime. The study measures where and by how much it fails.
- **The finding.** Abaqus already reports beam strains in the beam's own rotated axes. Our code rotated them again, which produced a false stress in the cross-section and inflated the errors, most of all in high-shear cases.
- **What this changes.** Two earlier explanations, "the error is caused by twist" and "smaller sections fail because of large rotations", were mostly this bug. Every earlier FEniCS result needs re-solving.
- **Where things stand (8 Oct 2026).** The full 1000-case run for the 1 × 1 section has finished, with Mises errors of 5.3–8.4%. The 0.5 × 0.5 run is next. The code is fixed and the changes are logged.

## The problem and the three models

The test object is a cantilever: a long square bar, 100 units long, clamped at one end and pushed on at the free end. Think of a diving board with a square cross-section, either 1 × 1 or 0.5 × 0.5.

The free end gets three loads: a sideways force in X, a sideways force in Y, and a bending moment about X. No twisting load is applied. All results are read at the middle of the bar (station Z = 50), where the three models are compared.

The same bar is modelled three ways:

| Model | Tool | What it is | What it gives us |
| --- | --- | --- | --- |
| 1D beam | Abaqus, B31 elements | The bar as a single line with a few numbers per point | Strains (SE, SK), forces (SF, SM) and rotations (UR) along the bar |
| 3D solid | Abaqus, C3D8R elements | The bar as a block of thousands of small bricks | Stress at every point of the mid-span section. This is the "truth" we compare against |
| Section solver | FEniCS, `warping_core.py` | A small 2D calculation on the cross-section only, driven by the 1D beam's strains and forces | A full stress picture of the section, including how it warps out of plane |

The section solver follows Arora, Kumar & Steinmann (2019). Their idea is a two-step reduction: 3D elasticity becomes a 1D rod (a Cosserat rod), and a rod whose strain varies slowly along its length reduces to a 2D problem on one cross-section (the Helical Cauchy–Born rule).

The central question: given only what the cheap 1D beam knows, how closely does the FEniCS section solver reproduce the 3D stresses, and under which loads does it fail?

Three ideas help throughout:

- **Warping.** A cross-section does not stay flat under torsion or shear. It bulges out of its own plane. Plain beam theory ignores this; the section solver computes it.
- **Material frame.** A rod's strains are measured in axes that rotate with the cross-section. This is what the subscript 0 in v0, k0, n0, m0 means in the paper.
- **Error measure.** Errors are reported as a relative L2 error: the root-mean-square difference between FEniCS and the 3D solid over the 121 section nodes, divided by a typical size of the 3D stress. 6% means the two stress fields differ by about 6% of their size.

## How one run works

Each load case goes through six steps, and two sets of gates decide whether it counts. The driver keeps drawing cases until 1000 are accepted.

&#91;embedded content: one run · 6 steps, 2 discard points\]

A case that fails either gate stage is dropped and replaced by the next of 300 spare loads. The highlighted step is where the bug sat (next section).

### The gates

| Stage | Gate | Rule | Why it exists |
| --- | --- | --- | --- |
| Before Abaqus | Strain | Estimated largest strain ≤ 0.005 (the gate table's formula) | Keeps strains small, the regime the theory assumes |
| Before Abaqus | Shear vs bending | Shear stays within 10% of bending (the table's ratio) | Beam theory suits bending-dominated cases |
| Before Abaqus | Mid-span moment | Bending moment at Z = 50 is at least 0.1 · M0 | Errors are relative, so they blow up where the moment is near zero |
| After Abaqus | 1D convergence | Abaqus residual ≤ 1% of the applied load | An unconverged beam run gives wrong strains |
| After Abaqus | Hourglass | Hourglass energy / internal energy ≤ 1% | The solid's bricks can bend in zero-energy modes that spoil results |
| After Abaqus | Admission | Bending dominates and largest 1D strain < 0.005 | Same intent as the first two gates, checked on the actual beam results |
| FEniCS | Convergence | The section solve must converge | A failed solve has nothing to compare |

### Loads in non-dimensional form

Loads are written as fractions of the beam's own natural scale: F0 = EI/L² for forces and M0 = EI/L for moments, where E is stiffness, I the section's moment of inertia and L the length. The two shear forces are drawn within ±0.6 F0 and the bending moment within ±0.7 M0. Equal fractions make the 1 × 1 and 0.5 × 0.5 sections bend by the same amount, so the two sizes can be compared fairly.

### Stiffness E depends on section size

E is set to 1.2×10⁵ / h⁴, where h is the section side. That keeps EI = 10⁴ for both sizes, so F0 = 1 and M0 = 100 everywhere, and one load file serves both runs. With the old E = 1 the loads were around 10⁻⁶, and Abaqus fell back to its default force tolerance (0.01). That tolerance is larger than the loads, so many 1D runs never really converged.

### FEniCS always solves at E = 1

The FEniCS section solve does not converge reliably with a large E. The material is linear elastic, so stress is simply proportional to E. The code therefore divides the forces by E before the solve and multiplies the stresses by E afterwards.

## What we tried and changed, in order

The study moved from "build the pipeline" to "fix the pipeline" to "find why the errors look the way they do". The table lists what changed and why, newest first.

| When | What changed | Why |
| --- | --- | --- |
| 8 Oct | Stopped rotating the Abaqus beam strains again before giving them to FEniCS. Added unit tests for the new gates (`test_gates.py`). | Root cause of the false in-plane stress and the inflated errors (next section) |
| 7 Oct | FEniCS now solves at E = 1; stresses are multiplied by E afterwards | The section solve did not converge with a large E |
| 7 Oct | Stress-component overlays restored for all six components | Needed to see the in-plane stress problem |
| 7 Oct | 1000 loads plus 300 spares generated exactly as the gate table says; E chosen per section so EI is fixed | The user's table fixed the loads, gates and stiffness; one load file serves both sizes |
| 7 Oct | New gates: 1D convergence, hourglass energy, mid-span moment floor | Old runs had unconverged 1D cases and near-zero-moment cases that dominated the error tails |
| 7 Oct | E raised from 1 to 1.2×10⁵ (a first try at 1×10⁶ was dropped) | Loads near 10⁻⁶ made Abaqus use its fallback tolerance, so 1D runs did not converge |
| 6 Oct | 3D time limit raised from 900 s to 2400 s; loads rescaled by h⁴ instead of h³ | 60 slow 3D cases were dropped though Abaqus finished them; h³ made the small section rotate twice as much |
| Late Sep | Half-size (0.5 × 0.5) section set up; loads can be replayed from a CSV; stale shear stiffness 0.4167 corrected to 0.1042; duplicate `*Static` card removed; quad meshes handled | Needed a second section size and comparable load sets |
| 23 Sep | Separate analysis: the induced torsion is real, and the 5–6% error "floor" is probably the 3D mesh itself | Written up in `twist_and_error_report.pdf` (see the caution below) |
| 7 Sep | Whole pipeline running in WSL with one driver, strain screen and per-case resume | Replaced a fragile Windows/WSL split |
| Aug–Jul | Loads non-dimensionalised; FEniCS solver verified against the paper and against Abaqus on I-beams, rectangles and pentagons | Baseline before the batch study |

**A caution about older reports.** The 23 Sep report and the first README concluded that twist drives the errors and that the small section fails because of larger rotations. Those readings were drawn from results that contained the double-rotation bug. The torsion itself is a real feature of the Abaqus runs: the tip moment keeps its direction in space while the beam tilts. What changed is that twist is not the main cause of the errors: it still goes with a small rise (r = 0.65), but only through hourglass energy. The remaining 6–7% level is still consistent with the 3D mesh being the limit, which the 23 Sep report estimated at about 7% for its bending-stress gradient. That part has not been tested yet.

## The big finding: a rotation applied twice

Abaqus already reports the beam's strains in the beam's own rotated axes, and our code rotated them back a second time before handing them to FEniCS. The result was a false stress in the cross-section.

**An analogy.** Suppose a map is already turned to face the way you are travelling. Turning it again "to face north" leaves it pointing the wrong way. The forces were handed over correctly; the strains were turned a second time. FEniCS received two inputs that no longer described the same state.

**The technical version.** In the Cosserat rod theory the subscript 0 in v0, k0, n0 and m0 means "already pulled back by the rotation R": these are components along the director axes d\_i = R·e\_i. Abaqus B31 output (SE, SK for strains; SF, SM for forces and moments) is in exactly that basis. The driver nevertheless applied Rᵀ to SE and SK, a double rotation. FEniCS enforces both the strains and the forces, through Lagrange multipliers. When the two disagree, the multipliers absorb the conflict as stress that does not exist.

**How it showed up.**

- FEniCS gave equal in-plane stresses S11 and S22 along the section centre line, growing with the shear force. The 3D solid had almost none.
- It was present in every earlier folder, and it hit I1 harder than Mises.
- The force comparison could not catch it, because FEniCS is forced to match the beam's forces by construction.

**How it was tracked down.** Two suspects were ruled out first: the shear stiffness (a corrected shear strain changed nothing) and the constraint that fixes in-plane rotation (variants gave identical results). Re-reading the paper's definitions of the 0-quantities pointed to the strain frame. Giving FEniCS the unrotated strains removed the artifact.

**Before and after** (five pilot cases, 1 × 1 section, re-solved after the fix):

| Measure | Before the fix | After the fix |
| --- | --- | --- |
| In-plane stress in FEniCS, relative to axial stress | 0.009–0.086 | 0.002–0.006 (3D solid: 0.0007–0.0017) |
| Mises error | 6.4–11.3% | 6.3–7.4% |
| I1 error | up to 17% | 5.8–6.6%, close to the Mises error |
| S13 and S23 errors, relative to axial stress | 1–10% | S13 0.2–2.4%, S23 1.6–2.8% |
| In-plane error (S11, S22), relative to axial stress | not recorded | 0.25–0.74% |

In a trial without the force constraints, the moments still matched to within 1e-4 of M0, and the constraint multipliers were about 1e-6. They are now nearly inactive, so they stay on at no cost.

**Small leftovers.** The section model's shear stiffness (strains only) is 0.64–0.84 of G·A depending on the case, and the axial force is off by 0.09–0.45 from higher-order strain terms. The force constraints absorb both with negligible effect on stress.

**What this does not change.** The load cases and the pass/fail decisions are the same, because the gates use the raw Abaqus values. Only the FEniCS stresses change, so older folders are re-solved rather than re-run in Abaqus.

## Where the runs stand

The 1 × 1 run is finished and checked; the 0.5 × 0.5 run is still to start. The numbers below are from 8 Oct 2026.

### Current set-up

| Setting | Value |
| --- | --- |
| Beam | Length 100, mid-span station Z = 50, square section of side h = 1 or 0.5 |
| Stiffness | E = 1.2×10⁵ / h⁴ (1.2×10⁵ at h = 1, 1.92×10⁶ at h = 0.5), Poisson's ratio 0.3 |
| Loads | Shears within ±0.6 F0, moment within ±0.7 M0; F0 = 1 and M0 = 100 for both sizes |
| 1D model | Abaqus B31, 100 elements, nonlinear geometry, default transverse shear stiffness |
| 3D model | Abaqus C3D8R, 10 × 10 bricks across the section, nonlinear geometry |
| Section mesh for FEniCS | 10 × 10 quads, solved at E = 1 |
| Error measure | Relative L2 error against the 3D stress, taken over the whole section |

### Run inventory

| Folder | Section | Status |
| --- | --- | --- |
| `pilot_table_1x1` | 1 × 1 | 5 cases. FEniCS re-solved after the fix; results in the previous section |
| `pilot_table_half` | 0.5 × 0.5 | Obsolete: made with the old E. Needs re-running |
| `gate_test_1x1` | 1 × 1 | Test of the gates with deliberately tight limits. Prepared, not run yet |
| `cases_EI_1x1` | 1 × 1 | Full run, finished 8 Oct: 1000 accepted cases (results below) |
| `cases_EI_half` | 0.5 × 0.5 | Full run, to start after the 1 × 1 run |
| Older folders (`cases_full`, `cases_quadfix`, `cases_half_ar_*`) | both | Used the double rotation, and some used E = 1. Keep as history; re-solve with `--redo-fenics` to compare |

### Typical errors before the fix

These medians come from the 23 Sep analysis of 1000 cases (1 × 1 section) and include the bug. The error is each component's difference from the 3D solid, relative to that component's own size.

| Quantity | Median error | Reading |
| --- | --- | --- |
| S33 (axial stress) | 6.4% | The main bending stress. A floor of about 5–6%, not a long tail |
| Mises | 6.6% | Overall stress intensity |
| I1 | 7.5% | Sum of the normal stresses; the one the bug inflated most |
| J3 | 21% | About three times the S33 error, because J3 is cubic in stress |
| S13, S23 | 47%, 15% | Shear stresses are tiny, so relative errors look large; judge them against the axial stress instead |

### What we know about the gates

- The 1D-convergence, hourglass and mid-moment gates pass synthetic pass/fail tests (`test_gates.py`).
- Applied to an older E = 1 folder, the 1D gate would have discarded 411 of 1000 cases; 546 had hit Abaqus's fallback tolerance. Such cases fed wrong strains into FEniCS in the old runs.
- No real case has yet failed the hourglass gate. The pilot cases sit at 0.02–0.26%, well under the 1% limit.

## Results of the full 1 × 1 run

All 1000 cases of the 1 × 1 run finished on 8 Oct, and the error has a floor of about 6%, not a long tail. Mises error runs from 5.3% to 8.4%, where earlier runs reached 18–21%.

One case (469) was lost: both Abaqus jobs finished but nothing was post-processed, so the first spare took its place. No case failed the convergence or hourglass gates, and the in-plane stress ratio in FEniCS stayed at or below 0.0073.

| Measure (percent) | Lowest | Median | 99th percentile | Highest |
| --- | --- | --- | --- | --- |
| Mises | 5.3 | 6.4 | 7.8 | 8.4 |
| I1 | 5.4 | 6.2 | 7.4 | 8.2 |
| S33 (axial stress) | 5.4 | 6.4 | 7.5 | 8.3 |
| J3 | 2.6 | 3.5 | 4.2 | 4.8 |
| S13 (shear) | 0.01 | 0.9 | 2.8 | 3.2 |
| S23 (shear) | 0.1 | 1.6 | 3.5 | 4.0 |
| Strain energy | 0.01 | 1.4 | 4.4 | 5.9 |

Stresses are measured against the 3D axial stress; energy against its own size.

### Why the cases with larger errors are larger

- **The floor is probably the 3D model.** Each 3D brick holds one stress value, taken at its centre. A node on the surface only has a brick centre 0.45 from the middle, not 0.5, so it reads about 90% of the true bending stress. The data agree: the outer 3D nodes read 90.5% of the FEniCS value. Arithmetic with that single effect predicts an error of about 6.7%; we measure a median of 6.4%. Along the centre line the error is 6.6% with all 11 points and 1.0% without the two outermost. Against exact beam theory, FEniCS is within about 1% and the 3D solid about 5% short (41 clean bending cases). This strongly suggests the 3D mesh, but it is not proven. The test was a better brick element (C3D8I) on the same mesh, run on 8 Oct for 30 loads (29 compared): the median Mises error fell from 5.6% to 1.9%, and every case improved by 2.7 to 6.2 points. So most of the floor comes from the C3D8R brick. About 2% is still unexplained (perhaps the coarse mesh; 20 bricks across would test that), and the cases were hand-picked, so a random sample is still needed. In these cases the correlation of the error with tip moment fell from 0.86 to 0.16.
- **Cases above the floor** (Mises 7.5–8.4%: cases 769, 406, 141, 887, 450) have a large tip moment and a good deal of twist. The error grows with the tip moment (5.9% to 6.7%) and with twist (6.1% to 7.4% for the centre-line S33), but the rise with twist sits entirely in the two outer nodes of the 3D solid; without them it disappears. The beam's total tilt has no simple link to the error on its own (R² = 0.007 over 1000 cases), so large rotations are unlikely to be the cause; once the tip moment is taken into account it does carry some information (partial r = −0.55), so it is not fully ruled out. Over the 1000 cases the error tracks the hourglass energy ratio (r = 0.70), the tip moment (r = 0.69) and twist (r = 0.65); twist works only through hourglass energy (partial r = 0.05 once that is removed). Abaqus's shear stiffness is identical in every case, so it cannot explain the rise. Whether the mid-span moment matters is not established: its correlation with the error is weak (r = 0.21). Tip moment and twist together explain 57% of the Mises variation (R² = 0.57); adding the hourglass ratio takes it to 76%.
- **Largest shear errors** (S13 and S23 at 3–4%: cases 24, 27, 9, 438, 305) have the biggest shear forces (near the 0.6 limit) and the lowest bending dominance (about 4). Shear stress is simply a bigger share of the answer there.
- **Largest energy errors** (4–6%: cases 804, 382, 754, 687) have a small moment at mid-span. The stored energy is then small, so a fixed small miss looks large in relative terms.

## What is still open and what to do next

The next step is to finish the two full runs and confirm that the fix holds across 1000 cases at both sizes.

### Next steps

1. `S`tart `cases_EI_half` (the commands are in `PARAMETERS_AND_CHANGES.md`, section 4).
2. After each run, check three things: the in-plane stress ratio is below about 0.01 for nearly all cases; the Mises and I1 errors are close; and the 0.5 and 1 sections give similar errors for the same non-dimensional load.
3. Optionally run the gate test, with its deliberately tight limits, to see the discard and spare-load logic work end to end.
4. Re-solve the older folders with `--redo-fenics` to see how much of the old error tail and size effect remains.
5. Follow up the 8 Oct test, which re-ran 30 cases with a better brick element on the same mesh (median error 5.6% to 1.9%), by running a random sample of unselected cases (C3D8I, set with --element-3d; commands in the notes, section 9), or with 20 bricks across.

### Open questions

- **Shear stiffness anomaly.** Resolved in practice. With the default shear stiffness, the shear force divided by shear strain is the same in both directions (0.68 · G · A in all 995 cases with shear). The old 0.25 pattern belonged to the old templates and its cause was never found.
- **S13 and S23 gap.** Settled for the 1 × 1 section: errors are 0.01–4% of the axial stress. Whether a gap appears between sizes needs the 0.5 × 0.5 run.
- **3D first increment.** The two 3D templates start with different step sizes (1.0 and 0.01). Results should not depend on it, but the templates should match.
- **Older notes.** The first README and the 23 Sep report still say that twist drives the errors. They predate the fix. The original `warping_core.py` has a comment telling users to rotate the strains back; that advice is wrong for Abaqus beam output, and the file is deliberately untouched.

### Glossary

| Term | Plain meaning |
| --- | --- |
| Cantilever | A bar clamped at one end and loaded at the other |
| Mid-span (Z = 50) | The cross-section halfway along the bar, where all comparisons are made |
| Strain | How much the material stretches or shears, as a fraction |
| Stress | Force per unit area inside the material. S33 is the axial stress; S13 and S23 are shear stresses |
| Mises | One number that summarises how hard the material is being loaded overall |
| I1 | The sum of the three normal stresses (the average push or pull) |
| J3 | A cubic stress measure that tells tension-dominated from compression-dominated states |
| B31 | Abaqus's 1D beam element |
| C3D8R | Abaqus's 3D brick element with one integration point |
| Hourglass energy | Energy wasted in zero-stiffness wiggles of the brick elements; a quality warning |
| Warping | A cross-section bulging out of its own plane under twist or shear |
| Cosserat rod | A beam model where each point of the centre line also carries a rotating frame (the directors) |
| Material (director) frame | Axes that rotate with the cross-section; the 0 in v0, k0, n0, m0 |
| Lagrange multiplier | A helper variable that forces a condition (such as "match this force") to hold exactly |
| Non-dimensional load | A load written as a fraction of the beam's natural load scale, F0 or M0 |
| Gate | A pass/fail check that decides if a case is kept |
| Spare load | An extra load case held in reserve to replace a discarded one |
| Nonlinear geometry (nlgeom) | Abaqus tracks large rotations instead of assuming small ones |
| relL2 | Relative root-mean-square error between two fields, as a percentage |

### Sources

All read from the project's log folder and the batch\_test folder on the user's computer; no outside sources.

- `PARAMETERS_AND_CHANGES.md` (batch\_test and batch folders): parameters, change log and observations to 8 Oct 2026
- `README.md` in the batch\_test folder: tail analysis and size-effect analysis, 6 Oct
- `README.md` in the log folder: pipeline, checks, frames and section-size notes, 28 Sep
- `README_session_changes.md`: half-section debugging session notes
- `project_reference_summary.md` and `pipeline_explained.tex`: pipeline reference and beginner's guide, 7 Sep
- `twist_and_error_report.pdf`: twist mechanism and error decomposition, 23 Sep
- `beam_nondim.pdf`: load non-dimensionalisation derivation, 20 Aug
- `gowri_report.pdf`, `gowri_update.pptx` and the `3D_stress_recovery*.pptx` decks: the dimensional-reduction method (Arora, Kumar & Steinmann, 2019) and its verification
- `LOG.txt`: notes on the solver's functions
- `checks [Autosaved].pptx`: planned verification tests
- `lab_meeting_tex.pdf` and the zip archives were not read in detail

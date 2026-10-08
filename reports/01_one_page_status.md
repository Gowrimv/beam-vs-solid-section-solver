# Beam vs 3D Solid vs Section Solver: One-Page Status

Oct 8, 2026 · @user

After fixing a double-rotation bug, the FEniCS section solver matches the 3D Abaqus solid to a median 6.4% in Mises stress over 1000 load cases on the 1 × 1 section. The remaining \~6% floor appears to come from the 3D mesh, not the solver.

## What we did

- Compared three models of a clamped square cantilever (L = 100): an Abaqus 1D beam (B31), an Abaqus 3D solid (C3D8R, 10 × 10 section mesh), and a FEniCS cross-section warping solver (Arora–Kumar–Steinmann 2019).
- Generated 1000 random load cases (shears up to ±0.6 F0, tip moment up to ±0.7 M0), with non-dimensional loads F0 = EI/L², M0 = EI/L, and E = 1.2e5/h⁴ so EI = 1e4 for both section sizes.
- Compared stresses at mid-span (Z = 50): Mises/resultant error against the 3D solid.
- Added automatic gates: load-table checks before running, then 1D residual ≤ 1%, hourglass energy ≤ 1%, and FEniCS convergence after.

## What we found

- **A bug, now fixed.** The code rotated Abaqus beam strains a second time, though they were already in the beam's own axes. That created a false in-plane stress and 11–17% errors. After the fix the false stress is gone.
- **Full 1000-case run (1×1 section):** Mises error is 5.3–8.4% (median 6.4%), shear stresses are at most 4% of the axial stress, and the long error tail (18–21% before) has disappeared. The earlier shear-stiffness anomaly and the S13/S23 gap are gone too.
- **A steady \~6% floor.** The best explanation is the 3D reference itself: with one integration point per brick, the outer nodes read about 90% of the true bending stress. Recomputed on saved data, the centre-line S33 error is 6.6% with all 11 points and 1.0% with the 9 inner ones (nothing is removed from the runs). The surface stress of the 3D solid is 0.905 of FEniCS's; the 3D slope is 0.949 of beam theory against 0.988 for FEniCS (41 clean bending cases); slope shortfall and error correlate at −0.88 (rank). The 90% reading predicts 6.7%; observed 6.4%. Strong evidence, not proof.
- The error also rises from 5.9% to 6.7% from the lowest to the highest fifth of tip moments. Across the 1000 cases it tracks the hourglass energy ratio (r = 0.70), tip moment (r = 0.69) and twist (r = 0.65); twist acts only through hourglass energy (partial r = 0.05 once hourglass is removed). Shear stiffness is identical in every case, so it cannot explain the rise. Total rotation alone has R² = 0.007 but adds explanatory power next to moment, so I am not calling it ruled out.

## Next steps and caveats

- Run \~30 cases with a better brick (`--element-3d C3D8I`, same mesh). Prediction: floor falls from \~6% to 1–2%. If it doesn't, the floor explanation is wrong.
- Run the 0.5×0.5 section (`cases_EI_half`), and re-solve older folders with `--redo-fenics`.
- Caveats: the 0.5 section has not been run; the gate end-to-end test and case 469 recovery are pending; the floor explanation is unproven until the C3D8I test is done.

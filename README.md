# Beam vs 3D solid vs section solver

Comparison of an Abaqus 1D beam (B31), an Abaqus 3D solid (C3D8R by default) and a
FEniCS cross-section warping solver (Arora, Kumar, Steinmann 2019) for a clamped
square cantilever (L = 100, stresses compared at mid-span Z = 50).

## Layout

| Folder | Contents |
| --- | --- |
| `batch/` | The pipeline: `batch_driver.py` (main), `abaqus_io.py`, `fenics_solve.py`, `compare.py`, gates in `strain_checks.py`, tests in `test_gates.py`, plus `README.md` and `PARAMETERS_AND_CHANGES.md` (full record of parameters, changes, observations) |
| `inputs/` | Abaqus input templates (h = 1 and h = 0.5) and the section geometry |
| `loads/` | Load tables (`loads_table*.csv`, seed 123) and test subsets (`loads_gate_test.csv`, `loads_subset_elem.csv`, 30 cases for the element test) |
| `results/` | `summary.csv` and accepted-load tables of the 1000-case 1x1 run, the pilot, and the record of the C3D8I run settings |

Large Abaqus output (.odb, .dat, case folders), meshes and the dated `_backup_*`
folders are not in the repository (see `.gitignore`).

## Not included

`warping_core.py` from the FEniCS section-solver code (Arora, Kumar, Steinmann 2019)
is required but is not redistributed here. Put your copy in `FEniCSfiles/` as before;
this repository never edits it.

## Headline results (see `batch/PARAMETERS_AND_CHANGES.md` sections 7-9)

- Bug fixed: Abaqus beam strains were rotated twice (they are already in the director basis).
- Full 1x1 run (1000 cases): Mises error vs 3D solid 5.3-8.4%, median 6.4%.
- The ~6% floor is probably the 3D reference (C3D8R reads outer nodes at ~90%); strongly
  suggested, not proven. Test: `--element-3d C3D8I` on `loads/loads_subset_elem.csv`.
- Error rises slightly with tip moment (5.9% to 6.7%) and with the hourglass energy ratio
  (r = 0.70, 0.69); twist acts through hourglass; shear stiffness is identical in all cases.

## Run

```
cd batch
BATCH_SECTION_SIDE=1 python3 batch_driver.py --cases 1000 ...   # see batch/README.md
python3 test_gates.py                                           # gate tests
python3 batch_driver.py ... --element-3d C3D8I                  # other 8-node brick
```

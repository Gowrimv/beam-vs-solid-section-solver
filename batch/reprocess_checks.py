"""
reprocess_checks.py
===================
DEPRECATED SHIM -- kept only so any existing `import reprocess_checks`
keeps working. All real code now lives in strain_checks.py.

WHY THIS FILE EXISTS AS A SHIM
------------------------------
reprocess_checks.py used to be a byte-for-byte COPY of strain_checks.py
(same 22 functions, same docstrings, no __main__). The two then drifted:
a fix to rotated_beam_forces() -- correcting the section-moment ordering
from [SM3, SM2, SM1] to [SM1, SM2, SM3] -- was applied HERE, in the copy,
while batch_driver.py and compare.py both import strain_checks. The fix
therefore never executed, and every force_moment_1d3d_error.csv on disk
was written with the transposed ordering (Mx_1d holding SM3, the torque,
and Mz_1d holding SM1, a bending moment), producing reported errors of
-96.8% and +2974% where the true agreement is ~1%.

Keeping a second copy of an 841-line module is how that happened, so the
copy is gone rather than re-synced. Import from strain_checks instead.
"""
from strain_checks import *          # noqa: F401,F403
import strain_checks as _sc

# re-export the private helpers too, since `import *` skips leading-_ names
_solid_corner_rows = _sc._solid_corner_rows
_solid_centerline_row = _sc._solid_centerline_row
_root_strain_and_mid_bd_from_beam_csv = _sc._root_strain_and_mid_bd_from_beam_csv
_mid_strain_and_bd_from_beam_csv = _sc._mid_strain_and_bd_from_beam_csv

# module-level constants that `import *` would otherwise skip only if
# __all__ were defined; re-bound explicitly for safety.
MAX_POST_STRAIN_DEFAULT = _sc.MAX_POST_STRAIN_DEFAULT

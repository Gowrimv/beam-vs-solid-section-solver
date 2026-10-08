"""
config.py
=========
All fixed, cheap-to-import configuration for the batch driver: file
paths, section/material constants, sampling ranges, the non-dimensional
scales (eq. 9-12, see batch_driver.py's module docstring for the full
derivation/rationale) used to report every "post" (comparison-stage)
quantity as a dimensionless fraction, and the one shared logging helper
(tprint).

Deliberately dependency-light (stdlib `os`/`threading` only) so ANY
other module can import it safely — including plot_results.py, which is
meant to regenerate plots WITHOUT Abaqus, FEniCS, or pyvista installed
at all.
"""
import os
import threading

# ── shared logging ──────────────────────────────────────────────────────────
_print_lock = threading.Lock()


def tprint(*args):
    """Thread-safe print (holds a lock so concurrent worker threads in
    a ThreadPoolExecutor don't interleave partial lines with each
    other)."""
    with _print_lock:
        print(*args, flush=True)

# ── paths ────────────────────────────────────────────────────────────────────
HERE       = os.path.dirname(os.path.abspath(__file__))
PROJECT    = '/mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test'
ABAQUS_BAT = '/mnt/c/SIMULIA/Commands/abaqus.bat'

# INP_1D / INP_3D / MESH_XMF are chosen from SECTION_SIDE below.
DEFAULT_CASES_DIR = os.path.join(PROJECT, 'cases')
EXTRACT_SCRIPT    = os.path.join(HERE, 'extract_odb.py')
FREEEBODY_SCRIPT = os.path.join(HERE, 'extract_freeBody.py')
# ── section / material ─────────────────────────────────────────────────────
# E is scaled to 1e6 (was 1.0) so Abaqus's absolute 'zero force' criterion
# never triggers: with E = 1 the loads were ~1e-7 and Abaqus accepted 1D
# solutions with residual forces up to 300% of the applied shear. The model
# is linear elastic and every output is normalised (F0/M0, rel-L2), so this
# changes no result. The *_E1e6.inp templates carry the same E.
# Backup of the previous config: _backup_Escale_20261007/
# 7 Oct 2026: E = 1.2e5 per the gate table, so F0 = EI/L^2 = 1 and M0 = 100
# for h = 1 (h = 0.5: F0 = 0.0625, M0 = 6.25) -- still far above Abaqus's
# zero-force level. Was 1e6 (backup: _backup_E1p2e5_20261007/).
E_REF = 1.2e5;  NU = 0.3;  L = 100.0
E = E_REF   # replaced below once SECTION_SIDE is known
# FEniCS is always solved with E = FENICS_E (= 1): its nonlinear solve does
# not converge reliably with large E. fenics_solve.py divides the force /
# moment targets by E/FENICS_E before the solve and multiplies every
# stress-unit output (S11..S23, Mises, inv_1 in the VTUs and
# *_section_X0.csv; n, m, Phi) by the same factor afterwards. Strains and
# displacements do not depend on E. (added 7 Oct 2026)
FENICS_E = 1.0
# Shear strains passed to FEniCS (added 7 Oct 2026). Abaqus's B31 default
# transverse shear stiffness is softened (SF/SE = 0.68 G A in the pilot), so
# its shear strains disagree with the FEniCS section's own shear stiffness
# while FEniCS is ALSO given the shear forces; the Lagrange multipliers then
# add spurious in-plane stresses (transverse_ratio_fen up to 0.086). Instead,
# pass SE_shear = SF_shear / (SHEAR_K_FACTOR * G * A). Default: Cowper's
# factor for a solid rectangle, 10(1+nu)/(12+11 nu) = 0.8497 at nu = 0.3.
# Tune after the pilot so transverse_ratio_fen ~ transverse_ratio_3d.
# Set CORRECT_SHEAR_STRAINS = False to pass Abaqus's SE2/SE3 unchanged.
CORRECT_SHEAR_STRAINS = True
SHEAR_K_FACTOR = 10.0 * (1 + NU) / (12.0 + 11.0 * NU)

# Square cross-section side length. CHANGE ONLY THIS -- A, I, c and the
# load ranges below are all derived from it. It must match the section in
# INP_1D (*Beam Section, RECT), the geometry in INP_3D and the FEniCS mesh
# MESH_XMF; batch_driver.py warns at start-up if INP_1D disagrees.
# Pick the section per run without editing this file:
#   BATCH_SECTION_SIDE=1   python3 batch_driver.py ...   (default)
#   BATCH_SECTION_SIDE=0.5 python3 batch_driver.py ...
SECTION_SIDE = float(os.environ.get('BATCH_SECTION_SIDE', '1'))
_TEMPLATES = {
    # side: (1D template, 3D template, FEniCS section mesh)
    1.0: ('1d_fx_fy_mx_E1p2e5.inp',         '3d_fx_fy_mx_E1p2e5.inp',               'slice_50_quad_native.xdmf'),
    0.5: ('1d_fx_fy_mx_half_ar_EI1e4.inp',  '3d_fx_fy_mx_half_ar_finer_EI1e4.inp', 'slice_ar_50_quad.xdmf'),
}
if SECTION_SIDE not in _TEMPLATES:
    raise SystemExit('config.py: no templates for SECTION_SIDE=%g' % SECTION_SIDE)
# 7 Oct 2026 (new table): E chosen per section so that E*I is the same for
# every size: E = E_REF / SECTION_SIDE**4 -> E I = 1e4, F0 = EI/L^2 = 1 and
# M0 = EI/L = 100 for h = 1, 0.5, 0.25 alike. The SAME physical loads
# (Vx = v_x, Vy = v_y, Mx = 100 m) then apply to every run.
#   h = 1: E = 1.2e5   h = 0.5: E = 1.92e6   h = 0.25: E = 3.072e7
E = E_REF / SECTION_SIDE**4
INP_1D   = os.path.join(PROJECT, _TEMPLATES[SECTION_SIDE][0])
INP_3D   = os.path.join(PROJECT, _TEMPLATES[SECTION_SIDE][1])
MESH_XMF = os.path.join(PROJECT, _TEMPLATES[SECTION_SIDE][2])
A  = SECTION_SIDE**2                 # area
I  = SECTION_SIDE**4 / 12.0          # second moment of area (either axis)
c  = SECTION_SIDE / 2.0              # centroid -> extreme fibre
G  = E / (2*(1+NU))
Z_MID = 50.0

# ── sampling ─────────────────────────────────────────────────────────────────
STRAIN_THRESHOLD = 5e-3
# Load ranges are defined for a 1 x 1 section and scaled by SECTION_SIDE**3
# (= (I/c) ratio), which keeps the BENDING STRAIN distribution -- and so the
# strain-screen pass rate (~62%) -- identical for any section size. Checked
# by Monte Carlo: 1x1 and 0.5x0.5 both give 62.4% pass, median strain
# 3.73e-3. Side effects of a smaller section at the same strain: shear
# strain scales by SECTION_SIDE (bending dominance ~ 1/SECTION_SIDE) and
# curvature/rotations by 1/SECTION_SIDE.
#   SECTION_SIDE = 1.0 -> V +-5e-6,    Mx +-5e-4
#   SECTION_SIDE = 0.5 -> V +-6.25e-7, Mx +-6.25e-5
# Loads are set directly in non-dimensional units: +-0.6 F0 and +-0.6 M0
# with F0 = E I / L^2, M0 = E I / L, for ANY section size and E. Equal
# non-dimensional loads -> equal bending rotations and twist across sizes
# (strain then scales with c/L). Same box as the old 1x1 ranges
# (V +-5e-6, Mx +-5e-4 at E = 1).
# 7 Oct 2026: per the gate table, m_x,tip ~ U(-0.7, 0.7) and
# v_x, v_y ~ U(-0.6, 0.6). (Replayed --loads-csv runs are unaffected.)
LOAD_ND   = 0.6     # shears, in F0 units
LOAD_ND_M = 0.7     # tip moment, in M0 units
VX_RANGE = (-LOAD_ND*E*I/L**2,    LOAD_ND*E*I/L**2)
VY_RANGE = (-LOAD_ND*E*I/L**2,    LOAD_ND*E*I/L**2)
MX_RANGE = (-LOAD_ND_M*E*I/L,     LOAD_ND_M*E*I/L)

# Gate-table thresholds (7 Oct 2026)
BENDING_GATE_MAX = 0.1    # pre-run: max transverse-shear / bending stress at Z_MID
HOURGLASS_MAX    = 0.01   # post-run: max ALLAE / ALLIE in the 3D solid
# Mid-span floor (gate table): (1/200)(|m_tip - v_y/2| + |v_x/2|) >= 0.1 eps_tol
# for h = 1, L = 100, i.e. |m_x,mid| + |m_y,mid| >= 0.1 * 0.005 * 200 = 0.1 in
# M0 units. Kept in M0 units (not strain) so the SAME draws are dropped at
# every section size and 1x1 / 0.5x0.5 runs stay paired. "Tune after a pilot."
MID_FLOOR_ND     = 0.1

# 3D solid element type (added 8 Oct 2026). None = keep whatever the 3D
# template says (C3D8R). Set with the BATCH_ELEMENT_3D environment variable
# or batch_driver.py --element-3d TYPE. Only 8-node brick types can be
# swapped in (C3D8R, C3D8, C3D8I, C3D8H, C3D8RH, ...): the template's
# connectivity is first-order, so quadratic bricks (C3D20R...) or
# tetrahedra need a new mesh and a matching FEniCS section mesh.
ELEMENT_3D = os.environ.get('BATCH_ELEMENT_3D') or None
N_TARGET     = 1000
MAX_PARALLEL = 3   # concurrent cases; 3 x 2 jobs x 5 tokens = 30 tokens

# ── POST-run admission gate — see LOAD ACCEPTANCE LOGIC in
# batch_driver.py's module docstring. Evaluated at Z_MID from the REAL
# Abaqus 1D section forces/moments (SF1-3, SM1-3) as a STRESS ratio:
#   sigma_max = |N|/A + (|M1|+|M2|) c/I
#   tau_max   = 1.5 |V|/A + |T|/(0.208 a^3)
# The case passes if sigma_max / tau_max >= this value. 1.0 means "the
# axial stress is the largest stress"; raise it for a stricter margin. ──
MIN_BENDING_DOMINANCE_POST = 1.0   # tune with --min-bending-dominance-post

# ── non-dimensional scales ────────────────────────────────────────────────────
# F0 (eq. 9):      force scale        = E I / L^2
# M0 (eq. 10):     moment scale       = F0 * L = E I / L
# sigma0 (eq. 11): stress scale       = M0 * c / I = E * c / L
#                  The bending stress a unit non-dimensional moment
#                  (M = M0) produces, so S/sigma0 = m for pure bending at
#                  ANY section size. Was E * STRAIN_THRESHOLD until 7 Oct
#                  2026: that left a factor c/L, so at equal non-dimensional
#                  loads the 0.5x0.5 stresses came out half the 1x1 ones.
#                  Backup: _backup_stressscale_20261007/
# U0 (eq. 12):     energy scale       = F0 = E I / L^2
#                  (energy per unit length has force units, so it shares
#                  F0's scale — not a coincidence, see eq. 12 note)
FORCE_SCALE  = E * I / L**2        # F0
MOMENT_SCALE = FORCE_SCALE * L     # M0 = F0 * L
STRESS_SCALE = MOMENT_SCALE * c / I   # sigma0 = M0 c / I = E c / L
ENERGY_SCALE = FORCE_SCALE         # U0 = F0
J3_SCALE     = STRESS_SCALE**3     # J3 = det(deviatoric stress) has
                                    # stress^3 units, so it gets its own
                                    # scale (sigma0^3), not sigma0 itself.

# Stress-unit (degree-1) components/invariants that get divided by sigma0
# in plots. J3 is degree-3 in stress and is normalized separately by
# J3_SCALE = sigma0**3 wherever it appears (see compare.py).
STRESS_UNIT_COMPS = {'S11','S22','S33','S12','S13','S23','Mises','inv_1'}
# Per-case FEniCS matplotlib plots (e.g. *_Stress_centerline.png). They are
# drawn INSIDE the FEniCS lock, so with several workers they hold up every
# other case's FEniCS solve. False = skip them (VTU + CSV outputs unchanged).
FENICS_PLOTS = False   # 7 Oct 2026: off -- FEniCS draws them at E = FENICS_E (1), before the stress rescale, so their stress values would be wrong. compare.py's contours/centerline plots are correct.

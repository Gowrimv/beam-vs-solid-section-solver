"""
strain_checks.py
=================
Pure numpy/pandas post-processing math for the batch pipeline: the
pre-run analytical (Mohr's-circle) screen, the post-run real-Abaqus
strain/bending-dominance checks, the frame-rotation helpers, the
two-part post-Abaqus admission check (strain_admit(): bending dominance
+ 1D small strain), the 1D-vs-3D force/moment cross-check, and the
stress invariants used by compare.py. Depends on nothing but numpy/pandas/config -- NO Abaqus,
FEniCS, or pyvista here -- so this module (like plotting.py) can be
imported anywhere those aren't installed, e.g. for offline re-analysis
of an existing summary.csv/beam_whole.csv/solid_mid.csv without needing
the full simulation stack.

Split out of batch_driver.py so the strain/screening logic -- which is
reused by run_case(), load_existing(), AND compare.py's contour/overlay
plots -- has one home instead of being duplicated or creating an import
cycle between run_case() and the comparison stage.
"""
import os
import math
import numpy as np
import pandas as pd

import config
from config import tprint


# ════════════════════════════════════════════════════════════════════════════
# POST-ABAQUS ADMISSION CHECK — defaults
# ════════════════════════════════════════════════════════════════════════════
# A case is admitted only if BOTH hold, evaluated from REAL Abaqus 1D
# (beam) output -- the 3D solid is NOT used for admission:
#   1. axial stress dominates -> bending_dominance_mid_abaqus = sigma_max/tau_max
#      (1D SF/SM at Z_MID, torsion included) >= config.MIN_BENDING_DOMINANCE_POST
#      (or --min-bending-dominance-post).
#   2. the 1D post-run principal strain, max over every beam_whole.csv
#      station, is under MAX_POST_STRAIN_DEFAULT (or --max-post-strain).
# See strain_admit() below for the actual combination logic.
MAX_POST_STRAIN_DEFAULT   = 5e-3   # "max strains are small (under 5e-3)"

def strain_screen(vx, vy, mx):
    """PRE-run pass criterion. Evaluated in main() on raw (Vx,Vy,Mx) draws,
    BEFORE any Abaqus job exists — there is no FE output to check yet, so
    this has to be an analytical (closed-form beam-theory / Mohr's-circle)
    estimate of the max principal strain at the clamped root. It exists to
    decide which random draws are worth spending Abaqus licence time on,
    and is deliberately a conservative worst-case bound. It is NOT a
    validity check on simulation results — see
    `principal_strain_from_beam_row` / `root_strain_under_threshold` below for the
    post-run counterpart that uses real Abaqus output."""
    # CONSERVATIVE BOUND, not the signed moment. This function is a GATE:
    # it must never UNDER-estimate, or a too-severe draw slips through.
    # The true root moment is mx - vy*L (see moment_at_station), but that
    # can be smaller than either term when they cancel, so using it here
    # under-estimates for some sign combinations. |mx| + |vy|*L is >= the
    # true |moment| for every sign combination, and matches the form
    # range_safety_check() already uses for its exact box certification.
    # Measured against the real Abaqus root strain over 999 cases:
    #   mx + vy*L   (old)   rho 0.042, under-estimates 582/999  <- inverted
    #   mx - vy*L   (true)  rho 0.9999, under-estimates  62/999 <- not a bound
    #   |mx|+|vy*L| (this)  rho 0.822,  under-estimates   0/999 <- safe
    mx_root = abs(mx) + abs(vy) * config.L;  my_root = abs(vx) * config.L
    eps = abs(mx_root)*config.c/(config.E*config.I) + abs(my_root)*config.c/(config.E*config.I)
    gam = (abs(vx) + abs(vy)) / config.A / config.G
    return float(eps/2 + np.sqrt((eps/2)**2 + (gam/2)**2))

def moment_at_station(vx, vy, mx, z):
    """Internal bending moments M_x(z), M_y(z) at station z (measured
    from the clamped root, z=0, same convention as everywhere else in
    this file) for a cantilever loaded only at the free tip (z=L) —
    static equilibrium of the free portion beyond z, no distributed
    load, so shear force is CONSTANT along the whole beam while moment
    grows linearly with the remaining lever arm (L - z):
        M_x(z) = mx - vy * (L - z)      (z=0 -> mx_root formula above)
        M_y(z) = vx * (L - z)           (z=0 -> my_root formula above)

    SIGN (fixed): the tip force Fy contributes r x F with r = (0,0,L-z)
    and F = (0,Fy,0), giving (-(L-z)*Fy, 0, 0) -- it SUBTRACTS from a
    positive applied mx. Verified against Abaqus SM1 at mid-span over
    1000 cases: rho(mx - 50*vy, SM1) = 0.99968 and median ratio 1.020,
    versus rho = 0.645 / ratio 0.513 for the previous '+' sign.
    """
    return mx - vy * (config.L - z), vx * (config.L - z)

def strain_screen_at_station(vx, vy, mx, z):
    """Same closed-form Mohr's-circle principal-strain formula as
    strain_screen(), but evaluated at an ARBITRARY station z instead of
    being hard-wired to the root (z=0).

    These two are deliberately NOT identical at z=0. strain_screen() is a
    GATE and uses the conservative bound |mx| + |vy|*L, which never
    under-estimates. This function is a COMPARISON quantity, so it uses
    the TRUE signed moment from moment_at_station() (mx - vy*(L-z)), which
    tracks the real value at rho = 0.9999 but is not a bound.

    Used for eps_mid_analytical (z = Z_MID), which main() compares against
    the post-run Abaqus eps_1d_mid_post and eps_3d_corner_post at the SAME
    station at the end of the batch. Comparing eps_root (root) against a
    Z_MID value would not be a fair check: the root usually has the larger
    bending strain."""
    mxz, myz = moment_at_station(vx, vy, mx, z)
    eps = (abs(mxz) + abs(myz)) * config.c / (config.E*config.I)
    gam = (abs(vx) + abs(vy)) / config.A / config.G
    return float(eps/2 + np.sqrt((eps/2)**2 + (gam/2)**2))

def bending_dominance(vx, vy, mx, z=None):
    """How strongly BENDING (curvature-driven axial/normal strain,
    varying with the lever arm to the free tip) dominates over
    TRANSVERSE SHEAR (constant along the whole beam) at station z —
    i.e. is this case still "thin-beam" (Euler-Bernoulli-like) behavior,
    or has it drifted toward "thick-beam"/shear-dominated behavior at
    the point where beam and solid are actually being compared.

    Returns eps_bend(z) / gam_shear, the same two ingredients
    strain_screen() combines via Mohr's circle, but reported as a RATIO
    instead of a combined principal strain, and evaluated at station z
    (default Z_MID, the comparison station) rather than only at the
    root. This matters because bending dominance is NOT worst (smallest)
    at the root — the lever arm (L - z), and therefore the bending
    contribution, is actually LARGEST at the root and shrinks toward the
    free tip (where a tip-applied moment alone remains but the lever arm
    vanishes) while shear stays constant everywhere. So Z_MID, being
    between root and tip, is a MORE demanding check than the root: a
    case can pass strain_screen()'s root-based bound comfortably while
    still being less bending-dominated at Z_MID than at the root.
    Returns +inf if there's no shear at all (gam == 0).

    NOTE: this is a DIAGNOSTIC value by default — see LOAD ACCEPTANCE
    LOGIC in the module docstring. It only becomes a real accept/reject
    gate when --min-bending-dominance is passed on the CLI."""
    if z is None:
        z = config.Z_MID
    mxz, myz = moment_at_station(vx, vy, mx, z)
    eps_bend = (abs(mxz) + abs(myz)) * config.c / (config.E*config.I)
    gam_shear = (abs(vx) + abs(vy)) / config.A / config.G
    if gam_shear <= 0:
        return float('inf')
    return float(eps_bend / gam_shear)

def range_safety_check(verbose=True):
    """EXACT (not Monte-Carlo) certification of whether the WHOLE sampling
    box (every (Vx,Vy,Mx) triple VX_RANGE x VY_RANGE x MX_RANGE can
    produce) already satisfies strain_screen()'s Mohr's-circle bound at
    BOTH ends of the beam -- root (z=0) and tip (z=L) -- the only two
    stations that can ever be governing, since |M(z)| is a convex
    (affine-in-z, absolute-valued) function of z and a convex function's
    max over an interval is always at an endpoint, never interior.

    Why this can be done exactly, with no sampling:
      TIP (z=L): the lever arm (L-z) vanishes, so eps_tip depends ONLY on
      Mx; gam depends ONLY on Vx,Vy. Disjoint variables -> their
      individual range-maxima can co-occur in one draw -> plugging the
      range extrema straight into the formula IS the tip's exact global
      worst case, not merely a conservative estimate of it.

      ROOT (z=0): mx_root = Mx - Vy*L involves Mx and Vy TOGETHER, but
      they're still independent variables, so |mx_root| is maximized by
      driving both to their range extrema with OPPOSING sign (always
      achievable -- neither constrains the other). The arithmetic below
      is unchanged by the sign fix: it already takes absolute values
      before combining, so (|Mx|max + |Vy|max*L) remains the exact
      worst case either way. my_root = Vx*L depends
      only on Vx, maximized independently at its own extreme, and since
      Vx doesn't appear in mx_root that choice doesn't conflict with the
      first. So both terms hit their max at ONE corner of the box
      simultaneously too -- same style of exact result as the tip.
      (This assumes ranges are symmetric about 0, e.g. VX_RANGE =
      (-4e-6, 4e-6), which is the current convention -- for a symmetric
      range the "matching sign" extremum is always reachable. For a
      one-sided/asymmetric range this formula is still a valid, safe
      UPPER bound, just not guaranteed exact.)

    A Monte Carlo pass-rate can never substitute for this: the exact
    worst-case corner is a single point in a continuous box, so random
    sampling has probability ~0 of ever landing on it, no matter how many
    draws you take or how close to 100% they pass.

    Run this any time you change VX_RANGE/VY_RANGE/MX_RANGE -- it's a
    handful of lines of arithmetic, not a simulation.

    Returns (root_safe, tip_safe): each True if that endpoint's exact
    worst-case principal strain stays under STRAIN_THRESHOLD."""
    def _abs_max(rng):
        return max(abs(rng[0]), abs(rng[1]))

    mx_max = _abs_max(config.MX_RANGE)
    vx_max = _abs_max(config.VX_RANGE)
    vy_max = _abs_max(config.VY_RANGE)
    gam_max = (vx_max + vy_max) / config.A / config.G

    def _principal(eps_max):
        return eps_max/2 + ((eps_max/2)**2 + (gam_max/2)**2)**0.5

    # ROOT: mx and vy combine (matching sign) in mx_root; vx alone drives
    # my_root. Both maximized simultaneously at one corner.
    eps_root_max = ((mx_max + vy_max*config.L) * config.c/(config.E*config.I)
                    + vx_max*config.L * config.c/(config.E*config.I))
    p_root_max = _principal(eps_root_max)
    root_margin = config.STRAIN_THRESHOLD - p_root_max
    root_safe = root_margin > 0

    # TIP: eps depends only on mx (lever arm vanishes at z=L).
    eps_tip_max = mx_max * config.c / (config.E * config.I)
    p_tip_max = _principal(eps_tip_max)
    tip_margin = config.STRAIN_THRESHOLD - p_tip_max
    tip_safe = tip_margin > 0

    if verbose:
        tprint('Range safety check (exact worst case over the WHOLE '
               'sampling box, not Monte Carlo):')
        tprint('  ROOT (z=0): p_root_max=%.4e vs T=%.4e -> %s (margin '
               '%.4e, %.1f%% of T)' % (
                   p_root_max, config.STRAIN_THRESHOLD,
                   'SAFE' if root_safe else 'UNSAFE', root_margin,
                   100.0*root_margin/config.STRAIN_THRESHOLD))
        tprint('  TIP  (z=L): p_tip_max=%.4e  vs T=%.4e -> %s (margin '
               '%.4e, %.1f%% of T)' % (
                   p_tip_max, config.STRAIN_THRESHOLD,
                   'SAFE' if tip_safe else 'UNSAFE', tip_margin,
                   100.0*tip_margin/config.STRAIN_THRESHOLD))
        if root_safe and tip_safe:
            tprint('  -> Both endpoints SAFE: EVERY (Vx,Vy,Mx) in the '
                   'current ranges is provably admissible. The pre-run '
                   'strain screen and the N_DRAW top-up margin are no '
                   'longer necessary for correctness (though harmless to '
                   'keep as a no-cost safety net).')
        else:
            bad_end = [] if root_safe else ['ROOT']
            bad_end += [] if tip_safe else ['TIP']
            tprint('  -> WARNING: %s breaches STRAIN_THRESHOLD at the '
                   'exact worst-case corner even though a per-case '
                   'screen could still discard that specific case. '
                   'Shrink the offending range(s), or keep relying on '
                   'the per-case screen (it is still doing real work).'
                   % ' and '.join(bad_end))
    return root_safe, tip_safe

def principal_strain_from_beam_row(row):
    """POST-run pass criterion. Max principal strain at a single station,
    built from Abaqus's OWN B31 section-strain output — SE1 (axial),
    SE2/SE3 (transverse shear), SK1/SK2 (bending curvature), SK3 (twist)
    — instead of the closed-form vx/vy/mx hand formula strain_screen()
    uses. Same combination logic (extreme-fiber bending strain =
    curvature * c, added to axial strain; shear from SE2/SE3; combined
    via the same worst-case Mohr's-circle construction as strain_screen),
    different inputs: this one is real FE output, evaluated after Abaqus
    has actually run, at whichever row of beam_whole.csv is passed in
    (call it with the row nearest Z=0 to check the clamped root).
    """
    se1 = float(row['SE1'])
    # Abaqus B31: SK1 and SK2 are the TWO BENDING curvatures; SK3 is the
    # TWIST rate. Confirmed independently: SM3/SK3 is constant to 8
    # significant figures and equals GJ (for the 1x1 section: implied
    # J = 0.14083 vs the closed-form 0.1406*a^4). SK3 therefore produces SHEAR
    # strain, not extreme-fibre normal strain, and must NOT be multiplied
    # by c and added into eps. (Previously this used SK2+SK3, which both
    # omitted the dominant bending term SK1 and wrongly included the
    # twist -- underestimating eps by a median factor of 2.4, worst 54x.)
    sk1, sk2 = float(row['SK1']), float(row['SK2'])
    se2, se3 = float(row['SE2']), float(row['SE3'])
    eps = abs(se1) + (abs(sk1) + abs(sk2)) * config.c
    gam = abs(se2) + abs(se3)
    return float(eps/2 + np.sqrt((eps/2)**2 + (gam/2)**2))

def rotated_beam_components(row, R_local):
    """[CHANGED 2026-10-08] Returns the B31 section strains/curvatures in the
    beam-local director basis AS ABAQUS REPORTS THEM (no rotation): SE/SK
    already are v0/k0 of Arora et al. The old behaviour, documented below,
    applied R_local.T, which rotated them twice. The *_rot / *_rot_global
    names are kept only for file compatibility; they now hold the
    director-basis values. Old description:

        sv_sp  = [SE3, SE2, SE1];  sv_rot = R_local.T @ sv_sp
        sk_sp  = [SK1, SK2, SK3];  sk_rot = R_local.T @ sk_sp
        axial         = sv_rot[2]   (from SE1)
        shearX        = sv_rot[0]   (from SE3)
        shearY        = sv_rot[1]   (from SE2)
        SK1_rot       = sk_rot[0]   (BENDING curvature)
        SK2_rot       = sk_rot[1]   (BENDING curvature)
        SK3_twist_rot = sk_rot[2]   (TWIST rate about the beam axis -- NOT a
                                     bending curvature; excluded from the
                                     extreme-fibre normal-strain term)

    This is the same transform used for sv_mat/sk_mat passed to FEniCS in
    batch_driver.run_case() (the pull-back described in warping_core's
    FRAMES note: sv_mat = R_local.T @ sv_spatial). Forces/moments, by
    contrast, use R_local (rotated_beam_forces()) and are passed to FEniCS
    unrotated as targets.

    Returns a dict of the individual components (floats). Used for logging
    (beam_1d_mid.csv, eps_1d_mid_post); not used by any admission gate."""
    sv_sp = np.array([float(row['SE3']), float(row['SE2']), float(row['SE1'])])
    sk_sp = np.array([float(row['SK1']), float(row['SK2']), float(row['SK3'])])
    # 2026-10-08: SE/SK are already director-basis (material) components;
    # no R_local.T pull-back (it was a double rotation). R_local kept in the
    # signature for compatibility. Keys/column names kept unchanged.
    sv_rot = sv_sp
    sk_rot = sk_sp
    return {
        'axial':          float(sv_rot[2]),
        'shearX':         float(sv_rot[0]),
        'shearY':         float(sv_rot[1]),
        'SK1_rot':        float(sk_rot[0]),   # BENDING curvature
        'SK2_rot':        float(sk_rot[1]),   # BENDING curvature
        'SK3_twist_rot':  float(sk_rot[2]),   # TWIST rate (beam axis)
    }

def rotated_beam_forces(row, R_local):
    """Rotate the B31 section forces/moments (SF1/SF2/SF3, SM1/SM2/SM3) from
    the beam-LOCAL frame into the GLOBAL frame with R_local (NO transpose):

        n_sp = [SF3, SF2, SF1];  n_rot = R_local @ n_sp
        m_sp = [SM1, SM2, SM3];  m_rot = R_local @ m_sp
        Fx_1d = n_rot[0];  Fy_1d = n_rot[1];  N_1d  = n_rot[2]
        Mx_1d = m_rot[0];  My_1d = m_rot[1];  Mz_1d = m_rot[2]  (torque)

    Axis "1" is the beam axis, "2"/"3" the section's transverse axes.

    VERIFIED direction (beam_mid.csv / test/test_200_summary.csv, Z=50,
    |UR| = 21.3 deg): R_local @ n gives Fz = 2.5e-8 (it must be ~0, no axial
    load is applied) and reproduces the FEniCS global Fx..Mz to 7 digits;
    R_local.T @ n gives Fz = 2.8e-6. Same direction is used for abq_n/abq_m
    in batch_driver.py and recompare.py."""
    n_sp = np.array([float(row['SF3']), float(row['SF2']), float(row['SF1'])])
    m_sp = np.array([float(row['SM1']), float(row['SM2']), float(row['SM3'])])
    n_rot = R_local @ n_sp
    m_rot = R_local @ m_sp
    return {
        'Fx_1d': float(n_rot[0]), 'Fy_1d': float(n_rot[1]), 'N_1d': float(n_rot[2]),
        'Mx_1d': float(m_rot[0]), 'My_1d': float(m_rot[1]), 'Mz_1d': float(m_rot[2]),
    }

def force_moment_1d3d_error(beam_row, R_local, resultant_row, min_ref=1e-9):
    """Compare the 1D beam's section forces/moments (rotated into the
    global frame via rotated_beam_forces()) against the 3D solid's
    NFORC-integrated resultants (solid_resultants.csv, at the SAME
    station) -- N/Vx/Vy/Mx/My, an independent numerical check on
    whether beam theory and the solved 3D stress field actually agree,
    same spirit as the existing 1D-vs-analytical / 3D-vs-analytical
    comparisons but now 1D-vs-3D directly, and at EVERY station (not
    just Z_MID).

    NOTE ON SIGN: this does NOT assume the two sides already share a
    sign convention -- SF/SM's sign convention and solid_resultants.csv
    NFORC's sign convention (internal reaction across a cut, see chat)
    were never verified to match a priori. A consistent ~200% relative
    error / opposite sign across every component is itself a
    diagnostic (points at a sign-convention mismatch, not necessarily a
    real modeling disagreement) -- don't just read the % blindly, check
    whether flipping one side's sign brings it back near 0%.

    Returns a dict with each pair (1d, 3d, diff, pct_diff) for N/Vx/Vy/
    Mx/My, using min_ref as the smallest denominator allowed (percent
    error is meaningless/unstable near a true zero)."""
    f1d = rotated_beam_forces(beam_row, R_local)
    pairs = {
        'N':  (f1d['N_1d'],  float(resultant_row['N'])),
        'Vx': (f1d['Fx_1d'], float(resultant_row['Vx'])),
        'Vy': (f1d['Fy_1d'], float(resultant_row['Vy'])),
        'Mx': (f1d['Mx_1d'], float(resultant_row['Mx'])),
        'My': (f1d['My_1d'], float(resultant_row['My'])),
    }
    # Mz (torque about the beam axis) -- only comparable if the 3D
    # NFORC extraction actually computed it (extract_odb.py's Mz column
    # was added after Mx/My/Vx/Vy/N; an OLDER solid_resultants.csv won't
    # have it). Skip the pair rather than crash/compare-against-nothing
    # if it's absent, same tolerance pattern used elsewhere in this file
    # for optional columns.
    if 'Mz' in resultant_row.index and resultant_row['Mz'] != '' \
            and not pd.isna(resultant_row['Mz']):
        pairs['Mz'] = (f1d['Mz_1d'], float(resultant_row['Mz']))
    # RESULTANT SCALES. A per-component percentage divides by that
    # component's own 3D value, which is meaningless wherever the true
    # value approaches zero -- and several of these do so BY CONSTRUCTION:
    # My = Vx*(L-z) and Mz both vanish at the free tip, and N is
    # identically zero because no axial load is applied. Measured at the
    # tip that yields "errors" of 124%, 115% and 4132% on quantities whose
    # correct value is zero, while the underlying vector agrees to a few
    # percent. Normalising each component by the magnitude of its own
    # resultant instead gives a denominator that cannot collapse, so the
    # number stays interpretable everywhere along the span.
    f3 = math.sqrt(sum(v3d * v3d for lab, (v1d, v3d) in pairs.items()
                       if lab in ('N', 'Vx', 'Vy')))
    m3 = math.sqrt(sum(v3d * v3d for lab, (v1d, v3d) in pairs.items()
                       if lab in ('Mx', 'My', 'Mz')))
    f3 = max(f3, min_ref)
    m3 = max(m3, min_ref)

    out = {}
    for label, (v1d, v3d) in pairs.items():
        diff = v1d - v3d
        ref = max(abs(v3d), min_ref)
        scale = f3 if label in ('N', 'Vx', 'Vy') else m3
        out['%s_1d' % label] = v1d
        out['%s_3d' % label] = v3d
        out['%s_diff' % label] = diff
        out['%s_pct_diff' % label] = 100.0 * diff / ref
        # preferred: cannot blow up where the component itself vanishes
        out['%s_pct_of_resultant' % label] = 100.0 * diff / scale
        # non-dimensionalised by the physical scales, matching the
        # *_norm convention already used in force_error.csv and
        # <case>_resultants.csv (F0 = EI/L^2, M0 = F0*L).
        sc0 = (config.FORCE_SCALE if label in ('N', 'Vx', 'Vy')
               else config.MOMENT_SCALE)
        out['%s_1d_norm' % label] = v1d / sc0
        out['%s_3d_norm' % label] = v3d / sc0
        out['%s_diff_norm' % label] = diff / sc0
    out['F_resultant_3d'] = f3
    out['M_resultant_3d'] = m3
    out['F_resultant_3d_norm'] = f3 / config.FORCE_SCALE
    out['M_resultant_3d_norm'] = m3 / config.MOMENT_SCALE
    return out

def write_force_moment_1d3d_csv(case_dir, beam_csv, resultants_csv, rotation_fn):
    """Loop over every Z station present in BOTH beam_whole.csv and
    solid_resultants.csv, compute force_moment_1d3d_error() at each,
    and write case_dir/force_moment_1d3d_error.csv (one row per
    station). rotation_fn(ur1, ur2, ur3) -> R_local is supplied by the
    caller (rotation_matrix_from_ur from fenics_solve.py in practice)
    rather than imported here, so this module keeps its documented
    zero-Abaqus/FEniCS/pyvista-dependency promise (see module
    docstring) -- callers that already have FEniCS loaded (batch_
    driver.py) or that load it just for this helper (reprocess_
    checks.py) both work without this file needing to know which.

    Returns the DataFrame written, or None if either CSV is missing/
    unusable."""
    try:
        df_b = pd.read_csv(beam_csv)
        df_b.columns = df_b.columns.str.strip()
        df_r = pd.read_csv(resultants_csv)
        df_r.columns = df_r.columns.str.strip()
    except Exception:
        return None
    needed_b = ['Z', 'UR1', 'UR2', 'UR3', 'SF1', 'SF2', 'SF3', 'SM1', 'SM2', 'SM3']
    needed_r = ['Z_snapped', 'N', 'Vx', 'Vy', 'Mx', 'My']
    if not all(c in df_b.columns for c in needed_b):
        return None
    if not all(c in df_r.columns for c in needed_r):
        return None

    rows = []
    for _, rrow in df_r.iterrows():
        z = float(rrow['Z_snapped'])
        idx = (df_b['Z'] - z).abs().idxmin()
        brow = df_b.loc[idx]
        if abs(float(brow['Z']) - z) > 1e-6 * max(config.L, 1.0):
            continue   # no matching beam station -- skip rather than mismatch
        ur1, ur2, ur3 = float(brow['UR1']), float(brow['UR2']), float(brow['UR3'])
        R_local = rotation_fn(ur1, ur2, ur3)
        err = force_moment_1d3d_error(brow, R_local, rrow)
        err['Z'] = z
        rows.append(err)

    if not rows:
        return None
    out_df = pd.DataFrame(rows).sort_values('Z')
    out_path = os.path.join(case_dir, 'force_moment_1d3d_error.csv')
    out_df.to_csv(out_path, index=False, float_format='%.8e')
    return out_df

def principal_strain_from_beam_row_rotated(row, R_local):
    """Same worst-case Mohr's-circle principal strain as
    principal_strain_from_beam_row(), but built from the components
    returned by rotated_beam_components() instead of the raw beam-local
    ones.

    Combines axial + extreme-fibre bending + shear into ONE worst-case
    number -- it is NOT the pure axial strain. Reported as eps_1d_mid_post
    for reference/comparison only; the admission gate (small_ok in
    strain_admit()) uses the UNROTATED max_beam_principal_strain() instead."""
    c = rotated_beam_components(row, R_local)
    eps = abs(c['axial']) + (abs(c['SK1_rot']) + abs(c['SK2_rot'])) * config.c
    gam = abs(c['shearX']) + abs(c['shearY'])
    return float(eps/2 + np.sqrt((eps/2)**2 + (gam/2)**2))

def bending_dominance_strain_from_beam_row(row):
    """REFERENCE ONLY (was the admission gate up to rule version 2).

        eps_bend  = (|SK1| + |SK2|) * c     extreme-fibre bending strain
        gam_shear = |SE2| + |SE3|           section-average shear strain
        returns eps_bend / gam_shear  (+inf if gam_shear ~ 0)

    A STRAIN ratio, ignores axial strain SE1 and torsion SK3. Kept so it
    can still be logged next to the stress-based check below."""
    sk1, sk2 = float(row['SK1']), float(row['SK2'])
    se2, se3 = float(row['SE2']), float(row['SE3'])
    eps_bend  = (abs(sk1) + abs(sk2)) * config.c
    gam_shear = abs(se2) + abs(se3)
    if gam_shear <= 1e-20:
        return float('inf')
    return float(eps_bend / gam_shear)


# Saint-Venant torsion of a square a x a section: tau_max = T / (0.208 a^3)
_TORSION_COEF_SQUARE = 0.208

def bending_dominance_from_beam_row(row):
    """ADMISSION GATE (rule version 3): is the peak AXIAL (normal) stress
    larger than the peak SHEAR stress at this beam station?

    Uses the Abaqus B31 section forces/moments of the row (the local
    section frame -- magnitudes, so no rotation needed):
        SF1 = axial force N,   SF2/SF3 = transverse shears,
        SM1/SM2 = bending moments,   SM3 = torque.

        sigma_max = |SF1|/A + (|SM1| + |SM2|) * c / I
                    (axial + both bending stresses at the worst corner)
        tau_max   = 1.5 * sqrt(SF2^2 + SF3^2) / A        (parabolic shear)
                  + |SM3| / (0.208 * a^3)                 (square torsion)
        returns sigma_max / tau_max  (+inf if tau_max ~ 0)

    Both are upper-bound peaks (the maxima do not sit at the same point),
    so the ratio is conservative. bending_ok = ratio >= k with
    k = config.MIN_BENDING_DOMINANCE_POST (1.0 -> "axial stress is the
    largest"). Falls back to the old strain ratio if the row has no
    SF/SM columns (very old extractions)."""
    keys = ('SF1', 'SF2', 'SF3', 'SM1', 'SM2', 'SM3')
    try:
        sf1, sf2, sf3, sm1, sm2, sm3 = [float(row[k]) for k in keys]
    except (KeyError, TypeError, ValueError):
        return bending_dominance_strain_from_beam_row(row)
    a = config.SECTION_SIDE
    sigma_max = abs(sf1) / config.A + (abs(sm1) + abs(sm2)) * config.c / config.I
    tau_max = (1.5 * math.sqrt(sf2**2 + sf3**2) / config.A
               + abs(sm3) / (_TORSION_COEF_SQUARE * a**3))
    if tau_max <= 1e-30:
        return float('inf')
    return float(sigma_max / tau_max)

def _solid_corner_rows(solid_csv):
    """Locate the extreme-fiber CORNER nodes in solid_mid.csv -- the
    nodes farthest from the section centroid (max sqrt(X^2+Y^2)), i.e.
    the (+/-Xmax, +/-Ymax) corners of a rectangular section. Robust to
    exact mesh coordinates (doesn't assume nodes sit at precisely
    config.c) by taking whichever nodes actually have the maximum
    radius present in the mesh, with a small relative tolerance so all
    ~4 corners are picked up even with float noise.

    Returns a DataFrame of the corner row(s) (usually 4, one per
    corner, at whatever Z is nearest config.Z_MID if solid_mid.csv ever
    has multiple Z's). Raises ValueError if solid_csv is empty/missing
    columns."""
    df = pd.read_csv(solid_csv)
    df.columns = df.columns.str.strip()
    if not {'X', 'Y', 'Z'}.issubset(df.columns) or len(df) == 0:
        raise ValueError('solid_csv missing X/Y/Z columns or empty: %s' % solid_csv)
    if df['Z'].nunique() > 1:
        z_pick = (df['Z'] - config.Z_MID).abs().idxmin()
        df = df[df['Z'] == df.loc[z_pick, 'Z']]
    r = np.sqrt(df['X'].values**2 + df['Y'].values**2)
    r_max = r.max()
    if r_max <= 0:
        raise ValueError('All nodes at r=0 in %s -- no corner to find' % solid_csv)
    mask = r >= 0.999 * r_max
    return df.loc[mask]

def principal_strain_from_solid_extreme_fiber(solid_csv):
    """Principal strain of the 3D solid at the section's EXTREME-FIBRE CORNER
    nodes at Z_MID (from solid_mid.csv), using LE33 (normal) and
    |LE13|+|LE23| (engineering shear) in the Mohr's-circle form. Returns
    the max over the corners found.

    REFERENCE ONLY -- not used for admission (strain_admit() gates on the
    1D beam alone). Reported as eps_3d_corner_post and compared with
    eps_mid_analytical at the end of the batch; the corner is the same
    physical point the beam-theory extreme-fibre prediction refers to.
    LE components are in GLOBAL axes, so under large section rotation this
    is not the section-frame strain.

    Returns (max_eps, x, y) for the governing corner, or (None, None, None)
    if solid_csv is unusable or missing LE columns."""
    try:
        corners = _solid_corner_rows(solid_csv)
    except Exception:
        return None, None, None
    if not {'LE33', 'LE13', 'LE23'}.issubset(corners.columns):
        return None, None, None
    best_eps, best_x, best_y = None, None, None
    for _, row in corners.iterrows():
        eps = abs(float(row['LE33']))
        gam = abs(float(row['LE13'])) + abs(float(row['LE23']))
        p = float(eps/2 + np.sqrt((eps/2)**2 + (gam/2)**2))
        if best_eps is None or p > best_eps:
            best_eps, best_x, best_y = p, float(row['X']), float(row['Y'])
    return best_eps, best_x, best_y

def principal_strain_from_solid_centerline(cl_row):
    """Principal strain of the 3D solid at the section centreline node
    (X=0, Y=0, located by _solid_centerline_row()), from LE33 and
    |LE13|+|LE23| in the Mohr's-circle form.

    REFERENCE ONLY (eps_3d_centreline_post). The centreline sits at/near the
    neutral axis, where bending strain is ~0, so this is NOT the worst 3D
    strain and is not comparable to extreme-fibre values."""
    eps = abs(float(cl_row['LE33']))
    gam = abs(float(cl_row['LE13'])) + abs(float(cl_row['LE23']))
    return float(eps/2 + np.sqrt((eps/2)**2 + (gam/2)**2))

def _solid_centerline_row(solid_csv, tol=1e-9):
    """Locate the exact X=0,Y=0 centreline node in solid_mid.csv (exact
    match within tol, not nearest-neighbour). If several Z layers are
    present, the one nearest Z_MID is used. Raises ValueError if no such
    node exists."""
    df = pd.read_csv(solid_csv)
    df.columns = df.columns.str.strip()
    mask = (df['X'].abs() < tol) & (df['Y'].abs() < tol)
    if not mask.any():
        raise ValueError('No node at X=0,Y=0 (tol=%.1e) in %s' % (tol, solid_csv))
    matches = df[mask]
    if len(matches) > 1:
        idx = (matches['Z'] - config.Z_MID).abs().idxmin()
        return matches.loc[idx]
    return matches.iloc[0]

def max_beam_principal_strain(beam_csv):
    """Max 1D principal strain over EVERY row of beam_whole.csv, using the
    RAW (beam-local, unrotated) principal_strain_from_beam_row() at each
    row. This is the value the small_ok admission gate uses.

    beam_whole.csv holds one row per Z_STATIONS entry passed to
    extract_odb.py -- currently 0, 10, ..., 100 (see _z_stations_arg() in
    batch_driver.py), so the root, Z_MID and the tip are all covered.
    Older case folders extracted before Z_STATIONS existed may hold only
    the Z_MID row; n_rows_scanned in the return value shows which.

    Returns (max_eps, z_at_max, n_rows_scanned) -- (None, None, 0) if
    beam_csv is unreadable or missing the needed columns."""
    try:
        df = pd.read_csv(beam_csv)
    except Exception:
        return None, None, 0
    df.columns = df.columns.str.strip()
    needed = ['Z', 'SE1', 'SE2', 'SE3', 'SK1', 'SK2', 'SK3']
    if not all(col in df.columns for col in needed) or len(df) == 0:
        return None, None, 0
    eps_vals = df.apply(principal_strain_from_beam_row, axis=1)
    idx = eps_vals.idxmax()
    return float(eps_vals.loc[idx]), float(df.loc[idx, 'Z']), int(len(df))

def strain_admit(beam_row_mid, beam_csv, solid_csv, R_local,
                  max_post_strain, bd_mid_abaqus=None,
                  min_bending_dominance_post=None):
    """The post-Abaqus admission check. A case is admitted
    (admit_all=True) only if BOTH hold:

      1. bending_ok -- bending_dominance_mid_abaqus = sigma_max/tau_max
         from the real Abaqus 1D SF/SM at Z_MID (via
         bending_dominance_from_beam_row(); axial+bending vs
         shear+torsion stress) is >= min_bending_dominance_post: the
         peak axial stress is the largest stress in the ACTUAL solved
         model. The old strain ratio is returned as
         bending_dominance_strain_mid_abaqus for reference only.
      2. small_ok   -- the maximum principal strain ANYWHERE beam_csv
         has a row (via max_beam_principal_strain(), NOT just the single
         Z_MID row) is under max_post_strain. 1D ONLY -- the 3D
         strains (eps_3d_corner_post / eps_3d_centreline_post) are
         returned for reference but are NOT gated on. The "anywhere in
         the beam" part deliberately does NOT use the combined Z_MID-
         only eps_1d_mid_post -- see max_beam_principal_strain()'s docstring
         for the caveat that this is only as thorough as what beam_csv
         actually contains.

    REMOVED: the former 1D/3D axial-agreement check (agree_ok /
    centerline_strain_agreement()) -- was a diagnostic comparison
    between beam-theory axial strain and the solid's true local LE33,
    never actually validating strain magnitude or bending/shear
    character, and is no longer part of admission. axial_1d_post is
    still returned for reference/logging only.

    Returns a dict with eps_1d_mid_post/eps_3d_centreline_post (full Z_MID principal
    strain, kept for reference/logging), max_eps_1d_beam/z_at_max_eps_1d
    /n_beam_rows_scanned (the beam-wide max actually used by small_ok's
    1D side), axial_1d_post (pure axial, reference only), bd_mid_abaqus,
    the two booleans, and admit_all. Never raises on a missing/
    unreadable solid_csv -- eps_3d_centreline_post comes back None (not applicable)
    instead."""
    beam_rot = rotated_beam_components(beam_row_mid, R_local)
    axial_1d = beam_rot['axial']
    eps_1d = principal_strain_from_beam_row_rotated(beam_row_mid, R_local)
    max_eps_1d_beam, z_at_max_eps_1d, n_beam_rows = max_beam_principal_strain(beam_csv)
    if max_eps_1d_beam is None:
        # beam_csv unreadable/missing columns -- fall back to the single
        # Z_MID row rather than silently passing an unknown max.
        max_eps_1d_beam, z_at_max_eps_1d, n_beam_rows = eps_1d, None, 0

    eps_3d = None          # centerline — reference only, not gated on
    eps_3d_corner = None   # extreme-fiber — reference only, not gated on
    corner_x = corner_y = None
    try:
        cl_row = _solid_centerline_row(solid_csv)
        eps_3d = principal_strain_from_solid_centerline(cl_row)
    except Exception:
        pass
    eps_3d_corner, corner_x, corner_y = principal_strain_from_solid_extreme_fiber(solid_csv)

    if bd_mid_abaqus is None:
        try:
            bd_mid_abaqus = bending_dominance_from_beam_row(beam_row_mid)
        except Exception:
            bd_mid_abaqus = None
    try:
        bd_strain_ref = bending_dominance_strain_from_beam_row(beam_row_mid)
    except Exception:
        bd_strain_ref = None
    bd_known = bd_mid_abaqus is not None and not (
        isinstance(bd_mid_abaqus, float) and np.isnan(bd_mid_abaqus))
    thresh = (min_bending_dominance_post if min_bending_dominance_post
              is not None else config.MIN_BENDING_DOMINANCE_POST)
    bending_ok = bool(bd_known and bd_mid_abaqus >= thresh) if bd_known else None

    # 1D-only gate: the 3D strains above are still computed and returned
    # for reference/logging, but no longer affect admission.
    small_ok = bool(max_eps_1d_beam < max_post_strain)

    checks = [bending_ok, small_ok]
    admit_all = all(c is not None and c for c in checks)

    return {
        'eps_1d_mid_post': eps_1d,
        'eps_3d_centreline_post': eps_3d,                 # centerline, reference only
        'eps_3d_corner_post': eps_3d_corner,    # extreme-fiber, reference only
        'corner_x': corner_x, 'corner_y': corner_y,
        'max_eps_1d_beam': max_eps_1d_beam,
        'z_at_max_eps_1d': z_at_max_eps_1d,
        'n_beam_rows_scanned': n_beam_rows,
        'axial_1d_post': axial_1d,
        # 'bd_mid_abaqus' removed: it was a byte-identical duplicate of
        # bending_dominance_mid_abaqus, which load_existing() also sets.
        'bending_dominance_mid_abaqus': bd_mid_abaqus,   # stress ratio (gate)
        'bending_dominance_strain_mid_abaqus': bd_strain_ref,  # old, reference
        'bending_ok': bending_ok,
        'small_ok': small_ok,
        'admit_all': admit_all,
    }

def write_beam_1d_row_csv(case_dir, row, R_local, eps_1d_mid_post_rotated=None):
    """Write case_dir/beam_1d_mid.csv: a one-row dump of the Abaqus 1D
    (beam) SE/SK/UR/SF/SM values at the Z_MID row used for the post-run
    checks, so each case folder records what went into them.

    Contains the raw beam-local values exactly as Abaqus reported them,
    plus the R_local.T-transformed values from rotated_beam_components()
    (*_rot_global columns) and eps_1d_mid_post_rotated.

    Raises KeyError if a required column is missing; callers wrap it in
    try/except and log a warning."""
    rot = rotated_beam_components(row, R_local)
    out = {
        'Z': float(row['Z']),
        # raw beam-local (as Abaqus reported them)
        'SE1': float(row['SE1']), 'SE2': float(row['SE2']), 'SE3': float(row['SE3']),
        'SK1': float(row['SK1']), 'SK2': float(row['SK2']), 'SK3': float(row['SK3']),
        'UR1': float(row['UR1']), 'UR2': float(row['UR2']), 'UR3': float(row['UR3']),
        'SF1': float(row.get('SF1', 0.0)), 'SF2': float(row.get('SF2', 0.0)),
        'SF3': float(row.get('SF3', 0.0)),
        'SM1': float(row.get('SM1', 0.0)), 'SM2': float(row.get('SM2', 0.0)),
        'SM3': float(row.get('SM3', 0.0)),
        # director-basis values (no rotation since 2026-10-08; names kept)
        'axial_rot_global':  rot['axial'],    # rotated SE1, <-> LE33
        'shearX_rot_global': rot['shearX'],   # rotated SE3, <-> LE13
        'shearY_rot_global': rot['shearY'],   # rotated SE2, <-> LE23
        'SK1_bend_rot_global':  rot['SK1_rot'],        # rotated SK1 (BENDING)
        'SK2_bend_rot_global':  rot['SK2_rot'],        # rotated SK2 (BENDING)
        'SK3_twist_rot_global': rot['SK3_twist_rot'],  # rotated SK3 (TWIST)
        'eps_1d_mid_post_rotated': eps_1d_mid_post_rotated,
    }
    pd.DataFrame([out]).to_csv(os.path.join(case_dir, 'beam_1d_mid.csv'),
                               index=False, float_format='%.10e')


# ════════════════════════════════════════════════════════════════════════════
# Stress invariants
# ════════════════════════════════════════════════════════════════════════════

def invariants(s):
    """I1 (hydrostatic), Mises (~J2), J3 (determinant-type third invariant)
    — the standard practical invariant triplet used across yield/failure
    criteria (Bai-Wierzbicki, Drucker-Prager, ductile fracture literature,
    etc.):
      - I1 = tr(sigma), the FULL stress tensor's first invariant.
        Controls pressure-sensitivity / triaxiality (hydrostatic part).
      - Mises = sqrt(3*J2), J2 the DEVIATORIC tensor's second invariant.
        Controls distortion energy (the "size" of the stress state).
      - J3 = det(dev(sigma)), the DEVIATORIC tensor's third invariant —
        distinguishes WHERE on the yield surface the state sits
        (tension-meridian vs compression-meridian) at the same Mises
        level: J3 > 0 for tensile-dominated states, J3 < 0 for
        compression-dominated states, J3 = 0 for pure shear.

    Replaced the old Lode-angle formula (theta = arccos(cos3theta)/3):
    J3 needs no arccos/clipping, no divide-by-J2 guard for the J2~0
    hydrostatic-state case (det() is defined everywhere), and it keeps
    I1/J2/J3 all from a consistent convention rather than mixing in
    I3 = det(full sigma) (a different, less common invariant). Being a
    determinant of a stress tensor, J3 has stress^3 units, so it needs
    its OWN physical scale (config.J3_SCALE, not STRESS_SCALE).
    """
    s11,s22,s33 = s['S11'],s['S22'],s['S33']
    s12,s13,s23 = s['S12'],s['S13'],s['S23']
    I1 = s11+s22+s33
    mises = np.sqrt(0.5*((s11-s22)**2+(s22-s33)**2+(s33-s11)**2)
                    +3*(s12**2+s13**2+s23**2))
    # deviatoric diagonal components (off-diagonals unchanged by trace shift)
    d11, d22, d33 = s11-I1/3, s22-I1/3, s33-I1/3
    # J3 = det(dev(sigma)), via cofactor expansion of the deviatoric
    # stress tensor [[d11,s12,s13],[s12,d22,s23],[s13,s23,d33]].
    J3 = (d11*(d22*d33-s23**2) - s12*(s12*d33-s23*s13)
          + s13*(s12*s23-d22*s13))
    return {'Mises':mises, 'inv_1':I1, 'J3':J3}

def _root_strain_and_mid_bd_from_beam_csv(beam_csv):
    """Load beam_whole.csv and return (eps_root_abaqus, z_root_used,
    bd_mid_abaqus) -- the ROOT-strain row and the Z_MID-bending-dominance
    row, which are DIFFERENT rows and must be selected independently.

    FIXED (was `_mid_strain_and_bd_from_beam_csv`): the old version
    picked ONE row (nearest Z=0) and used it for BOTH eps AND bd, then
    called the bd result "bd_mid_abaqus" even though it was evaluated at
    the root, not Z_MID. This disagreed with the fresh-run path in
    batch_driver.py's run_case(), which has always correctly used two
    separate rows (idx_root for eps, idx nearest config.Z_MID for bd) --
    so a case's bending-dominance number could differ depending on
    whether it was freshly run or reloaded from disk. This version picks
    each row independently, matching run_case()'s logic exactly.

    REQUIRES beam_whole.csv to actually contain a row near Z=0 (root) AND
    a row near config.Z_MID -- i.e. Z_STATIONS must include both when
    extract_odb.py is called (see _z_stations_arg() in batch_driver.py).
    If beam_whole.csv only has ONE row (e.g. an older run, or Z_STATIONS
    wasn't passed), idxmin/idx-nearest both degenerate to that same row,
    and z_root_used will NOT actually be near 0 -- check z_root_used
    against 0.0 if you need to confirm this ran correctly.

    Returns (None, None, None) if beam_whole.csv is unavailable or missing
    required columns (e.g. an older run)."""
    try:
        df = pd.read_csv(beam_csv)
    except Exception:
        return None, None, None
    df.columns = df.columns.str.strip()
    needed = ['Z','SE1','SE2','SE3','SK1','SK2','SK3']
    if not all(col in df.columns for col in needed) or len(df) == 0:
        return None, None, None

    idx_root = df['Z'].abs().idxmin()
    row_root = df.iloc[idx_root]
    z_root_used = float(row_root['Z'])
    eps = principal_strain_from_beam_row(row_root)
    if abs(z_root_used) > 1e-6 * max(config.L, 1.0):
        tprint('WARNING: %s has no row near Z=0 -- z_root_used=%.4f is '
               'not the true root (beam_whole.csv likely only has the '
               'Z_MID station). eps_root_abaqus below is NOT a root '
               'check.' % (beam_csv, z_root_used))

    idx_mid = (df['Z'] - config.Z_MID).abs().idxmin()
    row_mid = df.iloc[idx_mid]
    try:
        bd = bending_dominance_from_beam_row(row_mid)
    except Exception:
        bd = None
    return eps, z_root_used, bd


# Backward-compat alias -- old name, now delegates to the fixed function.
# Remove once nothing else imports the old name.
def _mid_strain_and_bd_from_beam_csv(beam_csv):
    return _root_strain_and_mid_bd_from_beam_csv(beam_csv)


# ── added 7 Oct 2026 ─────────────────────────────────────────────────────────
def mid_moment_nd(vx, vy, mx):
    """Non-dimensional bending moment at Z_MID, |m_x| + |m_y| in M0 units,
    from the tip loads: Mx(z) = Mx_tip - Vy (L - z), My(z) = Vx (L - z).
    Equal at every section size for equal non-dimensional loads. Small
    values mean the moment crosses zero near mid-span, where relative
    stress/energy errors blow up."""
    arm = config.L - config.Z_MID
    return (abs(mx - vy * arm) + abs(vx * arm)) / config.MOMENT_SCALE


def read_case_diagnostics(case_dir):
    """Scalars from convergence_1d.csv and section_frame_check.csv, if
    present, for summary.csv. Never raises."""
    import os as _os
    res = {}
    for fn in ('convergence_1d.csv', 'section_frame_check.csv',
               'solid_energy_history.csv', 'shear_strain_correction.csv'):
        p = _os.path.join(case_dir, fn)
        if _os.path.isfile(p):
            try:
                row = pd.read_csv(p).iloc[0].to_dict()
                res.update({k: v for k, v in row.items()})
            except Exception:
                pass
    if res.get('ALLIE'):
        try:
            res['ALLAE_over_ALLIE'] = abs(float(res.get('ALLAE', float('nan')))) / abs(float(res['ALLIE']))
        except Exception:
            pass
    return res


# ── gate-table pre-run checks (added 7 Oct 2026) ─────────────────────────────
# Written in non-dimensional loads m = M/M0, v = V/F0, valid for any section
# size: bending strain = m c / L, and transverse shear / bending stress at a
# station = 1.5 (v F0 / A) / (m M0 c / I) = v h / (4 L m).
def _nd_loads(vx, vy, mx):
    return (vx / config.FORCE_SCALE, vy / config.FORCE_SCALE,
            mx / config.MOMENT_SCALE)


def table_strain_gate(vx, vy, mx):
    """Max bending strain along the span: (c/L) max(|m_tip - v_y| + |v_x|,
    |m_tip|). |m(xi)| is linear in xi, so the max is at the root or the tip."""
    vxn, vyn, m = _nd_loads(vx, vy, mx)
    return (config.c / config.L) * max(abs(m - vyn) + abs(vxn), abs(m))


def table_bending_ratio(vx, vy, mx, z=None):
    """Transverse shear / bending stress at station z (default Z_MID):
    sqrt(vx^2+vy^2) h / (4 L (|m_x(z)| + |m_y(z)|)). = v/(400 m) for h=1,
    L=100. inf if the moment is exactly zero."""
    z = config.Z_MID if z is None else z
    vxn, vyn, m = _nd_loads(vx, vy, mx)
    arm = (config.L - z) / config.L
    mm = abs(m - vyn * arm) + abs(vxn * arm)
    v = (vxn**2 + vyn**2) ** 0.5
    return float('inf') if mm == 0 else v * config.SECTION_SIDE / (4.0 * config.L * mm)

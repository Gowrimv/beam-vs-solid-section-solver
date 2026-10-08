"""
batch_driver.py  —  WSL-native
================================
Run from WSL:  python3 batch_driver.py [--cases 1000] [--seed 42] [--dry-run]

Optimisations vs previous version:
  1. FEniCS runs in-process (mesh + problem built once, reused every case)
  2. 1D and 3D Abaqus jobs run simultaneously per case (Popen + wait)
  3. Up to MAX_PARALLEL cases run concurrently via ThreadPoolExecutor
  4. Completed cases are detected and skipped on resume

Architecture:
  - Everything runs in WSL Python 3 (FEniCS, pyvista, pandas all available)
  - Abaqus called via /mnt/c/... path — works from WSL
  - abaqus python extract_odb.py  — Windows subprocess, no licence consumed
  - FEniCS forward_solve + postprocess — direct function call, no subprocess
  - compare.py's in-process logic — pyvista available in WSL

Token budget:  MAX_PARALLEL=3 cases × 2 Abaqus jobs × 5 tokens = 30 tokens
               (you have 38 available — safe margin of 8)

FILE LAYOUT
----------------------------------------------------------------------
  config.py         — paths, section/material constants, sampling ranges,
                       the eq.9-12 non-dimensional scales, and the shared
                       tprint() logging helper. Stdlib only — nothing
                       here ever needs Abaqus, FEniCS, or pyvista, so
                       ANY other module can import it for free.
  plotting.py        — the three end-of-batch summary plots
                       (histograms/energy-scatter/strain-screen-audit).
                       Depends on nothing but matplotlib/pandas/numpy/
                       config — no pyvista, no dolfin — so it can
                       regenerate plots from an existing summary.csv
                       without re-running the batch at all.
  plot_results.py    — standalone CLI built on plotting.py: `python3
                       plot_results.py --out-dir <cases dir>` regenerates
                       histograms.png/energy_scatter.png/
                       strain_screen_audit.png from that directory's
                       summary.csv alone, in seconds, with no Abaqus/
                       FEniCS/pyvista involved.
  strain_checks.py   — pure numpy/pandas post-processing math: the
                       pre-run analytical (Mohr's-circle) screen, the
                       post-run real-Abaqus strain/bending-dominance
                       checks, the frame-rotation helpers, the
                       post-Abaqus admission check (strain_admit), and the
                       stress invariants used by compare.py. No Abaqus/
                       FEniCS/pyvista dependency, same reasoning as
                       plotting.py — reused by run_case(), load_existing(),
                       AND compare.py, so it needs one home rather than
                       being duplicated or creating an import cycle.
  abaqus_io.py       — Abaqus/WSL job I/O: path translation, .inp
                       editing, job submission (stdin starved, hard
                       timeout), and the leftover-output-file cleanup.
                       Stdlib only (os/re/subprocess/threading).
  fenics_solve.py    — the in-process FEniCS warping solve: builds the
                       mesh/problem ONCE (init_fenics(), from main())
                       and re-solves it per case (solve_and_postprocess())
                       under a module-level lock, since ThreadPoolExecutor
                       workers share this process's memory. Needs
                       FEniCS/dolfin — one of the two import-heavy
                       modules (with compare.py).
  compare.py         — the in-process FEniCS-vs-Abaqus-3D comparison
                       stage: stress invariants/errors.csv, the six-
                       component contour plots, the stress AND (new)
                       strain centerline overlays, and force_error.csv.
                       Needs pyvista + matplotlib.
  batch_driver.py    — (this file) orchestration only now: the
                       resume/reload logic (load_existing), the per-case
                       pipeline (run_case — Abaqus I/O, extraction,
                       FEniCS solve, comparison, and the post-run gates,
                       all delegated to the modules above), and the
                       CLI/screening/top-up loop (main()).

LOAD ACCEPTANCE LOGIC — what actually gates a case vs what only reports
------------------------------------------------------------------------
A random (Vx,Vy,Mx) draw passes these checkpoints before it becomes a
row in summary.csv. Every discard triggers a top-up re-draw.

  1. PRE-RUN, ALWAYS ON — strain_screen(vx,vy,mx):
     Closed-form Mohr's-circle principal strain at the clamped ROOT,
     using the conservative moment bound |Mx| + |Vy|*L. A draw with
     strain_screen(...) > STRAIN_THRESHOLD is discarded before any
     Abaqus job is submitted.

  2. PRE-RUN, OFF BY DEFAULT — bending_dominance(vx,vy,mx,Z_MID):
     Closed-form bending/shear ratio at Z_MID. Only rejects a draw if
     --min-bending-dominance VALUE is passed. Otherwise
     bending_dominance_mid/_root are just written to summary.csv.

  3. RUN-TIME, ALWAYS ON — Abaqus / extraction / FEniCS failure:
     run_case() returns its OWN reason string for each failure type
     (DISCARD_ABAQUS_FAILED, DISCARD_EXTRACT_FAILED,
     DISCARD_FENICS_NOT_BUILT, DISCARD_FENICS_DIVERGED,
     DISCARD_FENICS_ERROR -- see the constants below the imports).
     main() counts and logs each reason separately.

  4. POST-ABAQUS, PRE-FEniCS, ALWAYS ON — strain_admit()
     (strain_checks.py), from the Abaqus 1D BEAM output only:
       bending_ok -- bending_dominance_mid_abaqus = sigma_max/tau_max
                     at Z_MID from 1D SF/SM:
                     sigma_max = |SF1|/A + (|SM1|+|SM2|)c/I,
                     tau_max = 1.5|V|/A + |SM3|/(0.208 a^3)
                     >= --min-bending-dominance-post
                     (default config.MIN_BENDING_DOMINANCE_POST).
       small_ok   -- max 1D principal strain over EVERY beam_whole.csv
                     station (0, 10, ..., L) < --max-post-strain
                     (default 5e-3).
     A failing case is dropped BEFORE the FEniCS solve, returns
     DISCARD_REJECTED_ADMISSION, and gets a REJECTED_admission.json
     marker (thresholds + ADMISSION_RULE_VERSION) so a rerun with the
     same settings skips it (DISCARD_REJECTED_PREV). The 3D solid
     strains (eps_3d_corner_post, eps_3d_centreline_post) are computed
     and reported but NOT gated on.

  5. POST-RUN, ALWAYS ON — root_strain_under_threshold (in main()):
     eps_root_abaqus = 1D principal strain at the beam_whole.csv row
     nearest Z=0 (the real clamped root, since Z_STATIONS includes 0),
     must be <= STRAIN_THRESHOLD. Normally already implied by small_ok;
     it only bites if --max-post-strain is set above STRAIN_THRESHOLD.
     main() also re-checks bending_dominance_mid_abaqus and admit_all,
     so a case reloaded from disk (load_existing()) is re-gated exactly
     like a fresh one. Failures count as n_post_gate_failed.

  REMOVED: the former 1D/3D axial-agreement check (25% limit).

  Summary table:
      check                                    default  gate flag
      ---------------------------------------  -------  ---------------------------
      strain_screen (root, analytical)         ON       (always)
      bending_dominance (Z_MID, theory)        OFF      --min-bending-dominance
      Abaqus / extract / FEniCS failure        ON       (always)
      bending_ok (Z_MID, Abaqus 1D)            ON       --min-bending-dominance-post
      small_ok (all stations, Abaqus 1D)       ON       --max-post-strain
      root_strain_under_threshold (Z=0, 1D)    ON       (always, STRAIN_THRESHOLD)

RESUME LOGIC — what actually makes run_case() treat a case as "already
done" vs redo it
------------------------------------------------------------------------
This is a SEPARATE decision from LOAD ACCEPTANCE LOGIC above (that's
about whether a case's numbers are good; this is about whether run_case()
does any work at all for a given case_dir on this invocation). On a
plain run with NO flags (--redo-fenics/--redo-compare both off, the
normal case), exactly two things are checked, independently, and only
one of them can produce a full skip:

  A. odbs_exist = both <tag>_1d.odb and <tag>_3d.odb are present.
     This ONLY skips the Abaqus call itself (write_inp + run_abaqus_case).
     It does NOT skip extract_odb.py, the FEniCS solve, or the
     comparison stage — those still run unconditionally whenever this is
     the only thing satisfied. In the log this looks like:
         [NNNN] ODBs exist — skipping Abaqus
     immediately followed by a full re-extract/re-solve/re-compare for
     that case. Having the .odb files (or beam_whole.csv, solid_mid.csv,
     *_resultants.csv, *_ref.vtu — none of those are checked here) is
     NOT sufficient to skip the rest of the pipeline.

  B. os.path.isfile(case_dir/'errors.csv') — this is the ONLY check that
     can skip a case ENTIRELY (Abaqus + extraction + FEniCS + comparison
     all skipped) and go straight to load_existing(). Nothing else in
     the case directory is consulted for this decision: not the ODBs,
     not beam_whole.csv/solid_mid.csv, not *_resultants.csv/*_ref.vtu, not
     summary.csv. If errors.csv is missing for any reason (an earlier
     run was interrupted or raised partway through _compare_inprocess,
     which writes errors.csv itself partway through that function —
     contours.png/centerline_overlay.png/force_error.csv are written
     AFTER it, so a crash between those points can leave errors.csv
     either present-but-file-list-incomplete or absent-while-later-stage-
     files-exist depending on exactly where it died), the case is
     redone from wherever odbs_exist leaves off, with no other file
     given credit as evidence of a completed run.

  KNOWN GAP — the "fully complete" message is printed BEFORE the reload
  is verified to have worked:
      if os.path.isfile(errors.csv) and not dry and not redo_compare:
          tprint('[%04d] fully complete — reloading' % case_id)
          result = load_existing(case_dir)
          if result:                      # <- can be falsy/None
              ...
              return result
      # falls through here if result was empty/None — NO further message
      # explains why; the case is silently redone in full despite having
      # just been logged as "fully complete."

CHANGE LOG
---------------------------
1. Added a POST-run pass criterion (`principal_strain_from_beam_row` /
   `root_strain_under_threshold`).
2. Stresses and section forces/moments in the "post" outputs reported as
   dimensionless fractions of physical scales defined in config.py.
3. Pulled config/plotting/plot_results out into separate modules.
4. Made both POST-run checks (4)/(5) above into ALWAYS-ON hard gates.
5. Fixed frame bug in `centerline_strain_agreement()`.
6. Changed solid-side centerline node lookup to exact match on X=0,Y=0.
7. extract_odb.py now receives a Z_STATIONS argument (comma-separated
   list) so beam_whole.csv has one row per station instead of one row
   total, and solid_resultants.csv is written with NFORC-integrated
   cross-section resultants at every station. _z_stations_arg() builds
   the station string; it is called directly in the run() call inside
   run_case() and once in main() for logging — no changes to run_case()'s
   signature or ex.submit() call.
8. Admission small-strain gate is 1D-only (3D strains reference only);
   unused MAX_1D3D_DIFF_PCT_DEFAULT / SKIP_AGREE_CHECK_DEFAULT removed.
9. run_case() returns a distinct discard reason per failure type instead
   of the catch-all 'diverged'; REJECTED_admission.json carries
   ADMISSION_RULE_VERSION.
10. extract_odb.py section-point shear no longer doubled (Abaqus E12/E13
    are already engineering shear strains).
"""
import argparse
import os
import glob
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

import config
from config import tprint
from plotting import plot_histograms, plot_energy_scatter, plot_strain_screen_audit

from abaqus_io import (win, abq, write_inp, run, clean_job_files,
                       run_abaqus_case, job_completed, residual_check_1d)
from fenics_solve import init_fenics, solve_and_postprocess, rotation_matrix_from_ur, fenics_ready
from loads_io import export_accepted_loads, load_fixed_loads
from strain_checks import (
    MAX_POST_STRAIN_DEFAULT,
    strain_screen, strain_screen_at_station, bending_dominance,
    range_safety_check, principal_strain_from_beam_row,
    bending_dominance_from_beam_row,
    strain_admit, write_force_moment_1d3d_csv,
    write_beam_1d_row_csv, _root_strain_and_mid_bd_from_beam_csv,
    mid_moment_nd, read_case_diagnostics,
    table_strain_gate, table_bending_ratio,
)
from compare import _compare_inprocess


# ════════════════════════════════════════════════════════════════════════════
# run_case() discard reasons. run_case() returns one of these strings
# (instead of the old catch-all 'diverged') when a case is dropped, so
# main() can report WHY each case was discarded and count each reason
# separately. Any other return value is a result dict (or None on a hard
# error).
# ════════════════════════════════════════════════════════════════════════════
DISCARD_ABAQUS_FAILED      = 'abaqus_failed'       # job failed or timed out
DISCARD_EXTRACT_FAILED     = 'extract_failed'      # extract_odb.py rc != 0
DISCARD_FENICS_NOT_BUILT   = 'fenics_not_built'    # init_fenics() never ran
DISCARD_FENICS_DIVERGED    = 'fenics_diverged'     # FEniCS non-convergence
DISCARD_FENICS_ERROR       = 'fenics_error'        # other exception in the
                                                   # post-extract/FEniCS block
DISCARD_REJECTED_ADMISSION = 'rejected_admission'  # strain_admit() failed
DISCARD_HOURGLASS          = 'hourglass_energy'    # ALLAE/ALLIE > --max-hourglass (added 7 Oct 2026)
DISCARD_UNCONVERGED_1D     = 'unconverged_1d'      # 1D residual > --max-1d-residual (added 7 Oct 2026)
DISCARD_REJECTED_PREV      = 'rejected_previously' # REJECTED_admission.json
                                                   # exists, same thresholds
# Bump whenever strain_admit()'s rules change. REJECTED_admission.json
# markers written under a different (or missing) version are ignored and
# the case is re-evaluated. v2 = 1D-only small-strain gate (3D dropped).
ADMISSION_RULE_VERSION = 3   # 3: stress-based sigma_max/tau_max bending check

DISCARD_REASONS = (DISCARD_ABAQUS_FAILED, DISCARD_EXTRACT_FAILED,
                   DISCARD_FENICS_NOT_BUILT, DISCARD_FENICS_DIVERGED,
                   DISCARD_FENICS_ERROR, DISCARD_UNCONVERGED_1D, DISCARD_HOURGLASS, DISCARD_REJECTED_ADMISSION,
                   DISCARD_REJECTED_PREV)
DISCARD_LABELS = {
    DISCARD_ABAQUS_FAILED:      'Abaqus failed/timed out',
    DISCARD_EXTRACT_FAILED:     'ODB extraction failed',
    DISCARD_FENICS_NOT_BUILT:   'FEniCS problem not built',
    DISCARD_FENICS_DIVERGED:    'FEniCS did not converge',
    DISCARD_FENICS_ERROR:       'FEniCS/post-extract error',
    DISCARD_UNCONVERGED_1D:     '1D Abaqus not converged (residual > limit)',
    DISCARD_HOURGLASS:          '3D hourglass energy too large (ALLAE/ALLIE)',
    DISCARD_REJECTED_ADMISSION: 'rejected at admission (axial-stress dominance/strain)',
    DISCARD_REJECTED_PREV:      'rejected at admission (earlier run)',
}


# ════════════════════════════════════════════════════════════════════════════
# CHANGE 7 — Z-station string builder
# Called directly inside run_case()'s run() call and once in main() for
# logging. No changes to run_case()'s signature or ex.submit().
# ════════════════════════════════════════════════════════════════════════════

# Files extract_odb.py must produce; a case is only kept if all exist.
_EXTRACT_REQUIRED = ('beam_whole.csv', 'solid_mid.csv', 'solid_resultants.csv')

# Stop the batch instead of burning through case IDs when the disk is
# full or every case is failing the same way.
MIN_FREE_GB            = 5.0   # abort before a new round if less is free
MAX_CONSECUTIVE_HARD   = 10    # abort after this many hard errors in a row

def _free_gb(path):
    import shutil
    try:
        return shutil.disk_usage(path).free / 1024.0**3
    except OSError:
        return float('inf')

def _z_stations_arg(L, step=10.0):
    """Build the comma-separated Z_STATIONS string passed to extract_odb.py
    as argv[5]. Produces evenly-spaced stations from 0 to L inclusive,
    always including config.Z_MID even when it doesn't land on a multiple
    of step."""
    stations = []
    z = 0.0
    while z <= L + 1e-9:
        stations.append(round(z, 6))
        z += step
    if config.Z_MID not in stations:
        stations.append(round(config.Z_MID, 6))
        stations.sort()
    return ','.join('%.6g' % z for z in stations)


# ════════════════════════════════════════════════════════════════════════════
# Per-case pipeline
# ════════════════════════════════════════════════════════════════════════════

def check_section_matches_inp(inp_path=None):
    """Start-up guard: read the *Beam Section, SECTION=RECT dimensions from
    the 1D template .inp and compare them with config.SECTION_SIDE. A
    mismatch means config.py (A, I, c, load ranges, scales) describes a
    different section than Abaqus will actually run -- every check and
    normalisation would silently be wrong. Returns True (match), False
    (mismatch) or None (couldn't read/parse -- e.g. dry run off-machine).
    The 3D .inp and the FEniCS mesh are NOT checked here; keep them in
    step by hand (see README, section 7)."""
    inp_path = inp_path or config.INP_1D
    try:
        with open(inp_path) as f:
            lines = f.read().splitlines()
    except Exception:
        tprint('Section check: could not read %s -- skipped.' % inp_path)
        return None
    for i, ln in enumerate(lines):
        low = ln.strip().lower().replace(' ', '')
        if low.startswith('*beamsection') and 'section=rect' in low:
            for nxt in lines[i+1:]:
                s = nxt.strip()
                if not s or s.startswith('**'):
                    continue
                try:
                    dims = [float(t) for t in s.split(',')[:2]]
                except ValueError:
                    break
                ok = all(abs(d - config.SECTION_SIDE) <= 1e-9 * max(1.0, d)
                         for d in dims)
                tprint('Section check: %s RECT %s vs config.SECTION_SIDE=%g -> %s'
                       % (os.path.basename(inp_path), dims,
                          config.SECTION_SIDE, 'OK' if ok else 'MISMATCH'))
                return ok
            break
    tprint('Section check: no *Beam Section, SECTION=RECT found in %s -- '
           'skipped.' % inp_path)
    return None


def load_existing(case_dir):
    """Re-read scalars from a previously completed case (the errors.csv-
    exists resume path in run_case()). Re-derives everything main()'s
    post-run gates need (eps_root_abaqus / root_strain_under_threshold,
    bending_dominance_mid_abaqus, and the full strain_admit() block, with
    R_local rebuilt from the Z_MID row's own UR1/UR2/UR3), so a reloaded
    case is re-gated exactly like a freshly-run one. Also rewrites
    beam_1d_mid.csv. Returns None if nothing could be read."""
    result = {}
    err_path = os.path.join(case_dir, 'errors.csv')
    if os.path.isfile(err_path):
        edf = pd.read_csv(err_path).set_index('component')
        for comp in ['Mises','inv_1','J3','S33','S13','S23','N_resultant','M_resultant','Energy']:
            if comp in edf.index:
                result['relL2_'+comp] = float(edf.loc[comp,'relL2'])
                result['L2_'+comp]    = float(edf.loc[comp,'L2'])
                # accept either spelling: errors.csv files written before the
                # rename still carry 'L2_normalized'. Same quantity either way
                # (the L2 error divided by the physical scale), so read whichever
                # is present and always WRITE the unambiguous new name.
                _sc = ('L2_over_scale' if 'L2_over_scale' in edf.columns
                       else ('L2_normalized' if 'L2_normalized' in edf.columns else None))
                if _sc is not None:
                    result['L2_over_scale_' + comp] = float(edf.loc[comp, _sc])
                if 'relL2_resultant' in edf.columns:
                    result['relL2_resultant_' + comp] = float(edf.loc[comp, 'relL2_resultant'])
                if 'relL2_Szz' in edf.columns:
                    result['relL2_Szz_' + comp] = float(edf.loc[comp, 'relL2_Szz'])
    energy_path = os.path.join(case_dir, 'solid_energy.csv')
    if os.path.isfile(energy_path):
        result['ELSE_per_length'] = float(
            pd.read_csv(energy_path)['U_per_length_ELSE'].iloc[0])
    resultants_path = os.path.join(case_dir, glob.glob(
        os.path.join(case_dir, '*_resultants.csv'))[0]) if glob.glob(
        os.path.join(case_dir, '*_resultants.csv')) else None
    beam_csv = os.path.join(case_dir, 'beam_whole.csv')
    if os.path.isfile(beam_csv):
        eps_abq, z_used, bd_abq = _root_strain_and_mid_bd_from_beam_csv(beam_csv)
        if eps_abq is not None:
            result['eps_root_abaqus'] = eps_abq
            result['z_root_used'] = z_used
            result['root_strain_under_threshold'] = eps_abq <= config.STRAIN_THRESHOLD
        if bd_abq is not None:
            result['bending_dominance_mid_abaqus'] = bd_abq

        solid_csv_cl = os.path.join(case_dir, 'solid_mid.csv')
        if os.path.isfile(solid_csv_cl):
            try:
                df_b = pd.read_csv(beam_csv)
                df_b.columns = df_b.columns.str.strip()
                needed = ['Z','SE1','SE2','SE3','UR1','UR2','UR3']
                if all(col in df_b.columns for col in needed) and len(df_b):
                    idx = (df_b['Z'] - config.Z_MID).abs().idxmin()
                    row = df_b.iloc[idx]
                    ur1, ur2, ur3 = (float(row['UR1']), float(row['UR2']),
                                     float(row['UR3']))
                    result['UR1'], result['UR2'], result['UR3'] = ur1, ur2, ur3
                    R_local = rotation_matrix_from_ur(ur1, ur2, ur3)
                    result.update(strain_admit(
                        row, beam_csv, solid_csv_cl, R_local,
                        MAX_POST_STRAIN_DEFAULT, bd_mid_abaqus=bd_abq))
                    write_beam_1d_row_csv(case_dir, row, R_local,
                                          result.get('eps_1d_mid_post'))
            except Exception as ex:
                tprint('WARNING: %s: strain admission check failed on '
                       'reload: %s' % (case_dir, ex))
    if result:
        result.update(read_case_diagnostics(case_dir))
    return result if result else None


def run_case(case_id, vx, vy, mx, case_dir, dry=False, redo_fenics=False,
             redo_compare=False, max_post_strain=MAX_POST_STRAIN_DEFAULT,
             min_bending_dominance_post=None, max_1d_residual=0.01,
             max_hourglass=config.HOURGLASS_MAX):
    os.makedirs(case_dir, exist_ok=True)
    tag = 'case_%04d' % case_id

    # ── REJECTED-case marker: a case that failed the pre-FEniCS admission
    # check (strain3['admit_all'] is False) never produces *_resultants.csv
    # or *_ref.vtu, since solve_and_postprocess() is never reached. Without
    # this marker, --redo-compare's use_redo_compare check below always
    # sees those files missing and falls back to a full run, which just
    # re-derives the SAME rejection from scratch (re-running extract_odb.py
    # and the admission check pointlessly) -- and a plain rerun (no
    # --redo-compare) does the same thing implicitly, since errors.csv was
    # never written either (see RESUME LOGIC note above). Check this BEFORE
    # any of that logic runs, and skip straight to DISCARD_REJECTED_PREV if the
    # rejection thresholds haven't changed since it was written -- if they
    # HAVE changed (e.g. a looser --max-post-strain this run), fall through
    # and let the case be genuinely re-evaluated, since it might pass now.
    rejected_marker = os.path.join(case_dir, 'REJECTED_admission.json')
    if os.path.isfile(rejected_marker) and not dry and not redo_fenics:
        try:
            import json
            with open(rejected_marker) as f:
                prev = json.load(f)
            same_thresholds = (
                prev.get('max_post_strain') == max_post_strain and
                prev.get('min_bending_dominance_post') == min_bending_dominance_post and
                prev.get('admission_rule') == ADMISSION_RULE_VERSION)
            if same_thresholds:
                tprint('[%04d] previously rejected at admission (bending_ok=%s '
                       'small_ok=%s, same thresholds) — skipping re-derivation'
                       % (case_id, prev.get('bending_ok'), prev.get('small_ok')))
                return DISCARD_REJECTED_PREV
            else:
                tprint('[%04d] previously rejected, but thresholds or '
                       'admission rules changed — re-evaluating' % case_id)
        except Exception:
            pass  # unreadable/corrupt marker -- fall through, re-evaluate normally

    # ── redo_fenics: wipe FEniCS/comparison outputs only, keep the .odb ────
    if redo_fenics and not dry:
        for pat in ['errors.csv', 'force_error.csv', 'contours.png',
                    'centerline_overlay.png', '*_ref.vtu', '*_deformed.vtu',
                    '*_resultants.csv', '*_section_X0.csv',
                    '*_Stress_centerline.png', '*_Strain_centerline.png',
                    'REJECTED_admission.json']:
            for p in glob.glob(os.path.join(case_dir, pat)):
                try: os.remove(p)
                except OSError: pass
        tprint('[%04d] redo-fenics: cleared FEniCS/comparison outputs, '
               'keeping .odb files' % case_id)

    inp1d = os.path.join(case_dir, tag+'_1d.inp')
    inp3d = os.path.join(case_dir, tag+'_3d.inp')
    odb1d_wsl = os.path.join(case_dir, tag+'_1d.odb')
    odb3d_wsl = os.path.join(case_dir, tag+'_3d.odb')
    odbs_exist = (os.path.isfile(odb1d_wsl) and os.path.isfile(odb3d_wsl)
                  and job_completed(case_dir, tag+'_1d')
                  and job_completed(case_dir, tag+'_3d'))

    beam_csv        = os.path.join(case_dir, 'beam_whole.csv')
    solid_csv       = os.path.join(case_dir, 'solid_mid.csv')
    resultants_path = os.path.join(case_dir, tag + '_resultants.csv')
    ref_vtus        = glob.glob(os.path.join(case_dir, tag + '_ref.vtu'))
    use_redo_compare = (redo_compare and not dry
                         and os.path.isfile(beam_csv)
                         and os.path.isfile(solid_csv)
                         and os.path.isfile(resultants_path)
                         and len(ref_vtus) > 0)
    if redo_compare and not dry and not use_redo_compare:
        tprint('[%04d] --redo-compare: prior outputs incomplete (need '
               'beam_whole.csv, solid_mid.csv, %s_resultants.csv, '
               '%s_ref.vtu) — this case needs a normal run (or '
               '--redo-fenics) first. Falling back to a full run.'
               % (case_id, tag, tag))

    if (os.path.isfile(os.path.join(case_dir, 'errors.csv')) and not dry
            and not redo_compare):
        tprint('[%04d] fully complete — reloading' % case_id)
        result = load_existing(case_dir)
        if result:
            result.update({'case':case_id,'Vx':vx,'Vy':vy,'Mx':mx,
                           'eps_root':strain_screen(vx,vy,mx),
                           'eps_mid_analytical':
                               strain_screen_at_station(vx, vy, mx, config.Z_MID)})
            return result

    eps_root_abaqus = None
    z_root_used = None
    strain_agree = {}
    strain3 = {}
    bd_mid_abaqus = None
    ur_vals = (None, None, None)
    R_local = None

    if use_redo_compare:
        # ── FAST PATH: --redo-compare ────────────────────────────────────
        tprint('[%04d] --redo-compare: skipping Abaqus/extract/FEniCS-solve, '
               'recomputing comparison only' % case_id)
        df_b = pd.read_csv(beam_csv)
        df_b.columns = df_b.columns.str.strip()

        if all(col in df_b.columns for col in
               ['Z','SE1','SE2','SE3','SK1','SK2','SK3']) and len(df_b):
            idx_root = df_b['Z'].abs().idxmin()
            row_root = df_b.iloc[idx_root]
            z_root_used = float(row_root['Z'])
            eps_root_abaqus = principal_strain_from_beam_row(row_root)
            if eps_root_abaqus > config.STRAIN_THRESHOLD:
                tprint('[%04d] POST-SCREEN: Abaqus root strain %.3e > '
                       'threshold %.3e at Z=%.3f (pre-run analytical '
                       'screen said %.3e)' % (
                    case_id, eps_root_abaqus, config.STRAIN_THRESHOLD,
                    z_root_used, strain_screen(vx, vy, mx)))

        idx = (df_b['Z'] - config.Z_MID).abs().idxmin()
        row = df_b.iloc[idx]

        try:
            bd_mid_abaqus = bending_dominance_from_beam_row(row)
            bd_mid_theory = bending_dominance(vx, vy, mx, config.Z_MID)
            tprint('[%04d] bending_dominance_mid: theory=%.4g  '
                   'abaqus=%.4g' % (case_id, bd_mid_theory, bd_mid_abaqus))
        except Exception as ex:
            tprint('[%04d] WARNING: bending_dominance_from_beam_row failed: '
                   '%s' % (case_id, ex))

        ur1, ur2, ur3 = float(row['UR1']), float(row['UR2']), float(row['UR3'])
        ur_vals = (ur1, ur2, ur3)
        R_local = rotation_matrix_from_ur(ur1, ur2, ur3)

        _solid_csv_cl = os.path.join(case_dir, 'solid_mid.csv')

        strain3 = strain_admit(
            row, beam_csv, _solid_csv_cl, R_local, max_post_strain,
            bd_mid_abaqus=bd_mid_abaqus,
            min_bending_dominance_post=min_bending_dominance_post)
        tprint('[%04d] admission check: max_eps_1d_beam=%.3e (n_rows=%d, z=%s)  '
               'eps_3d_ref=%s  bending_ok=%s  small_ok(1D)=%s  -> admit=%s'
               % (case_id,
                  strain3['max_eps_1d_beam'], strain3['n_beam_rows_scanned'],
                  ('%.2f' % strain3['z_at_max_eps_1d']) if strain3['z_at_max_eps_1d'] is not None else 'n/a',
                  ('%.3e' % strain3['eps_3d_corner_post']) if strain3['eps_3d_corner_post'] is not None else 'n/a',
                  strain3['bending_ok'], strain3['small_ok'],
                  strain3['admit_all']))
        if strain3['n_beam_rows_scanned'] <= 1:
            tprint('[%04d] NOTE: beam_whole.csv only has %d row(s) -- '
                   'max_eps_1d_beam is really just the Z_MID value, not a '
                   'true whole-beam maximum.'
                   % (case_id, strain3['n_beam_rows_scanned']))
        try:
            write_beam_1d_row_csv(case_dir, row, R_local, strain3['eps_1d_mid_post'])
        except Exception as ex:
            tprint('[%04d] WARNING: write_beam_1d_row_csv failed: %s' % (case_id, ex))

        sf1 = float(row.get('SF1', 0.0)); sf2 = float(row.get('SF2', 0.0))
        sf3 = float(row.get('SF3', 0.0)); sm1 = float(row.get('SM1', 0.0))
        sm2 = float(row.get('SM2', 0.0)); sm3 = float(row.get('SM3', 0.0))
        n_t_sp = np.array([sf3, sf2, sf1])
        m_t_sp = np.array([sm1, sm2, sm3])
        abq_n = R_local @ n_t_sp
        abq_m = R_local @ m_t_sp

        res_row = pd.read_csv(resultants_path).iloc[0]
        Phi   = float(res_row['Phi'])
        n_out = np.array([float(res_row['n1']), float(res_row['n2']), float(res_row['n3'])])
        m_out = np.array([float(res_row['m1']), float(res_row['m2']), float(res_row['m3'])])

    elif not dry:
        # 1. write INPs and run Abaqus (skipped if ODBs already exist)
        if odbs_exist:
            tprint('[%04d] ODBs exist — skipping Abaqus' % case_id)
        else:
            write_inp(config.INP_1D, inp1d, vx, vy, mx)
            write_inp(config.INP_3D, inp3d, vx, vy, mx)

            job1d = tag+'_1d';  job3d = tag+'_3d'
            clean_job_files(case_dir, job1d)
            clean_job_files(case_dir, job3d)

            tprint('[%04d] launching Abaqus 1D + 3D simultaneously' % case_id)
            _ta = time.time()
            ok1, ok2 = run_abaqus_case(case_id, job1d, inp1d, job3d, inp3d, case_dir)
            tprint('[%04d] timing: Abaqus 1D+3D %.0fs' % (case_id, time.time()-_ta))
            if not (ok1 and ok2):
                tprint('[%04d] Abaqus failed or timed out (1D ok=%s, 3D ok=%s) '
                       '— discarding case' % (case_id, ok1, ok2))
                return DISCARD_ABAQUS_FAILED

        # 2. extract from ODBs
        # CHANGE 7: pass Z_STATIONS as argv[5] to extract_odb.py so it
        # snaps beam nodes at every station and writes solid_resultants.csv.
        # _z_stations_arg() is called here directly — no signature change
        # to run_case() or ex.submit().
        odb1d = win(odb1d_wsl)
        odb3d = win(odb3d_wsl)
        _extract_outputs = [os.path.join(case_dir, f) for f in _EXTRACT_REQUIRED]
        if not dry:
            for f in _extract_outputs:
                if os.path.exists(f):
                    try: os.remove(f)
                    except OSError: pass
        _te = time.time()
        rc = run(abq('python', win(config.EXTRACT_SCRIPT), odb1d, odb3d,
                     win(case_dir), str(config.Z_MID),
                     _z_stations_arg(config.L, step=10.0)),   # <-- only new arg
                 cwd=case_dir, dry=dry)
        # abaqus python via cmd.exe can return rc=0 after a traceback
        # (seen: corrupt 3D ODB -> OdbError, rc 0), so also require every
        # expected output file to have been (re)written.
        tprint('[%04d] timing: extract_odb %.0fs' % (case_id, time.time()-_te))
        missing = [os.path.basename(f) for f in _extract_outputs
                   if not os.path.isfile(f)] if not dry else []
        if rc != 0 or missing:
            tprint('[%04d] extract_odb failed (rc=%s, missing=%s) — '
                   'discarding case' % (case_id, rc, missing or 'none'))
            return DISCARD_EXTRACT_FAILED

        # 1D-vs-3D force/moment resultant cross-check, every station --
        # see force_moment_1d3d_error()/write_force_moment_1d3d_csv() in
        # strain_checks.py. Independent of the strain checks below;
        # never discards a case, purely diagnostic output.
        if not dry:
            solid_resultants_csv = os.path.join(case_dir, 'solid_resultants.csv')
            try:
                fm_df = write_force_moment_1d3d_csv(
                    case_dir, beam_csv, solid_resultants_csv,
                    rotation_matrix_from_ur)
                if fm_df is not None:
                    mz_str = (' Mz=%.1f' % fm_df['Mz_pct_diff'].abs().median()) \
                             if 'Mz_pct_diff' in fm_df.columns else ' Mz=n/a'
                    tprint('[%04d] force_moment_1d3d_error.csv: %d station(s), '
                           'median |%%diff|: N=%.1f Vx=%.1f Vy=%.1f Mx=%.1f My=%.1f%s'
                           % (case_id, len(fm_df),
                              fm_df['N_pct_diff'].abs().median(),
                              fm_df['Vx_pct_diff'].abs().median(),
                              fm_df['Vy_pct_diff'].abs().median(),
                              fm_df['Mx_pct_diff'].abs().median(),
                              fm_df['My_pct_diff'].abs().median(),
                              mz_str))
                else:
                    tprint('[%04d] force_moment_1d3d_error.csv: skipped '
                           '(no matching stations / missing columns)' % case_id)
            except Exception as ex:
                tprint('[%04d] WARNING: force/moment 1D-vs-3D check failed: '
                       '%s' % (case_id, ex))

        # extract_freeBody.py removed (see chat: built on Abaqus scripting
        # API calls -- session.FreeBodyFromOdb, this viewport pattern --
        # that don't match the real interface and can't be fixed without
        # documentation this session doesn't have reliable access to.
        # solid_resultants.csv (NFORC-based, in extract_odb.py) already
        # covers the same data and is working.)

        # ── 1D convergence check (added 7 Oct 2026) ─────────────────────
        # Reject a case whose 1D Abaqus solution was accepted with a large
        # residual (Abaqus's zero-force fallback). Written to
        # convergence_1d.csv either way so summary.csv can report it.
        conv1d = residual_check_1d(case_dir, tag + '_1d', vx, vy, mx, config.L)
        try:
            pd.DataFrame([conv1d]).to_csv(
                os.path.join(case_dir, 'convergence_1d.csv'), index=False)
        except Exception:
            pass
        if (max_1d_residual is not None and conv1d.get('res_ratio_1d') is not None
                and conv1d['res_ratio_1d'] > max_1d_residual):
            tprint('[%04d] 1D Abaqus NOT converged: residual = %.3g of applied '
                   'load (limit %.3g, time-avg fallback=%s) -- discarding case'
                   % (case_id, conv1d['res_ratio_1d'], max_1d_residual,
                      conv1d.get('tavg_fallback_1d')))
            return DISCARD_UNCONVERGED_1D

        # ── 3D hourglass-energy check (added 7 Oct 2026) ────────────────
        try:
            _eh = os.path.join(case_dir, 'solid_energy_history.csv')
            if os.path.isfile(_eh):
                _e = pd.read_csv(_eh).iloc[0]
                if 'ALLAE' in _e and 'ALLIE' in _e and abs(float(_e['ALLIE'])) > 0:
                    _hg = abs(float(_e['ALLAE'])) / abs(float(_e['ALLIE']))
                    if max_hourglass is not None and _hg > max_hourglass:
                        tprint('[%04d] 3D hourglass energy ALLAE/ALLIE = %.3g > %.3g '
                               '-- discarding case' % (case_id, _hg, max_hourglass))
                        return DISCARD_HOURGLASS
                else:
                    tprint('[%04d] WARNING: ALLAE/ALLIE not in solid_energy_history.csv '
                           '-- hourglass check skipped' % case_id)
            else:
                tprint('[%04d] WARNING: no solid_energy_history.csv -- hourglass '
                       'check skipped' % case_id)
        except Exception as ex:
            tprint('[%04d] WARNING: hourglass check failed: %s' % (case_id, ex))

        if not fenics_ready():
            tprint('[%04d] FEniCS problem was never built — discarding.'
                   % case_id)
            return DISCARD_FENICS_NOT_BUILT

        try:
            df_b = pd.read_csv(beam_csv)
            df_b.columns = df_b.columns.str.strip()

            if all(col in df_b.columns for col in
                   ['Z','SE1','SE2','SE3','SK1','SK2','SK3']) and len(df_b):
                idx_root = df_b['Z'].abs().idxmin()
                row_root = df_b.iloc[idx_root]
                z_root_used = float(row_root['Z'])
                eps_root_abaqus = principal_strain_from_beam_row(row_root)
                if eps_root_abaqus > config.STRAIN_THRESHOLD:
                    tprint('[%04d] POST-SCREEN: Abaqus root strain %.3e > '
                           'threshold %.3e at Z=%.3f (pre-run analytical '
                           'screen said %.3e)' % (
                        case_id, eps_root_abaqus, config.STRAIN_THRESHOLD,
                        z_root_used, strain_screen(vx, vy, mx)))

            idx  = (df_b['Z'] - config.Z_MID).abs().idxmin()
            row  = df_b.iloc[idx]

            try:
                bd_mid_abaqus = bending_dominance_from_beam_row(row)
                bd_mid_theory = bending_dominance(vx, vy, mx, config.Z_MID)
                tprint('[%04d] bending_dominance_mid: theory=%.4g  '
                       'abaqus=%.4g' % (case_id, bd_mid_theory, bd_mid_abaqus))
            except Exception as ex:
                tprint('[%04d] WARNING: bending_dominance_from_beam_row '
                       'failed: %s' % (case_id, ex))

            # ── EARLY REJECTION before FEniCS ────────────────────────────
            # Run the admission check here so a shear-dominated or
            # strain-too-large case never reaches solve_and_postprocess.
            sv_sp = np.array([float(row['SE3']),
                              float(row['SE2']),
                              float(row['SE1'])])
            # ── shear strains consistent with the section (7 Oct 2026) ──
            # FEniCS gets both strains and forces; Abaqus's softened default
            # shear stiffness makes SE2/SE3 ~25% too large for the section,
            # which shows up as spurious in-plane stress. Recompute them from
            # the shear forces: SE = SF / (k G A). Axial SE1 is untouched.
            shear_corr = (np.nan, np.nan)
            if getattr(config, 'CORRECT_SHEAR_STRAINS', False):
                _K = config.SHEAR_K_FACTOR * config.G * config.A
                _se3_new = float(row['SF3']) / _K
                _se2_new = float(row['SF2']) / _K
                shear_corr = (
                    (float(row['SE3']) / _se3_new) if abs(_se3_new) > 0 else np.nan,
                    (float(row['SE2']) / _se2_new) if abs(_se2_new) > 0 else np.nan)
                sv_sp = np.array([_se3_new, _se2_new, float(row['SE1'])])
                tprint('[%04d] shear strains for FEniCS: SE3 %.4e -> %.4e, '
                       'SE2 %.4e -> %.4e (k=%.4f)'
                       % (case_id, float(row['SE3']), _se3_new,
                          float(row['SE2']), _se2_new, config.SHEAR_K_FACTOR))
                try:
                    pd.DataFrame([{'shear_k_factor': config.SHEAR_K_FACTOR,
                                   'SE3_abaqus_over_used': shear_corr[0],
                                   'SE2_abaqus_over_used': shear_corr[1]}]).to_csv(
                        os.path.join(case_dir, 'shear_strain_correction.csv'), index=False)
                except Exception:
                    pass
            sk_sp = np.array([float(row['SK1']),
                              float(row['SK2']),
                              float(row['SK3'])])
            ur1,ur2,ur3 = float(row['UR1']),float(row['UR2']),float(row['UR3'])
            ur_vals = (ur1, ur2, ur3)
            ur_fen  = np.array([ur1, ur2, ur3])
            R_local = rotation_matrix_from_ur(*ur_fen)
            # Abaqus B31 SE/SK are already components in the beam-local (director)
            # basis d_i = R e_i, i.e. exactly v0 = R^T v, k0 = R^T k of Arora et al.
            # (2019). Pulling back again with R_local.T was a double rotation
            # (fixed 2026-10-08, see PARAMETERS_AND_CHANGES.md 2.19). R_local is
            # still used on the OUTPUT side (postprocess, abq_n/abq_m).
            sv_mat  = sv_sp
            sk_mat  = sk_sp
            R_MAP   = np.eye(3)

            _solid_csv_cl = os.path.join(case_dir, 'solid_mid.csv')

            strain3 = strain_admit(
                row, beam_csv, _solid_csv_cl, R_local, max_post_strain,
                bd_mid_abaqus=bd_mid_abaqus,
                min_bending_dominance_post=min_bending_dominance_post)
            tprint('[%04d] admission check: max_eps_1d_beam=%.3e (n_rows=%d, z=%s)  '
                   'eps_3d_ref=%s  bending_ok=%s  small_ok(1D)=%s  -> admit=%s'
                   % (case_id,
                      strain3['max_eps_1d_beam'], strain3['n_beam_rows_scanned'],
                      ('%.2f' % strain3['z_at_max_eps_1d']) if strain3['z_at_max_eps_1d'] is not None else 'n/a',
                      ('%.3e' % strain3['eps_3d_corner_post']) if strain3['eps_3d_corner_post'] is not None else 'n/a',
                      strain3['bending_ok'], strain3['small_ok'],
                      strain3['admit_all']))
            if strain3['n_beam_rows_scanned'] <= 1:
                tprint('[%04d] NOTE: beam_whole.csv only has %d row(s) -- '
                       'max_eps_1d_beam is really just the Z_MID value.'
                       % (case_id, strain3['n_beam_rows_scanned']))
            try:
                write_beam_1d_row_csv(case_dir, row, R_local, strain3['eps_1d_mid_post'])
            except Exception as ex:
                tprint('[%04d] WARNING: write_beam_1d_row_csv failed: %s' % (case_id, ex))

            # reject before FEniCS if the admission check failed
            if strain3.get('admit_all') is False:
                tprint('[%04d] pre-FEniCS rejection: bending_ok=%s  '
                       'small_ok=%s — skipping FEniCS solve'
                       % (case_id, strain3['bending_ok'], strain3['small_ok']))
                try:
                    import json
                    with open(rejected_marker, 'w') as f:
                        json.dump({
                            'bending_ok': strain3.get('bending_ok'),
                            'small_ok': strain3.get('small_ok'),
                            'max_post_strain': max_post_strain,
                            'min_bending_dominance_post': min_bending_dominance_post,
                            'admission_rule': ADMISSION_RULE_VERSION,
                        }, f)
                except Exception as ex:
                    tprint('[%04d] WARNING: failed to write REJECTED_admission.json: '
                           '%s' % (case_id, ex))
                return DISCARD_REJECTED_ADMISSION

            tprint('[%04d] row Z=%.6f  SE1=%.10e SE2=%.10e SE3=%.10e  '
                   'SK1=%.10e SK2=%.10e SK3=%.10e' % (
                case_id, float(row['Z']),
                float(row['SE1']), float(row['SE2']), float(row['SE3']),
                float(row['SK1']), float(row['SK2']), float(row['SK3'])))
            tprint('[%04d] sv_sp=[%.10e %.10e %.10e]  sk_sp=[%.10e %.10e %.10e]'
                   % (case_id, sv_sp[0], sv_sp[1], sv_sp[2],
                      sk_sp[0], sk_sp[1], sk_sp[2]))
            tprint('[%04d] sv_mat=[%.10e %.10e %.10e]  sk_mat=[%.10e %.10e %.10e]'
                   % (case_id, sv_mat[0], sv_mat[1], sv_mat[2],
                      sk_mat[0], sk_mat[1], sk_mat[2]))

            sf1 = float(row.get('SF1', 0.0)); sf2 = float(row.get('SF2', 0.0))
            sf3 = float(row.get('SF3', 0.0)); sm1 = float(row.get('SM1', 0.0))
            sm2 = float(row.get('SM2', 0.0)); sm3 = float(row.get('SM3', 0.0))
            n_t_sp = np.array([sf3, sf2, sf1])
            m_t_sp = np.array([sm1, sm2, sm3])

            base = tag
            Phi, n_out, m_out = solve_and_postprocess(
                n_t_sp, m_t_sp, sv_mat, sk_mat, case_dir, base, R_local, R_MAP)
            abq_n = R_local @ n_t_sp
            abq_m = R_local @ m_t_sp

            tprint('[%04d] SF/SM raw (local) n=%s m=%s' % (case_id, n_t_sp, m_t_sp))
            tprint('[%04d] SF/SM rotated (global) n=%s m=%s' % (case_id, abq_n, abq_m))

            pd.DataFrame([{
                'Z_mid':config.Z_MID,
                'Phi':Phi, 'Phi_norm': Phi/config.ENERGY_SCALE,
                'n1':n_out[0],'n2':n_out[1],'n3':n_out[2],
                'n1_norm':n_out[0]/config.FORCE_SCALE,
                'n2_norm':n_out[1]/config.FORCE_SCALE,
                'n3_norm':n_out[2]/config.FORCE_SCALE,
                'm1':m_out[0],'m2':m_out[1],'m3':m_out[2],
                'm1_norm':m_out[0]/config.MOMENT_SCALE,
                'm2_norm':m_out[1]/config.MOMENT_SCALE,
                'm3_norm':m_out[2]/config.MOMENT_SCALE,
            }]).to_csv(os.path.join(case_dir, base+'_resultants.csv'),
                       index=False, float_format='%.8e')
        except Exception as ex:
            msg = str(ex)
            if 'did not converge' in msg or 'DIVERGED' in msg or 'nonlinear' in msg.lower():
                tprint('[%04d] FEniCS did not converge — discarding case' % case_id)
                return DISCARD_FENICS_DIVERGED
            tprint('[%04d] FEniCS/post-extract error (%s) — discarding case'
                   % (case_id, ex))
            return DISCARD_FENICS_ERROR
    else:
        Phi = np.nan; n_out = m_out = np.zeros(3)
        abq_n = abq_m = np.zeros(3)

    # 3. comparison (in-process — pyvista available in WSL)
    solid_csv = os.path.join(case_dir, 'solid_mid.csv')
    vtus = glob.glob(os.path.join(case_dir, '*_ref.vtu'))
    if not dry and vtus and os.path.isfile(solid_csv):
        try:
            _compare_inprocess(case_dir, solid_csv, vtus[0], n_out, m_out,
                               abq_n, abq_m, Phi, R=R_local)
        except Exception as ex:
            tprint('[%04d] WARNING: compare failed: %s' % (case_id, ex))

    # 4. collect scalars
    eps_analytical = strain_screen(vx, vy, mx)
    result = {'case':case_id,'Vx':vx,'Vy':vy,'Mx':mx,
              'eps_root':eps_analytical,
              'eps_root_abaqus': eps_root_abaqus,
              'z_root_used': z_root_used,
              'root_strain_under_threshold': (eps_root_abaqus <= config.STRAIN_THRESHOLD)
                                   if eps_root_abaqus is not None else None,
              'bending_dominance_mid': bending_dominance(vx, vy, mx, config.Z_MID),
              'bending_dominance_root': bending_dominance(vx, vy, mx, 0.0),
              'bending_dominance_mid_abaqus': bd_mid_abaqus,
              'eps_mid_analytical': strain_screen_at_station(vx, vy, mx, config.Z_MID),
              'UR1': ur_vals[0], 'UR2': ur_vals[1], 'UR3': ur_vals[2],
              'Phi_FEniCS': Phi if not dry else np.nan,
              'm_mid_nd': mid_moment_nd(vx, vy, mx)}
    result.update(strain_agree)
    result.update(strain3)

    err_path = os.path.join(case_dir, 'errors.csv')
    if os.path.isfile(err_path) and not dry:
        edf = pd.read_csv(err_path).set_index('component')
        for comp in ['Mises','inv_1','J3','S33','S13','S23','N_resultant','M_resultant','Energy']:
            if comp in edf.index:
                result['relL2_'+comp] = float(edf.loc[comp,'relL2'])
                result['L2_'+comp]    = float(edf.loc[comp,'L2'])
                _sc = ('L2_over_scale' if 'L2_over_scale' in edf.columns
                       else ('L2_normalized' if 'L2_normalized' in edf.columns
                             else None))
                if _sc is not None:
                    result['L2_over_scale_'+comp] = float(edf.loc[comp,_sc])
                if 'relL2_resultant' in edf.columns:
                    result['relL2_resultant_'+comp] = float(edf.loc[comp,'relL2_resultant'])
                if 'relL2_Szz' in edf.columns:
                    result['relL2_Szz_'+comp] = float(edf.loc[comp,'relL2_Szz'])
    elif dry:
        rng = np.random.default_rng(case_id)
        for comp in ['Mises','inv_1','J3','S33','S13','S23','N_resultant','M_resultant','Energy']:
            result['relL2_'+comp] = float(rng.uniform(0.01,0.15))
            result['L2_'+comp]    = float(rng.uniform(1e-8,1e-5))

    energy_path = os.path.join(case_dir, 'solid_energy.csv')
    if os.path.isfile(energy_path) and not dry:
        result['ELSE_per_length'] = float(
            pd.read_csv(energy_path)['U_per_length_ELSE'].iloc[0])

    if not dry:
        result.update(read_case_diagnostics(case_dir))
    return result


def main():
    t_main_start = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument('--cases',   type=int, default=config.N_TARGET)
    ap.add_argument('--seed',    type=int, default=42)
    ap.add_argument('--loads-csv', type=str, default=None,
                    help='Replay Vx/Vy/Mx from a file written by '
                         'loads_io.export_accepted_loads() instead of '
                         'drawing a fresh random sample.')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--workers', type=int, default=config.MAX_PARALLEL,
                    help='parallel cases (default %(default)s)')
    ap.add_argument('--redo-fenics', action='store_true')
    ap.add_argument('--redo-compare', action='store_true')
    ap.add_argument('--min-bending-dominance', type=float, default=None)
    ap.add_argument('--min-bending-dominance-post', type=float,
                    default=config.MIN_BENDING_DOMINANCE_POST)
    ap.add_argument('--max-post-strain', type=float,
                    default=MAX_POST_STRAIN_DEFAULT)
    ap.add_argument('--max-1d-residual', type=float, default=0.01,
                    help='discard a case if the 1D Abaqus residual force/moment '
                         'exceeds this fraction of the applied load (default '
                         '%(default)s; pass a negative value to disable)')
    ap.add_argument('--max-shear-bending', type=float, default=config.BENDING_GATE_MAX,
                    help='pre-run bending gate (gate table): skip draws whose '
                         'transverse-shear/bending stress ratio at Z_MID exceeds '
                         'this (default %(default)s; negative disables)')
    ap.add_argument('--max-hourglass', type=float, default=config.HOURGLASS_MAX,
                    help='discard a case if 3D ALLAE/ALLIE exceeds this '
                         '(default %(default)s; negative disables)')
    ap.add_argument('--min-mid-moment-nd', type=float, default=config.MID_FLOOR_ND,
                    help='pre-run: skip draws whose mid-span bending moment '
                         '|m_x|+|m_y| (in M0 units) is below this, i.e. moment '
                         'zero crossings near Z_MID. Size-independent. Default '
                         '%(default)s (gate table); negative disables.')
    ap.add_argument('--element-3d', type=str, default=None,
                    help='3D solid element type, e.g. C3D8I, C3D8, C3D8R '
                         '(only 8-node bricks; default: whatever the 3D '
                         'template has, normally C3D8R). Same as the '
                         'BATCH_ELEMENT_3D environment variable.')
    ap.add_argument('--out-dir', type=str, required=True,
                    help='folder for this batch\'s case_* subfolders (required)')
    args = ap.parse_args()
    N_DRAW = int(np.ceil(args.cases / 0.80))

    cases_dir = os.path.abspath(args.out_dir)
    if args.out_dir:
        tprint('Output directory overridden: %s' % cases_dir)

    os.makedirs(cases_dir, exist_ok=True)

    if args.element_3d:
        config.ELEMENT_3D = args.element_3d
    if config.ELEMENT_3D:
        from abaqus_io import set_solid_element
        set_solid_element('*Element, type=C3D8R\n', config.ELEMENT_3D)  # validates the name
        tprint('3D element type overridden: %s (template: %s)'
               % (config.ELEMENT_3D.upper(), os.path.basename(config.INP_3D)))
    try:
        with open(os.path.join(cases_dir, 'run_settings.txt'), 'a') as _rs:
            _rs.write('%s  section_side=%s  E=%s  inp_3d=%s  element_3d=%s  loads_csv=%s\n'
                      % (time.strftime('%Y-%m-%d %H:%M:%S'), config.SECTION_SIDE,
                         config.E, os.path.basename(config.INP_3D),
                         (config.ELEMENT_3D.upper() if config.ELEMENT_3D else 'template (C3D8R)'),
                         args.loads_csv))
    except OSError:
        pass

    range_safety_check(verbose=True)
    if check_section_matches_inp() is False and not args.dry_run:
        raise SystemExit('config.SECTION_SIDE does not match the beam section '
                         'in %s -- fix one of them before running (README, '
                         'section 7).' % config.INP_1D)
    # CHANGE 7: log the station list once at startup for auditability
    tprint('Z_STATIONS (step=10): %s' % _z_stations_arg(config.L, step=10.0))

    skip_fenics_init = args.redo_compare and not args.redo_fenics
    if not args.dry_run and not skip_fenics_init:
        init_fenics(config.MESH_XMF, cases_dir)
    elif args.redo_compare:
        tprint('--redo-compare: skipping FEniCS mesh/problem build.')

    if args.loads_csv:
        draws = load_fixed_loads(args.loads_csv)
    else:
        rng = np.random.default_rng(args.seed)
        draws = pd.DataFrame({
            'Vx': rng.uniform(*config.VX_RANGE, N_DRAW),
            'Vy': rng.uniform(*config.VY_RANGE, N_DRAW),
            'Mx': rng.uniform(*config.MX_RANGE, N_DRAW),
        })

    accepted_draws = []
    discarded = 0
    discarded_bending = 0
    discarded_midfloor = 0
    discarded_shearbend = 0
    for _, row in draws.iterrows():
        # With --loads-csv, screen EVERY row: the rows beyond --cases are kept
        # as spares and used (in file order) to replace cases rejected after
        # the run, before any random top-up draw (added 7 Oct 2026).
        if len(accepted_draws) >= args.cases and not args.loads_csv: break
        eps = strain_screen(row['Vx'], row['Vy'], row['Mx'])   # logged as eps_root
        # Strain gate per the gate table (7 Oct 2026): (c/L) max(|m-vy|+|vx|, |m|)
        if table_strain_gate(row['Vx'], row['Vy'], row['Mx']) > config.STRAIN_THRESHOLD:
            discarded += 1
            continue
        if (args.max_shear_bending >= 0 and
                table_bending_ratio(row['Vx'], row['Vy'], row['Mx']) > args.max_shear_bending):
            discarded_shearbend += 1
            continue
        if args.min_bending_dominance is not None:
            bd = bending_dominance(row['Vx'], row['Vy'], row['Mx'], config.Z_MID)
            if bd < args.min_bending_dominance:
                discarded_bending += 1
                continue
        if (args.min_mid_moment_nd is not None and args.min_mid_moment_nd >= 0 and
                mid_moment_nd(row['Vx'], row['Vy'], row['Mx']) < args.min_mid_moment_nd):
            discarded_midfloor += 1
            continue
        accepted_draws.append((len(accepted_draws),
                               row['Vx'], row['Vy'], row['Mx'], eps))

    tprint('Screened: %d accepted, %d discarded (strain threshold %.1e)%s'
           % (len(accepted_draws), discarded, config.STRAIN_THRESHOLD,
              ', %d discarded (bending_dominance < %.1f at Z_MID)'
              % (discarded_bending, args.min_bending_dominance)
              if args.min_bending_dominance is not None else ''))

    results      = []
    t_batch      = time.time()
    n_discard    = {r: 0 for r in DISCARD_REASONS}
    n_post_gate_failed = 0

    def _discard_str():
        parts = ['%s=%d' % (r, n) for r, n in n_discard.items() if n]
        return ', '.join(parts) if parts else 'none'
    spare_draws = []
    if args.loads_csv and len(accepted_draws) > args.cases:
        spare_draws = [d[1:] for d in accepted_draws[args.cases:]]   # (vx, vy, mx, eps)
        accepted_draws = accepted_draws[:args.cases]
        tprint('Loads file: %d spare screened rows kept for replacing rejected cases'
               % len(spare_draws))
    case_counter = len(accepted_draws)
    pending      = list(accepted_draws)

    n_hard = 0
    consecutive_hard = 0
    abort_reason = None

    while len(results) < args.cases and pending and abort_reason is None:
        free = _free_gb(cases_dir)
        if free < MIN_FREE_GB:
            abort_reason = ('only %.1f GB free on the output drive (< %.1f GB)'
                            % (free, MIN_FREE_GB))
            break
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futures = {}
            for case_id, vx, vy, mx, eps in pending:
                case_dir = os.path.join(cases_dir, 'case_%04d' % case_id)
                fut = ex.submit(run_case, case_id, vx, vy, mx,
                                case_dir, args.dry_run, args.redo_fenics,
                                args.redo_compare, args.max_post_strain,
                                args.min_bending_dominance_post,
                                (args.max_1d_residual if args.max_1d_residual >= 0 else None),
                                (args.max_hourglass if args.max_hourglass >= 0 else None))
                futures[fut] = (case_id, eps)
            pending = []

            for fut in as_completed(futures):
                case_id, eps = futures[fut]
                if fut.cancelled():      # dropped after an abort
                    continue
                try:
                    res = fut.result()
                except Exception as e:
                    tprint('[%04d] EXCEPTION: %s' % (case_id, e))
                    res = None
                    if isinstance(e, OSError) and getattr(e, 'errno', None) == 28:
                        abort_reason = 'disk full (No space left on device)'
                        for f in futures: f.cancel()   # drop queued cases
                if res is not None and not isinstance(res, str):
                    res['eps_root'] = eps
                    fail_reasons = []
                    # Root-strain re-check removed 7 Oct 2026: it only repeated
                    # small_ok (max 1D strain over ALL stations incl. Z=0, same
                    # 5e-3 limit), so it never rejected anything new.
                    # root_strain_under_threshold is still written to summary.csv.
                    bd_abq = res.get('bending_dominance_mid_abaqus')
                    bd_known = bd_abq is not None and not (
                        isinstance(bd_abq, float) and np.isnan(bd_abq))
                    if bd_known and bd_abq < args.min_bending_dominance_post:
                        fail_reasons.append(
                            'bending_dominance_mid_abaqus %.3g < %.3g '
                            '(peak axial stress not dominant in the actual '
                            'solved Abaqus output)'
                            % (bd_abq, args.min_bending_dominance_post))
                    if res.get('admit_all') is False:
                        parts = []
                        if res.get('bending_ok') is False:
                            parts.append('shear/torsion stress dominant '
                                         '(bending_dominance_mid_abaqus=%.3g)'
                                         % res.get('bending_dominance_mid_abaqus',
                                                   float('nan')))
                        if res.get('small_ok') is False:
                            parts.append('strain too large (max_eps_1d_beam=%.3e '
                                         'at Z=%s over %d row(s), limit %.1e; eps_3d_ref=%s not gated)'
                                         % (res.get('max_eps_1d_beam', float('nan')),
                                            ('%.2f' % res['z_at_max_eps_1d']) if res.get('z_at_max_eps_1d') is not None else 'n/a',
                                            res.get('n_beam_rows_scanned', 0),
                                            args.max_post_strain,
                                            ('%.3e' % res['eps_3d_corner_post']) if res.get('eps_3d_corner_post') is not None else 'n/a'))
                        fail_reasons.append('admission check failed: ' + '; '.join(parts))
                    consecutive_hard = 0
                    if fail_reasons:
                        n_post_gate_failed += 1
                        tprint('[%04d] discarded (post-run gate): %s'
                               % (case_id, '; '.join(fail_reasons)))
                    else:
                        results.append(res)
                        abq_str = ('%.3e' % res['eps_root_abaqus']
                                   if res.get('eps_root_abaqus') is not None else 'n/a')
                        bd_str = ('%.3g' % bd_abq if bd_known else 'n/a')
                        tprint('[%04d] done  eps=%.3e (analytical)  eps_abq=%s (post)  '
                               'bd_abq=%s  elapsed=%.0fs'
                               % (case_id, eps, abq_str, bd_str, time.time()-t_batch))
                elif res in n_discard:
                    consecutive_hard = 0
                    n_discard[res] += 1
                    tprint('[%04d] discarded (%s)'
                           % (case_id, DISCARD_LABELS.get(res, res)))
                else:
                    n_hard += 1
                    consecutive_hard += 1
                    tprint('[%04d] FAILED (hard error)' % case_id)
                    if (consecutive_hard >= MAX_CONSECUTIVE_HARD
                            and abort_reason is None):
                        abort_reason = ('%d hard errors in a row -- last: '
                                        'see [%04d] above' % (consecutive_hard, case_id))
                        for f in futures: f.cancel()

        if abort_reason is not None:
            break
        shortfall = args.cases - len(results)
        if shortfall > 0:
            tprint('Topping up: need %d more cases (discarded: %s; '
                   'post_gate_failed=%d so far)'
                   % (shortfall, _discard_str(), n_post_gate_failed))
            # spares from the loads file first, in file order
            while spare_draws and len(pending) < shortfall:
                svx, svy, smx, seps = spare_draws.pop(0)
                pending.append((case_counter, svx, svy, smx, seps))
                case_counter += 1
            if spare_draws or len(pending) >= shortfall:
                extra = pd.DataFrame({'Vx': [], 'Vy': [], 'Mx': []})
            else:
                tprint('Loads file has no spares left -- topping up with random draws')
                extra_draws = int(np.ceil(shortfall / 0.76)) + 10
                rng2 = np.random.default_rng(args.seed + case_counter)
                extra = pd.DataFrame({
                    'Vx': rng2.uniform(*config.VX_RANGE, extra_draws),
                    'Vy': rng2.uniform(*config.VY_RANGE, extra_draws),
                    'Mx': rng2.uniform(*config.MX_RANGE, extra_draws),
                })
            for _, row in extra.iterrows():
                if len(pending) >= shortfall: break
                eps = strain_screen(row['Vx'], row['Vy'], row['Mx'])
                if table_strain_gate(row['Vx'], row['Vy'], row['Mx']) > config.STRAIN_THRESHOLD:
                    continue
                if (args.max_shear_bending >= 0 and
                        table_bending_ratio(row['Vx'], row['Vy'], row['Mx']) > args.max_shear_bending):
                    continue
                if args.min_bending_dominance is not None:
                    bd = bending_dominance(row['Vx'], row['Vy'], row['Mx'], config.Z_MID)
                    if bd < args.min_bending_dominance:
                        continue
                if (args.min_mid_moment_nd is not None and args.min_mid_moment_nd >= 0 and
                        mid_moment_nd(row['Vx'], row['Vy'], row['Mx']) < args.min_mid_moment_nd):
                    continue
                pending.append((case_counter,
                                row['Vx'], row['Vy'], row['Mx'], eps))
                case_counter += 1

    tprint('\nBatch complete: %d/%d target  post_gate_failed=%d  '
           'elapsed=%.0fs'
           % (len(results), args.cases, n_post_gate_failed,
              time.time()-t_batch))
    if abort_reason is not None:
        tprint('*** BATCH ABORTED: %s. Fix it and rerun the same command; '
               'finished cases are reused. ***' % abort_reason)
    tprint('Hard errors (exceptions): %d' % n_hard)
    tprint('Discarded by reason:')
    for r in DISCARD_REASONS:
        tprint('  %-40s %d' % (DISCARD_LABELS[r], n_discard[r]))
    tprint('  %-40s %d' % ('pre-run: mid-span moment below floor', discarded_midfloor))
    tprint('  %-40s %d' % ('pre-run: shear/bending gate', discarded_shearbend))
    # Keep a record of what every gate dropped (added 7 Oct 2026), so a
    # biased sample (e.g. losing the large-rotation cases) is visible.
    try:
        rows_dc = [{'reason': r, 'label': DISCARD_LABELS[r], 'count': n_discard[r]}
                   for r in DISCARD_REASONS]
        rows_dc += [{'reason': 'pre_strain_screen', 'label': 'pre-run strain screen', 'count': discarded},
                    {'reason': 'pre_bending_dominance', 'label': 'pre-run bending dominance', 'count': discarded_bending},
                    {'reason': 'pre_shear_bending_gate', 'label': 'pre-run shear/bending gate (table)', 'count': discarded_shearbend},
                    {'reason': 'pre_mid_moment_floor', 'label': 'pre-run mid-span moment floor', 'count': discarded_midfloor},
                    {'reason': 'post_run_gate', 'label': 'post-run gates in main()', 'count': n_post_gate_failed},
                    {'reason': 'hard_error', 'label': 'exceptions', 'count': n_hard},
                    {'reason': 'accepted', 'label': 'accepted', 'count': len(results)}]
        pd.DataFrame(rows_dc).to_csv(os.path.join(cases_dir, 'discard_counts.csv'), index=False)
    except Exception as ex:
        tprint('WARNING: could not write discard_counts.csv: %s' % ex)

    if not results:
        tprint('No results.'); return

    df = pd.DataFrame(results).sort_values('case')
    df.to_csv(os.path.join(cases_dir,'summary.csv'),
              index=False, float_format='%.8e')
    tprint('Saved summary.csv')
    export_accepted_loads(os.path.join(cases_dir,'summary.csv'),
                           os.path.join(cases_dir,'accepted_loads.csv'))

    if {'eps_mid_analytical', 'eps_1d_mid_post'}.issubset(df.columns):
        cmp_cols = ['case', 'eps_mid_analytical', 'eps_1d_mid_post']
        if 'eps_3d_corner_post' in df.columns:
            cmp_cols.append('eps_3d_corner_post')
        cmp_df = df[cmp_cols].dropna(subset=['eps_mid_analytical', 'eps_1d_mid_post'])
        if len(cmp_df):
            err_1d = (cmp_df['eps_1d_mid_post'] - cmp_df['eps_mid_analytical']) / cmp_df['eps_mid_analytical']
            n_nonconservative_1d = int((err_1d > 0).sum())
            tprint('\nMohr\'s-circle analytical screen vs real Abaqus, AT THE SAME '
                   'STATION (Z_MID) — %d cases:' % len(cmp_df))
            tprint('  1D (beam): median=%.1f%%  max=%.1f%%  min=%.1f%%' %
                   (100*err_1d.median(), 100*err_1d.max(), 100*err_1d.min()))
            if n_nonconservative_1d:
                tprint('  WARNING: %d/%d cases have eps_1d_mid_post > eps_mid_analytical'
                       % (n_nonconservative_1d, len(cmp_df)))
            else:
                tprint('  -> Analytical screen was conservative for every case.')
            if 'eps_3d_corner_post' in cmp_df.columns:
                # FIX: compare against the EXTREME-FIBER corner value, not
                # the centerline -- both eps_mid_analytical and this are
                # now the same physical point (see chat: the old
                # centerline comparison always showed ~99% "error"
                # because it was near the neutral axis, not a real
                # accuracy signal).
                cmp_df3 = cmp_df.dropna(subset=['eps_3d_corner_post'])
                if len(cmp_df3):
                    err_3d = (cmp_df3['eps_3d_corner_post'] - cmp_df3['eps_mid_analytical']) / cmp_df3['eps_mid_analytical']
                    n_nonconservative_3d = int((err_3d > 0).sum())
                    tprint('  3D (solid, EXTREME FIBER): median=%.1f%%  max=%.1f%%  min=%.1f%%'
                           % (100*err_3d.median(), 100*err_3d.max(), 100*err_3d.min()))
                    if n_nonconservative_3d:
                        tprint('  WARNING: %d/%d cases have eps_3d_corner_post > eps_mid_analytical'
                               % (n_nonconservative_3d, len(cmp_df3)))

    ur_cols = [c for c in ['case', 'eps_1d_mid_post', 'UR1', 'UR2', 'UR3']
               if c in df.columns]
    if len(ur_cols) > 1:
        ur_path = os.path.join(cases_dir, 'accepted_1d_strain_ur.csv')
        df[ur_cols].to_csv(ur_path, index=False, float_format='%.8e')
        tprint('Saved %s (%d accepted cases)' % (ur_path, len(df)))

    rel_cols = [c for c in df.columns if c.startswith('relL2_')]
    if rel_cols:
        tprint('\nMedian rel-L2 errors [%]:')
        tprint((df[rel_cols].median()*100).round(3).to_string())

    if 'bending_dominance_mid' in df.columns:
        bd = df['bending_dominance_mid'].replace([np.inf, -np.inf], np.nan).dropna()
        if len(bd):
            tprint('\nBending dominance at Z_MID: min=%.1f  median=%.1f  max=%.1f'
                   % (bd.min(), bd.median(), bd.max()))
            n_weak = int((bd < 10).sum())
            if n_weak:
                tprint('WARNING: %d/%d cases have bending_dominance_mid < 10'
                       % (n_weak, len(bd)))

    if 'root_strain_under_threshold' in df.columns and not df['root_strain_under_threshold'].all():
        tprint('\nWARNING: %d accepted case(s) have root_strain_under_threshold == False '
               '— reloaded from a pre-gate run?'
               % int((df['root_strain_under_threshold'] == False).sum()))

    plot_histograms(df, cases_dir)
    plot_energy_scatter(df, cases_dir)
    plot_strain_screen_audit(df, cases_dir)
    tprint('\nDone. Results in: ' + cases_dir)

    total_elapsed = time.time() - t_main_start
    tprint('Total run time: %.0fs (%s)'
           % (total_elapsed, time.strftime('%H:%M:%S', time.gmtime(total_elapsed))))


if __name__ == '__main__':
    main()
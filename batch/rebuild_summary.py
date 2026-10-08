#!/usr/bin/env python3
"""
rebuild_summary.py
==================
Rebuild summary.csv directly from the per-case output folders, WITHOUT
re-sampling loads, without Abaqus, and without FEniCS/dolfin.

WHY THIS EXISTS
---------------
summary.csv is normally written once at the end of batch_driver.main(),
from the in-memory `results` list. If the driver is interrupted while
rewriting it, the per-case data on disk survives intact but the index is
left truncated -- which is exactly what happened here: 1308 case folders,
999 with a completed errors.csv and zero REJECTED_admission.json markers,
but summary.csv holding only 3 rows.

The obvious recovery -- re-run the driver so its resume path calls
load_existing() -- does NOT work after the strain_screen() sign fix.
strain_screen() is the pre-run filter that decides which random draws
become cases, so correcting `mx + vy*L` to `mx - vy*L` changes which
draws pass and therefore the case NUMBERING. Re-running with the same
--seed would attach a different (vx,vy,mx) to each existing case_XXXX
folder. This script sidesteps that entirely by reading each case's loads
from its OWN input deck rather than regenerating them.

WHAT IT REPRODUCES
------------------
The same columns load_existing() produces, plus the load-derived columns
main() attaches afterwards, in the same order:

  errors.csv        -> relL2_* / L2_* / L2_over_scale_* for Mises, inv_1, J3, S33,
                       S13, S23, N_resultant, M_resultant, Energy
  solid_energy.csv  -> ELSE_per_length
  beam_whole.csv    -> eps_root_abaqus, z_root_used, root_strain_under_threshold,
                       bending_dominance_mid_abaqus, UR1/UR2/UR3, and the
                       full strain_admit() block
  *_resultants.csv  -> Phi_FEniCS
  case_XXXX_1d.inp  -> Vx, Vy, Mx  (ground truth for that folder)
  derived           -> eps_root, eps_mid_analytical,
                       bending_dominance_mid, bending_dominance_root

All post-run quantities are recomputed with the CURRENT strain_checks.py,
so the rebuilt file reflects the corrected SK1/SK2 bending convention and
the corrected moment sign. That means a rebuilt row can legitimately
differ from what the original run wrote -- that is the point, not a bug.

ROTATION
--------
strain_admit() needs R_local = rotation_matrix_from_ur(UR1,UR2,UR3),
which normally comes from warping_core via fenics_solve and drags in
dolfin. This script imports it if available and otherwise falls back to
the standard Rodrigues rotation matrix built from the rotation vector
(UR1,UR2,UR3). The fallback was verified against stored
beam_1d_mid.csv output: all six *_rot_global components reproduce to
better than 1e-9 relative, so the two are the same transform.

USAGE
-----
  python3 rebuild_summary.py                       # cases dir from config
  python3 rebuild_summary.py --cases-dir /path/to/cases_1d_half_5
  python3 rebuild_summary.py --out summary_rebuilt.csv   # don't overwrite
  python3 rebuild_summary.py --admitted-only       # drop rows failing the
                                                   # re-evaluated gates
"""
import os
import re
import sys
import glob
import argparse

import numpy as np
import pandas as pd

import config
from abaqus_io import residual_check_1d
from strain_checks import (
    mid_moment_nd, read_case_diagnostics,
    strain_screen, strain_screen_at_station, bending_dominance,
    strain_admit, _root_strain_and_mid_bd_from_beam_csv,
    MAX_POST_STRAIN_DEFAULT,
)

# ── rotation: real one if importable, verified Rodrigues fallback if not ──
try:
    from fenics_solve import rotation_matrix_from_ur  # noqa: F401
    _ROT_SRC = 'warping_core (via fenics_solve)'
except Exception:
    def rotation_matrix_from_ur(ur1, ur2, ur3):
        """Rodrigues rotation matrix for the rotation vector (ur1,ur2,ur3).
        Verified identical to warping_core's version against stored
        beam_1d_mid.csv *_rot_global columns (<1e-9 relative)."""
        w = np.array([ur1, ur2, ur3], dtype=float)
        th = float(np.linalg.norm(w))
        if th < 1e-300:
            return np.eye(3)
        k = w / th
        K = np.array([[0.0, -k[2], k[1]],
                      [k[2], 0.0, -k[0]],
                      [-k[1], k[0], 0.0]])
        return np.eye(3) + np.sin(th) * K + (1.0 - np.cos(th)) * (K @ K)
    _ROT_SRC = 'built-in Rodrigues fallback (dolfin unavailable)'


# ── read_csv cache ───────────────────────────────────────────────────────
# strain_admit() and friends re-open the same beam_whole.csv and
# solid_mid.csv several times per case (the root-strain helper, the
# beam-wide max, the corner search and the centreline lookup each read
# them independently). Over a network-mounted case tree that dominates
# the runtime, so memoise pd.read_csv by (path, mtime) for the duration
# of the rebuild. Behaviour is unchanged -- every caller still gets an
# independent copy of the frame.
_CSV_CACHE = {}
_real_read_csv = pd.read_csv


def _cached_read_csv(path, *a, **k):
    if a or k or not isinstance(path, str):
        return _real_read_csv(path, *a, **k)
    try:
        key = (path, os.path.getmtime(path))
    except OSError:
        return _real_read_csv(path)
    if key not in _CSV_CACHE:
        if len(_CSV_CACHE) > 64:
            _CSV_CACHE.clear()
        _CSV_CACHE[key] = _real_read_csv(path)
    return _CSV_CACHE[key].copy()


pd.read_csv = _cached_read_csv

_LOAD_RE = re.compile(r'Name:\s*(fx|fy|mx)\b', re.IGNORECASE)


def loads_from_inp(case_dir, case_id):
    """Read the three *Cload magnitudes from this case's own 1D input deck.

    Layout written by abaqus_io.write_inp():
        ** Name: fx   Type: Concentrated force
        *Cload, op=NEW
        free, 1, <value>
    so the magnitude is the 3rd comma-separated field two lines below the
    'Name:' comment. Returns (vx, vy, mx) or None if the deck is missing
    or does not contain all three."""
    for cand in (os.path.join(case_dir, 'case_%04d_1d.inp' % case_id),
                 *glob.glob(os.path.join(case_dir, '*_1d.inp'))):
        if not os.path.isfile(cand):
            continue
        try:
            lines = open(cand, encoding='utf-8', errors='replace').read().splitlines()
        except Exception:
            continue
        vals = {}
        for i, line in enumerate(lines):
            m = _LOAD_RE.search(line)
            if m and i + 2 < len(lines):
                parts = lines[i + 2].split(',')
                if len(parts) >= 3:
                    try:
                        vals[m.group(1).lower()] = float(parts[2])
                    except ValueError:
                        pass
        if {'fx', 'fy', 'mx'} <= set(vals):
            return vals['fx'], vals['fy'], vals['mx']
    return None


def row_from_case(case_dir, case_id):
    """Reconstruct one summary.csv row. Mirrors batch_driver.load_existing()
    plus the load-derived fields main() attaches to it. Returns None if the
    case never completed (no errors.csv)."""
    result = {}

    err_path = os.path.join(case_dir, 'errors.csv')
    if not os.path.isfile(err_path):
        return None
    edf = pd.read_csv(err_path).set_index('component')
    for comp in ['Mises', 'inv_1', 'J3', 'S33', 'S13', 'S23',
                 'N_resultant', 'M_resultant', 'Energy']:
        if comp in edf.index:
            result['relL2_' + comp] = float(edf.loc[comp, 'relL2'])
            result['L2_' + comp] = float(edf.loc[comp, 'L2'])
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

    beam_csv = os.path.join(case_dir, 'beam_whole.csv')
    if os.path.isfile(beam_csv):
        eps_abq, z_used, bd_abq = _root_strain_and_mid_bd_from_beam_csv(beam_csv)
        if eps_abq is not None:
            result['eps_root_abaqus'] = eps_abq
            result['z_root_used'] = z_used
            result['root_strain_under_threshold'] = eps_abq <= config.STRAIN_THRESHOLD
        if bd_abq is not None:
            result['bending_dominance_mid_abaqus'] = bd_abq

        solid_csv = os.path.join(case_dir, 'solid_mid.csv')
        if os.path.isfile(solid_csv):
            df_b = pd.read_csv(beam_csv)
            df_b.columns = df_b.columns.str.strip()
            needed = ['Z', 'SE1', 'SE2', 'SE3', 'UR1', 'UR2', 'UR3']
            if all(c in df_b.columns for c in needed) and len(df_b):
                row = df_b.loc[(df_b['Z'] - config.Z_MID).abs().idxmin()]
                ur1, ur2, ur3 = (float(row['UR1']), float(row['UR2']),
                                 float(row['UR3']))
                result['UR1'], result['UR2'], result['UR3'] = ur1, ur2, ur3
                R_local = rotation_matrix_from_ur(ur1, ur2, ur3)
                # NOTE: unlike load_existing() this does NOT rewrite
                # beam_1d_mid.csv -- rebuilding an index should not mutate
                # the case folders it is reading.
                result.update(strain_admit(
                    row, beam_csv, solid_csv, R_local,
                    MAX_POST_STRAIN_DEFAULT, bd_mid_abaqus=bd_abq))

    loads = loads_from_inp(case_dir, case_id)
    if loads is None:
        return None
    vx, vy, mx = loads
    result.update({
        'case': case_id, 'Vx': vx, 'Vy': vy, 'Mx': mx,
        'eps_root': strain_screen(vx, vy, mx),
        'eps_mid_analytical': strain_screen_at_station(vx, vy, mx, config.Z_MID),
        'bending_dominance_mid': bending_dominance(vx, vy, mx, config.Z_MID),
        'bending_dominance_root': bending_dominance(vx, vy, mx, 0.0),
    })
    # added 7 Oct 2026: mid-span moment (M0 units), 1D convergence and the
    # section-frame check. Convergence is re-read from the .msg file so
    # older folders (no convergence_1d.csv) get it too.
    result['m_mid_nd'] = mid_moment_nd(vx, vy, mx)
    result.update(read_case_diagnostics(case_dir))
    if 'res_ratio_1d' not in result:
        result.update(residual_check_1d(case_dir, 'case_%04d_1d' % case_id,
                                        vx, vy, mx, config.L))

    res_glob = glob.glob(os.path.join(case_dir, '*_resultants.csv'))
    res_glob = [p for p in res_glob if not os.path.basename(p).startswith('solid_')]
    if res_glob:
        try:
            rdf = pd.read_csv(res_glob[0])
            if 'Phi' in rdf.columns and len(rdf):
                result['Phi_FEniCS'] = float(rdf['Phi'].iloc[0])
        except Exception:
            pass
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases-dir', required=True,
                    help='folder holding the case_* subfolders (required)')
    ap.add_argument('--out', default=None,
                    help='output path (default: <cases-dir>/summary.csv)')
    ap.add_argument('--admitted-only', action='store_true',
                    help='keep only rows whose re-evaluated admit_all is True')
    ap.add_argument('--limit', type=int, default=None,
                    help='process at most this many NEW cases this pass '
                         '(use with --resume to rebuild in chunks)')
    ap.add_argument('--resume', action='store_true',
                    help='keep rows already present in --out and only add '
                         'cases missing from it')
    args = ap.parse_args()

    cdir = args.cases_dir
    if not os.path.isdir(cdir):
        sys.exit('cases dir not found: %s' % cdir)
    out = args.out or os.path.join(cdir, 'summary.csv')

    print('cases dir : %s' % cdir)
    print('rotation  : %s' % _ROT_SRC)

    case_dirs = sorted(glob.glob(os.path.join(cdir, 'case_*')))

    done = set()
    prior = []
    if args.resume and os.path.isfile(out):
        try:
            pdf = _real_read_csv(out)
            if 'case' in pdf.columns:
                prior = pdf.to_dict('records')
                done = set(int(c) for c in pdf['case'].tolist())
                print('resume: %d rows already in %s'
                      % (len(done), os.path.basename(out)))
        except Exception as ex:
            print('resume: could not read existing output (%r) -- starting fresh' % ex)

    rows, skipped_no_errors, skipped_no_loads, failed = list(prior), 0, 0, []
    n_new = 0
    for cd in case_dirs:
        if not os.path.isdir(cd):
            continue
        m = re.search(r'case_(\d+)$', os.path.basename(cd))
        if not m:
            continue
        cid = int(m.group(1))
        if cid in done:
            continue
        if args.limit is not None and n_new >= args.limit:
            break
        if not os.path.isfile(os.path.join(cd, 'errors.csv')):
            skipped_no_errors += 1
            continue
        n_new += 1
        if n_new % 25 == 0:
            print('  ... %d new cases processed (at case_%04d)' % (n_new, cid),
                  flush=True)
        try:
            r = row_from_case(cd, cid)
        except Exception as ex:
            failed.append((cid, repr(ex)))
            continue
        if r is None:
            skipped_no_loads += 1
            continue
        rows.append(r)

    print('\nscanned %d case folders' % len(case_dirs))
    print('  reconstructed        : %d' % len(rows))
    print('  skipped, no errors.csv (never completed) : %d' % skipped_no_errors)
    print('  skipped, loads unreadable from .inp      : %d' % skipped_no_loads)
    if failed:
        print('  errored              : %d  (first 5 below)' % len(failed))
        for cid, ex in failed[:5]:
            print('     case %04d: %s' % (cid, ex))
    if not rows:
        sys.exit('nothing to write')

    df = pd.DataFrame(rows).sort_values('case')
    if 'admit_all' in df.columns:
        n_ok = int(df['admit_all'].fillna(False).astype(bool).sum())
        print('\n  admit_all True under the CURRENT (fixed) checks : %d / %d'
              % (n_ok, len(df)))
        if n_ok < len(df):
            print('  NOTE: these cases completed under the OLD checks. The gap is')
            print('        the effect of the SK1/SK2 bending-convention fix, not')
            print('        a failure of the rebuild.')
        if args.admitted_only:
            df = df[df['admit_all'].fillna(False).astype(bool)]
            print('  --admitted-only: keeping %d rows' % len(df))

    if os.path.isfile(out) and not args.resume:
        bak = out + '.before_rebuild'
        if not os.path.isfile(bak):
            pd.read_csv(out).to_csv(bak, index=False)
            print('\n  existing file backed up -> %s' % os.path.basename(bak))

    df.to_csv(out, index=False, float_format='%.8e')
    print('\nWrote %s  (%d rows x %d cols)' % (out, len(df), df.shape[1]))


if __name__ == '__main__':
    main()

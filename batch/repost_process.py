#!/usr/bin/env python3
"""
repost_process.py
=================
Regenerate the per-case derived CSVs using the CORRECTED strain_checks.py,
without re-running Abaqus or FEniCS. Companion to rebuild_summary.py.

WHAT IT REWRITES, per case folder
---------------------------------
  force_moment_1d3d_error.csv
      Was written with m_sp = [SM3, SM2, SM1], which put the TORQUE in the
      Mx_1d slot and a BENDING moment in the Mz_1d slot, so the file
      reported errors of -96.8% and +2974% where the true agreement is
      about 1%. Regenerated with m_sp = [SM1, SM2, SM3].

  beam_1d_mid.csv
      Column names bend1/bend2/bend3_rot_global asserted that all three
      rotated curvatures are bending. SK3 is the TWIST. Renamed to
      SK1_bend_rot_global / SK2_bend_rot_global / SK3_twist_rot_global.
      eps_1d_mid_post_rotated is also recomputed, since the underlying
      principal-strain formula changed (SK2+SK3 -> SK1+SK2).

WHAT IT DOES NOT TOUCH
----------------------
  Nothing that came out of Abaqus or FEniCS: no .odb, no beam_whole.csv,
  no solid_mid.csv, no solid_resultants.csv, no errors.csv. Every file
  written here is derived arithmetic that can be regenerated at will, so
  no backup is taken.

  It also cannot fix the load-space coverage gap. strain_screen() decides
  which draws ever became cases, so its sign fix only takes effect on a
  fresh run. Post-processing corrects how existing cases are MEASURED,
  not which ones exist.

IDEMPOTENT / RESUMABLE
----------------------
  A case whose beam_1d_mid.csv already carries the new column names is
  treated as done and skipped, so the script can be run in chunks with
  --limit over a slow mount and simply re-run until it reports 0 left.

USAGE
-----
  python3 repost_process.py --cases-dir ../cases_1d_half_5
  python3 repost_process.py --cases-dir ../cases_1d_half_5 --limit 300
  python3 repost_process.py --cases-dir ../cases_1d_half_5 --dry-run
"""
import os
import re
import sys
import glob
import argparse

import numpy as np
import pandas as pd

import config
from strain_checks import (write_force_moment_1d3d_csv, write_beam_1d_row_csv,
                           principal_strain_from_beam_row_rotated)

try:
    from fenics_solve import rotation_matrix_from_ur
    _ROT = 'warping_core (via fenics_solve)'
except Exception:
    def rotation_matrix_from_ur(ur1, ur2, ur3):
        """Rodrigues matrix from the rotation vector. Verified identical to
        warping_core's against stored *_rot_global output (<1e-9 rel)."""
        w = np.array([ur1, ur2, ur3], dtype=float)
        th = float(np.linalg.norm(w))
        if th < 1e-300:
            return np.eye(3)
        k = w / th
        K = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
        return np.eye(3) + np.sin(th) * K + (1.0 - np.cos(th)) * (K @ K)
    _ROT = 'built-in Rodrigues fallback (dolfin unavailable)'

NEW_COLS = {'SK1_bend_rot_global', 'SK2_bend_rot_global', 'SK3_twist_rot_global'}


def already_done(case_dir):
    p = os.path.join(case_dir, 'beam_1d_mid.csv')
    if not os.path.isfile(p):
        return False
    try:
        with open(p) as fh:
            header = set(fh.readline().strip().split(','))
        return NEW_COLS <= header
    except Exception:
        return False


def process(case_dir, dry=False):
    """Returns (did_fm, did_mid, note)."""
    beam_csv = os.path.join(case_dir, 'beam_whole.csv')
    res_csv = os.path.join(case_dir, 'solid_resultants.csv')
    if not os.path.isfile(beam_csv):
        return False, False, 'no beam_whole.csv'

    df = pd.read_csv(beam_csv)
    df.columns = df.columns.str.strip()
    if 'Z' not in df.columns or not len(df):
        return False, False, 'beam_whole.csv unusable'

    did_fm = False
    if os.path.isfile(res_csv):
        if not dry:
            out = write_force_moment_1d3d_csv(case_dir, beam_csv, res_csv,
                                              rotation_matrix_from_ur)
            did_fm = out is not None
        else:
            did_fm = True

    did_mid = False
    need = ['SE1', 'SE2', 'SE3', 'SK1', 'SK2', 'SK3', 'UR1', 'UR2', 'UR3']
    if all(c in df.columns for c in need):
        row = df.loc[(df['Z'] - config.Z_MID).abs().idxmin()]
        R = rotation_matrix_from_ur(float(row['UR1']), float(row['UR2']),
                                    float(row['UR3']))
        if not dry:
            eps = principal_strain_from_beam_row_rotated(row, R)
            write_beam_1d_row_csv(case_dir, row, R, eps)
        did_mid = True
    return did_fm, did_mid, ''


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases-dir', required=True,
                    help='folder holding the case_* subfolders (required)')
    ap.add_argument('--limit', type=int, default=None,
                    help='process at most this many cases this pass')
    ap.add_argument('--dry-run', action='store_true',
                    help='report what would change, write nothing')
    ap.add_argument('--force', action='store_true',
                    help='reprocess cases already carrying the new columns')
    args = ap.parse_args()

    cdir = args.cases_dir
    if not os.path.isdir(cdir):
        sys.exit('cases dir not found: %s' % cdir)
    print('cases dir : %s' % cdir)
    print('rotation  : %s' % _ROT)
    if args.dry_run:
        print('DRY RUN -- nothing will be written')

    dirs = sorted(d for d in glob.glob(os.path.join(cdir, 'case_*'))
                  if os.path.isdir(d))
    todo = [d for d in dirs if args.force or not already_done(d)]
    print('\n%d case folders, %d already done, %d to process'
          % (len(dirs), len(dirs) - len(todo), len(todo)))
    if args.limit:
        todo = todo[:args.limit]
        print('  --limit %d -> processing %d this pass' % (args.limit, len(todo)))
    if not todo:
        print('\nnothing left to do.')
        return

    n_fm = n_mid = 0
    skipped = []
    for i, d in enumerate(todo, 1):
        try:
            fm, mid, note = process(d, dry=args.dry_run)
        except Exception as ex:
            skipped.append((os.path.basename(d), repr(ex)))
            continue
        n_fm += bool(fm)
        n_mid += bool(mid)
        if note:
            skipped.append((os.path.basename(d), note))
        if i % 25 == 0:
            print('  ... %d/%d  (%s)' % (i, len(todo), os.path.basename(d)),
                  flush=True)

    print('\nforce_moment_1d3d_error.csv rewritten : %d' % n_fm)
    print('beam_1d_mid.csv rewritten             : %d' % n_mid)
    if skipped:
        print('skipped/errored                       : %d  (first 5)' % len(skipped))
        for nm, why in skipped[:5]:
            print('   %s: %s' % (nm, why))
    remaining = len([d for d in dirs if not already_done(d)])
    print('\nstill to do after this pass: %d' % remaining)
    if remaining:
        print('  re-run the same command until it reports 0.')


if __name__ == '__main__':
    main()

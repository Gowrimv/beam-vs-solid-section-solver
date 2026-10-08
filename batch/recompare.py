#!/usr/bin/env python3
"""
recompare.py
============
Regenerate errors.csv (and the contour/overlay plots, and force_error.csv)
for existing cases, WITHOUT Abaqus and WITHOUT re-solving FEniCS.

WHY THIS EXISTS
---------------
batch_driver.py's --redo-compare does the same recomputation, but it goes
through main(), which RE-DRAWS the loads every run and assigns case ids
sequentially over whichever draws survive strain_screen(). Since the screen
was changed (the mx +/- vy*L sign fix, and the switch to a conservative
bound), the accepted draw sequence is different -- measured with seed 42,
only 2 of 1000 case ids receive the same draw as before. --redo-compare
would therefore pair each existing .odb with a DIFFERENT (Vx,Vy,Mx) and
write those wrong loads into summary.csv, silently.

This script never re-draws. It reads each case's inputs from that case's
own folder, so the pairing is correct by construction -- the same reason
rebuild_summary.py and repost_process.py read loads from each case's .inp.

WHAT IT NEEDS (all already present in a completed case folder)
--------------------------------------------------------------
  solid_mid.csv           the 3D stress/strain field at Z_MID
  <case>_ref.vtu          the FEniCS recovered field (from the prior solve)
  <case>_resultants.csv   n1..n3, m1..m3, Phi  (FEniCS side)
  beam_whole.csv          SF/SM/UR at Z_MID    (Abaqus beam side)

WHAT IT PRODUCES
----------------
  errors.csv              now including the new relL2_resultant column
  force_error.csv
  contours.png, centerline_overlay.png, strain_overlay.png

REQUIREMENTS
------------
  pyvista + matplotlib + scipy (compare.py's imports). Run it in the same
  environment the pipeline normally runs in.

USAGE
-----
  python3 recompare.py --cases-dir /path/to/cases_1d_half_5
  python3 recompare.py --cases-dir ... --limit 200      # chunked
  python3 recompare.py --cases-dir ... --dry-run        # report only
"""
import os
import re
import sys
import glob
import argparse

import numpy as np
import pandas as pd

import config


def rotation_matrix_from_ur(ur1, ur2, ur3):
    """Rodrigues matrix from the rotation vector. Verified identical to
    warping_core's version against stored *_rot_global output (<1e-9 rel).
    Imported from fenics_solve when dolfin is available."""
    w = np.array([ur1, ur2, ur3], dtype=float)
    th = float(np.linalg.norm(w))
    if th < 1e-300:
        return np.eye(3)
    k = w / th
    K = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    return np.eye(3) + np.sin(th) * K + (1.0 - np.cos(th)) * (K @ K)


try:
    from fenics_solve import rotation_matrix_from_ur as _rmu
    rotation_matrix_from_ur = _rmu
    _ROT = 'warping_core (via fenics_solve)'
except Exception:
    _ROT = 'built-in Rodrigues fallback (dolfin unavailable)'


def inputs_for(case_dir, case_id):
    """Collect everything _compare_inprocess() needs, or None if the case
    is not complete enough to recompare."""
    tag = 'case_%04d' % case_id
    solid_csv = os.path.join(case_dir, 'solid_mid.csv')
    beam_csv = os.path.join(case_dir, 'beam_whole.csv')
    res_csv = os.path.join(case_dir, tag + '_resultants.csv')
    vtus = glob.glob(os.path.join(case_dir, tag + '_ref.vtu'))
    if not (os.path.isfile(solid_csv) and os.path.isfile(beam_csv)
            and os.path.isfile(res_csv) and vtus):
        return None

    r = pd.read_csv(res_csv)
    if not len(r):
        return None
    n_out = np.array([float(r['n1'].iloc[0]), float(r['n2'].iloc[0]),
                      float(r['n3'].iloc[0])])
    m_out = np.array([float(r['m1'].iloc[0]), float(r['m2'].iloc[0]),
                      float(r['m3'].iloc[0])])
    phi = float(r['Phi'].iloc[0]) if 'Phi' in r.columns else None

    b = pd.read_csv(beam_csv)
    b.columns = b.columns.str.strip()
    need = ['Z', 'SF1', 'SF2', 'SF3', 'SM1', 'SM2', 'SM3', 'UR1', 'UR2', 'UR3']
    if not all(c in b.columns for c in need) or not len(b):
        return None
    row = b.loc[(b['Z'] - config.Z_MID).abs().idxmin()]
    R = rotation_matrix_from_ur(float(row['UR1']), float(row['UR2']),
                                float(row['UR3']))
    # same index arrangement and same R (no transpose) as batch_driver.py:
    # axial and torque last, R_local @ v for local -> global. Verified
    # against the 3D NFORC resultants: the axial component lands at 3e-10
    # where it must be zero, versus 2.4e-06 with the transpose.
    abq_n = R @ np.array([float(row['SF3']), float(row['SF2']), float(row['SF1'])])
    abq_m = R @ np.array([float(row['SM1']), float(row['SM2']), float(row['SM3'])])
    return dict(solid_csv=solid_csv, vtu=vtus[0], n_out=n_out, m_out=m_out,
                abq_n=abq_n, abq_m=abq_m, phi=phi)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases-dir', required=True,
                    help='folder holding the case_* subfolders (required)')
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--force', action='store_true',
                    help='recompare even cases whose errors.csv already has '
                         'the relL2_resultant column')
    args = ap.parse_args()

    if not os.path.isdir(args.cases_dir):
        sys.exit('cases dir not found: %s' % args.cases_dir)
    print('cases dir : %s' % args.cases_dir)
    print('rotation  : %s' % _ROT)

    if not args.dry_run:
        try:
            from compare import _compare_inprocess
        except Exception as ex:
            sys.exit('cannot import compare.py (%r).\nrecompare needs pyvista, '
                     'matplotlib and scipy -- run this in the environment the '
                     'pipeline normally uses.' % ex)

    def done(cd):
        p = os.path.join(cd, 'errors.csv')
        if not os.path.isfile(p):
            return False
        try:
            with open(p) as fh:
                return 'relL2_resultant' in fh.readline()
        except Exception:
            return False

    dirs = sorted(d for d in glob.glob(os.path.join(args.cases_dir, 'case_*'))
                  if os.path.isdir(d))
    todo = [d for d in dirs if args.force or not done(d)]
    print('\n%d case folders, %d already have relL2_resultant, %d to do'
          % (len(dirs), len(dirs) - len(todo), len(todo)))
    if args.limit:
        todo = todo[:args.limit]
        print('  --limit %d -> %d this pass' % (args.limit, len(todo)))
    if not todo:
        print('\nnothing left to do.')
        return

    ok = incomplete = 0
    failed = []
    for i, cd in enumerate(todo, 1):
        m = re.search(r'case_(\d+)$', os.path.basename(cd))
        if not m:
            continue
        cid = int(m.group(1))
        try:
            arg = inputs_for(cd, cid)
        except Exception as ex:
            failed.append((os.path.basename(cd), repr(ex)))
            continue
        if arg is None:
            incomplete += 1
            continue
        if args.dry_run:
            ok += 1
            continue
        try:
            _compare_inprocess(cd, arg['solid_csv'], arg['vtu'],
                               arg['n_out'], arg['m_out'],
                               abq_n=arg['abq_n'], abq_m=arg['abq_m'],
                               phi_fenics=arg['phi'])
            ok += 1
        except Exception as ex:
            failed.append((os.path.basename(cd), repr(ex)))
        if i % 25 == 0:
            print('  ... %d/%d (%s)' % (i, len(todo), os.path.basename(cd)),
                  flush=True)

    print('\nrecompared          : %d' % ok)
    print('skipped, incomplete : %d  (missing *_ref.vtu / *_resultants.csv)'
          % incomplete)
    if failed:
        print('failed              : %d  (first 5)' % len(failed))
        for nm, why in failed[:5]:
            print('   %s: %s' % (nm, why))
    left = len([d for d in dirs if not done(d)])
    print('\nstill to do: %d' % left)


if __name__ == '__main__':
    main()

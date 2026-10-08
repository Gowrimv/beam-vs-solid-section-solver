#!/usr/bin/env python3
"""
reextract.py
============
Re-run extract_odb.py over existing case folders, so solid_resultants.csv
and beam_whole.csv are regenerated with the DEFORMED-CONFIGURATION moment
resultants. No solving: the .odb files are read only, so no Abaqus licence
is consumed and nothing is re-analysed.

WHY
---
NFORC gives nodal forces in the CURRENT configuration, but extract_odb.py
used to pair them with n.coordinates -- the UNDEFORMED node positions --
as lever arms, and take moments about the reference origin rather than the
deformed section centroid. Under nlgeom=YES that mixes two configurations.
The error is zero at the clamped root and grows with deformation.

Validated at the free tip, where the answer is known exactly because the
lever arm vanishes (case 0070, 18.2 deg rotation):

    applied Mx            = -1.61287000e-04
    reference config      = -1.53288787e-04    off by  -4.9590 %
    deformed  config      = -1.61286963e-04    off by  -0.0000 %

    My  reference = +1.1052e-05   deformed = +7.08e-12   1,560,792x closer to 0
    Mz  reference = -3.1376e-07   deformed = -1.03e-11      30,467x closer to 0

Mx/My/Mz in solid_resultants.csv are now the deformed-configuration values.
(The reference-configuration Mx_refcfg/My_refcfg/Mz_refcfg columns were
written for a while for differencing and have since been dropped.)

NOTE the forces (N, Vx, Vy) are unchanged: they are plain sums of nodal
forces with no lever arm, so they never depended on the coordinates.

USAGE (from WSL, where config.py's /mnt/c paths resolve)
--------------------------------------------------------
  python3 reextract.py --cases-dir /mnt/c/.../cases_1d_half_5
  python3 reextract.py --cases-dir ... --limit 50      # chunked
  python3 reextract.py --cases-dir ... --dry-run       # print commands only

Idempotent: a case whose solid_resultants.csv already carries the
moment_config column (or the older Mx_refcfg column) is treated as done
and skipped, so re-running resumes.
"""
import os
import re
import sys
import glob
import time
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

import config
from abaqus_io import win, abq, run


def z_stations_arg(L, step=10.0):
    """Evenly spaced stations 0..L inclusive, always including Z_MID."""
    stations, z = [], 0.0
    while z < L - 1e-9:
        stations.append(z)
        z += step
    stations.append(L)
    if not any(abs(s - config.Z_MID) < 1e-9 for s in stations):
        stations.append(config.Z_MID)
        stations.sort()
    return ','.join('%g' % s for s in stations)


def already_done(case_dir):
    p = os.path.join(case_dir, 'solid_resultants.csv')
    if not os.path.isfile(p):
        return False
    try:
        with open(p) as fh:
            hdr = fh.readline()
        # 'moment_config' = current extract_odb.py; 'Mx_refcfg' = the
        # earlier deformed-config version (moments already correct).
        return 'moment_config' in hdr or 'Mx_refcfg' in hdr
    except Exception:
        return False


def one(case_dir, case_id, zarg, dry):
    tag = 'case_%04d' % case_id
    odb1 = os.path.join(case_dir, tag + '_1d.odb')
    odb3 = os.path.join(case_dir, tag + '_3d.odb')
    if not (os.path.isfile(odb1) and os.path.isfile(odb3)):
        return case_id, 'missing .odb'
    rc = run(abq('python', win(config.EXTRACT_SCRIPT), win(odb1), win(odb3),
                 win(case_dir), str(config.Z_MID), zarg),
             cwd=case_dir, dry=dry)
    return case_id, ('ok' if rc == 0 else 'extract_odb rc=%s' % rc)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases-dir', required=True,
                    help='folder holding the case_* subfolders (required)')
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--workers', type=int, default=3,
                    help='concurrent abaqus python processes (default 3). '
                         'These only READ odbs, so no licence is taken.')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()

    if not os.path.isdir(args.cases_dir):
        sys.exit('cases dir not found: %s' % args.cases_dir)
    zarg = z_stations_arg(config.L, step=10.0)
    print('cases dir   : %s' % args.cases_dir)
    print('Z_STATIONS  : %s' % zarg)
    print('extract     : %s' % config.EXTRACT_SCRIPT)

    dirs = sorted(d for d in glob.glob(os.path.join(args.cases_dir, 'case_*'))
                  if os.path.isdir(d))
    todo = [d for d in dirs if args.force or not already_done(d)]
    print('\n%d case folders, %d already deformed-config, %d to do'
          % (len(dirs), len(dirs) - len(todo), len(todo)))
    if args.limit:
        todo = todo[:args.limit]
        print('  --limit %d -> %d this pass' % (args.limit, len(todo)))
    if not todo:
        print('\nnothing left to do.')
        return

    t0 = time.time()
    ok = 0
    bad = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {}
        for d in todo:
            m = re.search(r'case_(\d+)$', os.path.basename(d))
            if not m:
                continue
            futs[ex.submit(one, d, int(m.group(1)), zarg, args.dry_run)] = d
        for i, f in enumerate(as_completed(futs), 1):
            try:
                cid, status = f.result()
            except Exception as e:
                bad.append((os.path.basename(futs[f]), repr(e)))
                continue
            if status == 'ok':
                ok += 1
            else:
                bad.append(('case_%04d' % cid, status))
            if i % 20 == 0:
                el = time.time() - t0
                print('  ... %d/%d  %.0fs elapsed, ~%.0fs left'
                      % (i, len(futs), el, el / i * (len(futs) - i)), flush=True)

    print('\nre-extracted : %d' % ok)
    if bad:
        print('failed       : %d  (first 5)' % len(bad))
        for nm, why in bad[:5]:
            print('   %s: %s' % (nm, why))
    left = len([d for d in dirs if not already_done(d)])
    print('\nstill to do  : %d' % left)
    if left and not args.dry_run:
        print('  re-run the same command until it reports 0.')


if __name__ == '__main__':
    main()

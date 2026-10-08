#!/usr/bin/env python3
"""
rerun_dropped.py
================
Finish cases whose Abaqus jobs COMPLETED but were dropped by the batch
driver (e.g. the old 900 s ABAQUS_TIMEOUT_3D killed the wait, so
extraction / FEniCS / comparison never ran). No Abaqus analysis is
re-run: run_case() sees both ODBs complete ("ODBs exist - skipping
Abaqus") and only does extract_odb.py -> FEniCS -> compare. The same
post-run gates as batch_driver.main() are applied, and accepted cases are
merged into summary.csv (a timestamped backup is written first).

Candidates (auto-detected unless --cases is given): case_* folders with
  - both <tag>_1d.sta and <tag>_3d.sta reporting COMPLETED SUCCESSFULLY
  - no errors.csv (never finished) and no REJECTED_admission.json
Loads are read back from each case's own <tag>_1d.inp (*Cload values).

USAGE (from WSL, same environment as batch_driver.py)
  python3 rerun_dropped.py --out-dir /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/cases_half_ar_3d_finer_123 --dry-run
  python3 rerun_dropped.py --out-dir /mnt/c/.../cases_half_ar_3d_finer_123
  python3 rerun_dropped.py --out-dir ... --cases 30,56,73
Resumable: a case that already has errors.csv is not a candidate.
"""
import os, re, glob, time, shutil, argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd

import config
from config import tprint
from abaqus_io import job_completed, RE_1D_VX, RE_1D_VY, RE_1D_MX
import batch_driver as bd
from strain_checks import strain_screen, MAX_POST_STRAIN_DEFAULT
from fenics_solve import init_fenics
from loads_io import export_accepted_loads
from plotting import plot_histograms, plot_energy_scatter, plot_strain_screen_audit


def loads_from_inp(case_dir, tag):
    with open(os.path.join(case_dir, tag + '_1d.inp')) as f:
        t = f.read()
    return tuple(float(rx.search(t).group(2)) for rx in (RE_1D_VX, RE_1D_VY, RE_1D_MX))


def candidates(out_dir):
    out = []
    for d in sorted(glob.glob(os.path.join(out_dir, 'case_*'))):
        tag = os.path.basename(d)
        if not re.fullmatch(r'case_\d{4}', tag):
            continue
        if os.path.isfile(os.path.join(d, 'errors.csv')):
            continue
        if os.path.isfile(os.path.join(d, 'REJECTED_admission.json')):
            continue
        if job_completed(d, tag + '_1d') and job_completed(d, tag + '_3d') \
                and os.path.isfile(os.path.join(d, tag + '_1d.odb')) \
                and os.path.isfile(os.path.join(d, tag + '_3d.odb')):
            out.append(int(tag[5:]))
    return out


def gate(res, min_bd_post):
    """Same post-run gates as batch_driver.main(). Returns list of failures."""
    fails = []
    if res.get('root_strain_under_threshold') is False:
        fails.append('root strain over threshold')
    b = res.get('bending_dominance_mid_abaqus')
    if b is not None and not (isinstance(b, float) and np.isnan(b)) and b < min_bd_post:
        fails.append('bending_dominance_mid_abaqus %.3g < %.3g' % (b, min_bd_post))
    if res.get('admit_all') is False:
        fails.append('admission check failed (bending_ok=%s small_ok=%s)'
                     % (res.get('bending_ok'), res.get('small_ok')))
    return fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--cases', default=None, help='comma-separated ids; default = auto-detect')
    ap.add_argument('--workers', type=int, default=config.MAX_PARALLEL)
    ap.add_argument('--min-bending-dominance-post', type=float,
                    default=config.MIN_BENDING_DOMINANCE_POST)
    ap.add_argument('--max-post-strain', type=float, default=MAX_POST_STRAIN_DEFAULT)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--max-1d-residual', type=float, default=0.01,
                    help='same 1D convergence gate as batch_driver.py; negative disables')
    args = ap.parse_args()
    out_dir = os.path.abspath(args.out_dir)

    ids = ([int(x) for x in args.cases.split(',') if x.strip()] if args.cases
           else candidates(out_dir))
    tprint('%d case(s) to finish: %s' % (len(ids), ids))
    jobs = []
    for cid in ids:
        tag = 'case_%04d' % cid
        cdir = os.path.join(out_dir, tag)
        if not (job_completed(cdir, tag + '_1d') and job_completed(cdir, tag + '_3d')):
            tprint('[%04d] SKIP: Abaqus jobs not both complete (this script never runs Abaqus)' % cid)
            continue
        vx, vy, mx = loads_from_inp(cdir, tag)
        jobs.append((cid, vx, vy, mx, cdir))
        tprint('[%04d] Vx=%.6e Vy=%.6e Mx=%.6e' % (cid, vx, vy, mx))
    if args.dry_run or not jobs:
        tprint('Dry run / nothing to do.'); return

    init_fenics(config.MESH_XMF, out_dir)
    new_rows, t0 = [], time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(bd.run_case, cid, vx, vy, mx, cdir, False, False, False,
                          args.max_post_strain, args.min_bending_dominance_post,
                          (args.max_1d_residual if args.max_1d_residual >= 0 else None)):
                (cid, vx, vy, mx) for cid, vx, vy, mx, cdir in jobs}
        for f in as_completed(futs):
            cid, vx, vy, mx = futs[f]
            try:
                res = f.result()
            except Exception as e:
                tprint('[%04d] EXCEPTION: %s' % (cid, e)); continue
            if res is None or isinstance(res, str):
                tprint('[%04d] discarded (%s)' % (cid, bd.DISCARD_LABELS.get(res, res))); continue
            res['eps_root'] = strain_screen(vx, vy, mx)
            fails = gate(res, args.min_bending_dominance_post)
            if fails:
                tprint('[%04d] discarded (post-run gate): %s' % (cid, '; '.join(fails))); continue
            new_rows.append(res)
            tprint('[%04d] done  elapsed=%.0fs' % (cid, time.time() - t0))

    tprint('\nFinished %d/%d case(s).' % (len(new_rows), len(jobs)))
    if not new_rows:
        return
    summ = os.path.join(out_dir, 'summary.csv')
    old = pd.read_csv(summ) if os.path.isfile(summ) else pd.DataFrame()
    if len(old):
        bak = os.path.join(out_dir, 'summary_before_rerun_%s.csv' % time.strftime('%Y%m%d_%H%M%S'))
        shutil.copy2(summ, bak); tprint('Backed up summary.csv -> %s' % bak)
    new = pd.DataFrame(new_rows)
    if len(old):
        old = old[~old['case'].isin(new['case'])]
    df = pd.concat([old, new], ignore_index=True, sort=False).sort_values('case')
    if len(old):
        df = df[list(old.columns) + [c for c in df.columns if c not in old.columns]]
    df.to_csv(summ, index=False, float_format='%.8e')
    tprint('summary.csv now has %d cases (+%d)' % (len(df), len(new)))
    export_accepted_loads(summ, os.path.join(out_dir, 'accepted_loads.csv'))
    ur = [c for c in ['case', 'eps_1d_mid_post', 'UR1', 'UR2', 'UR3'] if c in df.columns]
    df[ur].to_csv(os.path.join(out_dir, 'accepted_1d_strain_ur.csv'), index=False, float_format='%.8e')
    plot_histograms(df, out_dir); plot_energy_scatter(df, out_dir); plot_strain_screen_audit(df, out_dir)
    tprint('Done. Plots/CSVs refreshed in %s' % out_dir)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
sweep_freebody.py -- STANDALONE.  Not part of the batch pipeline.

Runs freebody_cut.py over many case folders and concatenates the per-case
CSVs into one table, then prints how far apart the two routes ever got.

    python3 sweep_freebody.py --sample 40
    python3 sweep_freebody.py --limit 100 --workers 3
    python3 sweep_freebody.py                       # every case (slow)

BEFORE YOU RUN THE WHOLE FOLDER, read this
------------------------------------------
Route A here does exactly the summation reextract.py does, so sweeping all
~1000 cases with this script is a SECOND full pass over the same odbs for
numbers reextract.py is going to write into solid_resultants.csv anyway.
Each odb open costs tens of seconds; the full sweep is hours.

This script earns its keep as a VALIDATION, not as a production run: use
--sample to take cases spread evenly across the folder, confirm routes A
and B agree, and then trust reextract.py for the rest.  --sample 40 is
enough to see a systematic disagreement if one exists.

Idempotent: a case that already has freebody_cut.csv is skipped unless
--force, so an interrupted sweep just gets re-run.
"""
import os
import re
import sys
import csv
import glob
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor

# Defaults.  Every one of these can be overridden on the command line or
# through the environment, so nothing here needs editing to point the sweep
# at a different folder:
#     --cases-dir / first positional argument (required)
#     --batch-dir                             / $FB_BATCH_DIR
#     --abaqus                                / $FB_ABAQUS_BAT
BATCH_DIR = os.environ.get('FB_BATCH_DIR') or os.path.dirname(
    os.path.abspath(__file__))
ABAQUS_BAT = os.environ.get('FB_ABAQUS_BAT') or \
    '/mnt/c/SIMULIA/Commands/abaqus.bat'


def win(p):
    """WSL /mnt/X/... -> Windows X:\\...  (same rule as abaqus_io.win)"""
    p = os.path.abspath(p)
    m = re.match(r'^/mnt/([a-zA-Z])(/.*)?$', p)
    if m:
        return m.group(1).upper() + ':' + (m.group(2) or '').replace('/', '\\')
    return p


def wsl(p):
    """Accept a folder however the user happens to type it:
    'C:\\Users\\...', 'C:/Users/...', a WSL path, or a relative one.
    Windows drive paths are converted to /mnt/<drive>/... so the same
    argument works whether it was copied from Explorer or from a shell."""
    if p is None:
        return None
    p = p.strip().strip('"').strip("'")
    m = re.match(r'^([a-zA-Z]):[\\/](.*)$', p)
    if m:
        return '/mnt/' + m.group(1).lower() + '/' + m.group(2).replace('\\', '/')
    return os.path.abspath(os.path.expanduser(p))


def find_case_files(cdir):
    """<case>_3d.odb / .inp by the folder's own name, else whatever single
    *_3d.odb the folder holds -- so a folder that does not follow the
    case_XXXX convention still works."""
    cid = os.path.basename(cdir.rstrip('/'))
    odb = os.path.join(cdir, cid + '_3d.odb')
    inp = os.path.join(cdir, cid + '_3d.inp')
    if os.path.exists(odb) and os.path.exists(inp):
        return odb, inp
    cands = sorted(glob.glob(os.path.join(cdir, '*_3d.odb')))
    if len(cands) == 1:
        odb = cands[0]
        inp = odb[:-4] + '.inp'
        if os.path.exists(inp):
            return odb, inp
    return None, None


def run_case(args):
    cdir, batch_dir, stations, force, abaqus_bat = args
    cid = os.path.basename(cdir.rstrip('/'))
    odb, inp = find_case_files(cdir)
    out = os.path.join(cdir, 'freebody_cut.csv')
    if odb is None:
        return cid, 'no <case>_3d.odb + .inp pair in this folder'
    if os.path.exists(out) and not force:
        return cid, 'skipped'
    cmd = ['cmd.exe', '/c', win(abaqus_bat), 'python',
           win(os.path.join(batch_dir, 'freebody_cut.py')),
           win(odb), '--inp', win(inp), '--csv', win(out)]
    if stations:
        cmd += ['--stations', stations]
    try:
        r = subprocess.run(cmd, cwd=batch_dir, stdin=subprocess.DEVNULL,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=1800)
    except subprocess.TimeoutExpired:
        return cid, 'TIMEOUT'
    if r.returncode != 0 or not os.path.exists(out):
        tail = r.stdout.decode('utf-8', 'replace').strip().splitlines()[-3:]
        return cid, 'FAILED: ' + ' | '.join(tail)
    return cid, 'ok'


def main():
    ap = argparse.ArgumentParser(
        description='Sweep freebody_cut.py over a folder of case folders.',
        epilog='The folder may be given as the first positional argument, '
               'or with --cases-dir (one of them is required). Windows spellings '
               '(C:\\Users\\... or C:/Users/...) are accepted.')
    ap.add_argument('cases_dir_pos', nargs='?', default=None,
                    metavar='CASES_DIR',
                    help='folder holding the case_* subfolders')
    ap.add_argument('--cases-dir', default=None)
    ap.add_argument('--batch-dir', default=None,
                    help='folder holding freebody_cut.py '
                         '(default: this script\'s own folder)')
    ap.add_argument('--abaqus', default=None,
                    help='path to abaqus.bat (default: $FB_ABAQUS_BAT)')
    ap.add_argument('--stations', default=None,
                    help='comma-separated Z cuts; default is the script default')
    ap.add_argument('--sample', type=int, default=0,
                    help='take N cases spread evenly across the folder')
    ap.add_argument('--limit', type=int, default=0, help='first N cases only')
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--force', action='store_true', help='redo finished cases')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if not (a.cases_dir or a.cases_dir_pos):
        ap.error('give the cases folder, e.g.  python3 sweep_freebody.py '
                 '/mnt/c/.../batch_test/cases_half_ar_seed123')
    a.cases_dir = wsl(a.cases_dir or a.cases_dir_pos)
    a.batch_dir = wsl(a.batch_dir or BATCH_DIR)
    a.abaqus = wsl(a.abaqus or ABAQUS_BAT)
    if not os.path.isdir(a.cases_dir):
        sys.exit('cases dir not found: %s' % a.cases_dir)
    fb = os.path.join(a.batch_dir, 'freebody_cut.py')
    if not os.path.exists(fb):
        sys.exit('freebody_cut.py not found in %s' % a.batch_dir)
    if not os.path.exists(a.abaqus):
        sys.exit('abaqus.bat not found: %s  (set --abaqus or $FB_ABAQUS_BAT)'
                 % a.abaqus)
    print('cases  : %s' % a.cases_dir)
    print('script : %s' % fb)
    print('abaqus : %s' % a.abaqus)

    dirs = sorted(d for d in glob.glob(os.path.join(a.cases_dir, 'case_*'))
                  if os.path.isdir(d))
    if not dirs:                       # allow pointing straight at ONE case
        if find_case_files(a.cases_dir)[0]:
            dirs = [a.cases_dir]
        else:
            sys.exit('no case_* subfolders and no *_3d.odb in %s' % a.cases_dir)
    if a.sample and a.sample < len(dirs):
        step = len(dirs) / float(a.sample)
        dirs = [dirs[int(i * step)] for i in range(a.sample)]
    if a.limit:
        dirs = dirs[:a.limit]
    print('%d case(s), %d worker(s)' % (len(dirs), a.workers))

    jobs = [(d, a.batch_dir, a.stations, a.force, a.abaqus)
            for d in dirs]
    done = 0
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        for cid, status in ex.map(run_case, jobs):
            done += 1
            if status not in ('ok', 'skipped'):
                print('  [%4d/%4d] %s  %s' % (done, len(dirs), cid, status))
            elif done % 10 == 0 or done == len(dirs):
                print('  [%4d/%4d] %s' % (done, len(dirs), cid))

    # ---- concatenate -----------------------------------------------------
    out = a.out or os.path.join(a.cases_dir, 'freebody_cut_all.csv')
    n, hdr = 0, None
    with open(out, 'w', newline='') as fh:
        w = None
        for d in dirs:
            f = os.path.join(d, 'freebody_cut.csv')
            if not os.path.exists(f):
                continue
            with open(f) as g:
                rd = csv.reader(g)
                h = next(rd)
                if w is None:
                    hdr = ['case'] + h
                    w = csv.writer(fh)
                    w.writerow(hdr)
                for row in rd:
                    w.writerow([os.path.basename(d)] + row)
                    n += 1
    print('\nwrote %s  (%d rows)' % (out, n))

    # ---- how far apart did the two routes ever get? ----------------------
    if not n:
        return
    try:
        import pandas as pd
    except ImportError:
        print('(install pandas for the summary table)')
        return
    d = pd.read_csv(out)
    print('\nroute A minus route B, normalised by the per-case root scale')
    print('%6s %6s %12s %12s %12s' % ('z', 'n', 'median', 'p90', 'max'))
    for z, g in d.groupby('Z'):
        v = g.dM_resultant_norm.abs()
        print('%6.0f %6d %11.3e %11.3e %11.3e'
              % (z, len(g), v.median(), v.quantile(.9), v.max()))
    worst = d.loc[d.dM_resultant_norm.abs().idxmax()]
    print('\nworst row: %s at Z=%.0f, |dM|/scale = %.3e'
          % (worst.case, worst.Z, abs(worst.dM_resultant_norm)))


if __name__ == '__main__':
    main()

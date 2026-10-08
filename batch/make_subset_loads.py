#!/usr/bin/env python3
"""
make_subset_loads.py  (added 8 Oct 2026)
Pick a few cases from a finished batch and write a loads file that
batch_driver.py can replay with --loads-csv, e.g. to re-run only the worst
and the best cases with another 3D element type.

  python3 make_subset_loads.py --src .../cases_EI_1x1 --out .../loads_subset_C3D8I.csv \
        --cases 769,406,24,438,804 --worst 10 --best 10 --metric relL2_resultant_Mises

  --cases   explicit case numbers (always included)
  --worst N the N largest values of --metric (default Mises error)
  --best N  the N smallest values of --metric
Cases are written in the order: explicit, worst, best (no duplicates), as
rows load_id 0..n-1 with Vx, Vy, Mx and extra columns orig_case, group.
The replay run numbers its cases 0..n-1 in that order, so
compare_element_runs.py can pair them with the originals via orig_case.
Needs only pandas.
"""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--src', required=True, help='finished batch folder (needs summary.csv)')
ap.add_argument('--out', required=True, help='loads CSV to write')
ap.add_argument('--cases', default='', help='comma-separated case numbers')
ap.add_argument('--worst', type=int, default=0)
ap.add_argument('--best', type=int, default=0)
ap.add_argument('--metric', default='relL2_resultant_Mises')
a = ap.parse_args()

s = pd.read_csv(os.path.join(a.src, 'summary.csv'))
rows, seen = [], set()
def add(ids, group):
    for c in ids:
        c = int(c)
        if c in seen or c not in set(s['case']):
            if c not in seen:
                print('WARNING: case %d not in summary.csv, skipped' % c)
            continue
        seen.add(c); rows.append((c, group))
add([x for x in a.cases.split(',') if x.strip()], 'explicit')
add(s.nlargest(a.worst, a.metric)['case'] if a.worst else [], 'worst')
add(s.nsmallest(a.best, a.metric)['case'] if a.best else [], 'best')
if not rows:
    raise SystemExit('nothing selected: give --cases, --worst and/or --best')

sub = pd.DataFrame(rows, columns=['orig_case', 'group']).merge(
    s[['case', 'Vx', 'Vy', 'Mx']], left_on='orig_case', right_on='case').drop(columns='case')
sub.insert(0, 'load_id', range(len(sub)))
sub[['load_id', 'Vx', 'Vy', 'Mx', 'orig_case', 'group']].to_csv(a.out, index=False, float_format='%.10e')
print('wrote %d loads to %s' % (len(sub), a.out))
print(sub.groupby('group').size().to_string())
print('run with:  --cases %d --loads-csv %s' % (len(sub), a.out))

#!/usr/bin/env python3
"""
compare_element_runs.py  (added 8 Oct 2026)
Compare a replay run (e.g. C3D8I) with the original batch it was drawn from.

  python3 compare_element_runs.py --base .../cases_EI_1x1 --new .../cases_EI_1x1_C3D8I \
        --loads-csv .../loads_subset_C3D8I.csv

Pairs case i of --new with orig_case of row i in the loads CSV, prints the
errors side by side (median per group and per case) and writes
element_comparison.csv in --new. Needs only pandas/numpy.
"""
import argparse, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--base', required=True)
ap.add_argument('--new', required=True)
ap.add_argument('--loads-csv', required=True)
a = ap.parse_args()

L = pd.read_csv(a.loads_csv)
b = pd.read_csv(os.path.join(a.base, 'summary.csv')).set_index('case')
n = pd.read_csv(os.path.join(a.new, 'summary.csv')).set_index('case')
M = [('Mises', 'relL2_resultant_Mises'), ('I1', 'relL2_resultant_inv_1'), ('S33', 'relL2_Szz_S33'),
     ('S13', 'relL2_Szz_S13'), ('S23', 'relL2_Szz_S23'), ('Energy', 'relL2_Energy')]
rows = []
for r in L.itertuples():
    if r.load_id not in n.index or r.orig_case not in b.index:
        continue
    d = {'orig_case': r.orig_case, 'new_case': r.load_id, 'group': getattr(r, 'group', '')}
    d['twist'] = abs(b.loc[r.orig_case, 'UR3']) / np.hypot(b.loc[r.orig_case, 'UR1'], b.loc[r.orig_case, 'UR2'])
    # sanity: the loads must be identical
    d['load_mismatch'] = max(abs(n.loc[r.load_id, c] - b.loc[r.orig_case, c]) for c in ['Vx', 'Vy', 'Mx'])
    for lab, col in M:
        d[lab + '_base_%'] = 100 * b.loc[r.orig_case, col]
        d[lab + '_new_%'] = 100 * n.loc[r.load_id, col]
    rows.append(d)
D = pd.DataFrame(rows)
if not len(D):
    raise SystemExit('no paired cases found: is the new run finished (summary.csv)?')
out = os.path.join(a.new, 'element_comparison.csv'); D.to_csv(out, index=False)
pd.set_option('display.width', 220)
print('paired %d cases; largest load mismatch %.2e' % (len(D), D.load_mismatch.max()))
cols = []
for lab, _ in M:
    cols += [lab + '_base_%', lab + '_new_%']
print('\nmedian over cases (percent):')
print(D.groupby('group')[cols].median().round(2).T.to_string())
print('\nper case, Mises and S33 (percent):')
print(D[['orig_case', 'group', 'twist', 'Mises_base_%', 'Mises_new_%', 'S33_base_%', 'S33_new_%']].round(3).to_string(index=False))
hi, lo = D[D.twist > D.twist.median()], D[D.twist <= D.twist.median()]
print('\ntwist trend (median Mises error, percent): base low/high twist %.2f / %.2f ; new %.2f / %.2f'
      % (lo['Mises_base_%'].median(), hi['Mises_base_%'].median(), lo['Mises_new_%'].median(), hi['Mises_new_%'].median()))
print('written', out)

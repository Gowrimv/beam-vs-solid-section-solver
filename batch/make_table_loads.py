#!/usr/bin/env python3
"""
make_table_loads.py  (added 7 Oct 2026)
=======================================
Generate ONE set of non-dimensional loads exactly as the gate table says
(Case A: h = 1, L = 100), then write it in physical units for every section
size, so the 1x1 and 0.5x0.5 runs get identical non-dimensional problems.

  Generate:  m_x,tip ~ U(-0.7, 0.7);  v_x, v_y ~ U(-0.6, 0.6)
  Strain gate:   (1/200) max(|m - v_y| + |v_x|, |m|) <= eps_tol (0.005)
  Bending gate:  sqrt(v_x^2+v_y^2) / (400 (|m - v_y/2| + |v_x/2|)) <= 0.1
  Midspan floor: (1/200)(|m - v_y/2| + |v_x/2|) >= 0.1 eps_tol
  Keep: draw again until N samples pass all three.
  Apply: Mx = M0 m,  Vx = F0 v_x,  Vy = F0 v_y   (F0 = E I / L^2, M0 = F0 L)

The gates are evaluated for Case A (h = 1). At h = 0.5 the same
non-dimensional loads give half the strain and half the shear/bending
ratio, so every row also passes there.

Writes N + spares rows; batch_driver.py uses the rows beyond --cases to
replace cases rejected after the run (same rows at both sizes as long as
the same cases are rejected).

  python3 make_table_loads.py --n 1000 --spares 300 --seed 123 \\
      --out-dir /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test
"""
import argparse, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--n', type=int, default=1000)
ap.add_argument('--spares', type=int, default=300)
ap.add_argument('--seed', type=int, default=123)
ap.add_argument('--eps-tol', type=float, default=0.005)
ap.add_argument('--E', type=float, default=1.2e5)
ap.add_argument('--L', type=float, default=100.0)
ap.add_argument('--sides', type=str, default='1,0.5')
ap.add_argument('--out-dir', type=str, default='.')
a = ap.parse_args()

rng = np.random.default_rng(a.seed)
need = a.n + a.spares
rows, drawn = [], 0
while len(rows) < need:
    k = 4096
    m  = rng.uniform(-0.7, 0.7, k)
    vx = rng.uniform(-0.6, 0.6, k)
    vy = rng.uniform(-0.6, 0.6, k)
    drawn += k
    strain = (1/200.) * np.maximum(np.abs(m - vy) + np.abs(vx), np.abs(m))
    mmid   = np.abs(m - vy/2) + np.abs(vx/2)
    with np.errstate(divide='ignore'):
        bend = np.hypot(vx, vy) / (400. * mmid)
    ok = (strain <= a.eps_tol) & (bend <= 0.1) & ((1/200.) * mmid >= 0.1 * a.eps_tol)
    for i in np.where(ok)[0]:
        rows.append((m[i], vx[i], vy[i], strain[i], bend[i], mmid[i]))
        if len(rows) >= need:
            break
nd = pd.DataFrame(rows, columns=['m_tip', 'v_x', 'v_y', 'strain_caseA', 'shear_bending_caseA', 'm_mid'])
nd.insert(0, 'load_id', range(len(nd)))
p = os.path.join(a.out_dir, 'loads_table_nd.csv')
nd.to_csv(p, index=False, float_format='%.10e')
print('drawn %d, kept %d (pass rate %.1f%%) -> %s' % (drawn, len(nd), 100.*len(nd)/drawn, p))

# 7 Oct 2026: E per section = E_REF / h^4 (E*I fixed = 1e4), so F0 = 1 and
# M0 = 100 at every size and ONE physical loads file serves all runs.
# loads_table.csv is that file; the per-size files below are identical to it.
for s in [float(x) for x in a.sides.split(',')]:
    I = s**4 / 12.0
    E_s = a.E / s**4
    F0 = E_s * I / a.L**2
    M0 = F0 * a.L
    out = pd.DataFrame({'load_id': nd.load_id, 'Vx': nd.v_x * F0, 'Vy': nd.v_y * F0, 'Mx': nd.m_tip * M0})
    tag = '1x1' if s == 1 else ('half' if s == 0.5 else 's%g' % s)
    p = os.path.join(a.out_dir, 'loads_table_%s.csv' % tag)
    out.to_csv(p, index=False, float_format='%.10e')
    print('side %g: E=%.6g F0=%.6g M0=%.6g -> %s' % (s, E_s, F0, M0, p))
    if s == 1:
        out.to_csv(os.path.join(a.out_dir, 'loads_table.csv'), index=False, float_format='%.10e')

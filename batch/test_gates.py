#!/usr/bin/env python3
"""Tests for the three gates that need no strains (added 7-8 Oct 2026):

  1. 1D convergence   abaqus_io.residual_check_1d()  -> DISCARD_UNCONVERGED_1D
                      (driver: res_ratio_1d > --max-1d-residual, default 0.01)
  2. 3D hourglass     ALLAE/ALLIE from solid_energy_history.csv -> DISCARD_HOURGLASS
                      (driver: ratio > --max-hourglass, default config.HOURGLASS_MAX)
  3. mid-moment floor strain_checks.mid_moment_nd()  (pre-run, on the loads)
                      (driver: m_mid_nd < --min-mid-moment-nd, default config.MID_FLOOR_ND)

Part A  synthetic unit tests (no Abaqus, seconds): known-pass and known-fail
        inputs for each gate, including Abaqus's 1e-2 time-average fallback.
Part B  --scan FOLDER [...]: apply gates 1 and 2 to existing case folders
        (reads *_1d.msg tails and solid_energy_history.csv) and list the cases
        the current thresholds would discard. Older folders give real
        unconverged 1D cases to check gate 1 against.

Usage:
  python3 test_gates.py                       # Part A only
  python3 test_gates.py --scan DIR [DIR ...]  # Part A + scan
"""
import os, sys, glob, tempfile, argparse
import pandas as pd

import config
from abaqus_io import residual_check_1d
from strain_checks import mid_moment_nd, read_case_diagnostics, table_strain_gate, table_bending_ratio

MAX_RES = 0.01
MAX_HG = config.HOURGLASS_MAX
FLOOR = config.MID_FLOOR_ND
fails = 0

def check(name, cond, detail=''):
    global fails
    print('  [%s] %s %s' % ('PASS' if cond else 'FAIL', name, detail))
    if not cond:
        fails += 1

MSG = """
          AVERAGE FORCE                     %s       TIME AVG. FORCE        %s
          LARGEST RESIDUAL FORCE            %s   AT NODE        101   DOF  1
          AVERAGE MOMENT                      11.7       TIME AVG. MOMENT       11.7
          LARGEST RESIDUAL MOMENT           %s   AT NODE        101   DOF  5
"""

def part_a():
    print('Part A: synthetic tests')
    vx, vy, mx = -0.40755, 0.46625, 25.529          # pilot case 0 (F0=1, M0=100)
    with tempfile.TemporaryDirectory() as d:
        # 1. convergence
        for label, tavg, rf, rm, expect_fail in [
                ('converged (pilot case 0 values)', '0.292', '1.088E-04', '2.172E-03', False),
                ('1e-2 fallback, big residual',      '1.000E-02', '5.000E-02', '3.000E+00', True),
                ('residual just above 1% of |V|',    '0.292', '6.300E-03', '1.000E-03', True)]:
            with open(os.path.join(d, 't_1d.msg'), 'w') as f:
                f.write(MSG % (tavg, tavg, rf, rm))
            r = residual_check_1d(d, 't_1d', vx, vy, mx, config.L)
            ratio = r['res_ratio_1d']
            check('gate 1: ' + label, ratio is not None and ((ratio > MAX_RES) == expect_fail),
                  'ratio=%.3g fallback=%s -> %s' % (ratio, r['tavg_fallback_1d'],
                                                   'discard' if ratio > MAX_RES else 'keep'))
        os.remove(os.path.join(d, 't_1d.msg'))
        r = residual_check_1d(d, 't_1d', vx, vy, mx, config.L)
        check('gate 1: missing .msg -> no decision (None)', r['res_ratio_1d'] is None)

        # 2. hourglass
        for label, allae, allie, expect_fail in [
                ('pilot case 0 (0.17%)', 6.58e-3, 3.775, False),
                ('5% hourglass',          0.19,   3.775, True)]:
            pd.DataFrame([{'ALLAE': allae, 'ALLIE': allie, 'ALLSE': allie, 'ALLWK': allie}]
                         ).to_csv(os.path.join(d, 'solid_energy_history.csv'), index=False)
            hg = read_case_diagnostics(d).get('ALLAE_over_ALLIE')
            check('gate 2: ' + label, hg is not None and ((hg > MAX_HG) == expect_fail),
                  'ALLAE/ALLIE=%.3g -> %s' % (hg, 'discard' if hg > MAX_HG else 'keep'))

    # 3. mid-moment floor: Mx chosen so Mx(z=Z_MID) = Mx_tip - Vy (L - Z_MID) ~ 0
    arm = config.L - config.Z_MID
    vy = 0.4; vx = 0.01
    for label, mx, expect_fail in [('moment crosses zero at mid-span', vy * arm, True),
                                   ('pilot case 0', 25.529, False)]:
        vxx, vyy = (vx, vy) if expect_fail else (-0.40755, 0.46625)
        m = mid_moment_nd(vxx, vyy, mx)
        check('gate 3: ' + label, (m < FLOOR) == expect_fail,
              'm_mid_nd=%.3g (floor %.2g) -> %s' % (m, FLOOR, 'discard' if m < FLOOR else 'keep'))
    # also report whether that designed row would pass the other pre-run gates,
    # so it really reaches the floor check in the driver
    print('      designed floor row Vx=%.3g Vy=%.3g Mx=%.4g: strain gate %.3g (limit %.3g), '
          'shear/bending %.3g (limit %.3g)'
          % (vx, vy, vy * arm, table_strain_gate(vx, vy, vy * arm), config.STRAIN_THRESHOLD,
             table_bending_ratio(vx, vy, vy * arm), config.BENDING_GATE_MAX))

def loads_for(folder):
    for fn in ('summary.csv', 'accepted_loads.csv'):
        p = os.path.join(folder, fn)
        if os.path.isfile(p):
            s = pd.read_csv(p)
            if {'case', 'Vx', 'Vy', 'Mx'} <= set(s.columns):
                return s.set_index('case')[['Vx', 'Vy', 'Mx']]
    return None

def part_b(folders):
    print('\nPart B: scan of existing folders (limits: res %.3g, ALLAE/ALLIE %.3g)' % (MAX_RES, MAX_HG))
    for folder in folders:
        L = loads_for(folder)
        rows = []
        for cd in sorted(glob.glob(os.path.join(folder, 'case_*'))):
            try:
                cid = int(os.path.basename(cd).split('_')[1])
            except ValueError:
                continue
            if L is None or cid not in L.index:
                continue
            vx, vy, mx = L.loc[cid]
            r = residual_check_1d(cd, 'case_%04d_1d' % cid, vx, vy, mx, config.L)
            hg = read_case_diagnostics(cd).get('ALLAE_over_ALLIE')
            rows.append({'case': cid, 'res_ratio_1d': r['res_ratio_1d'],
                         'fallback': r['tavg_fallback_1d'], 'ALLAE_over_ALLIE': hg})
        df = pd.DataFrame(rows)
        if not len(df):
            print('  %s: no cases with loads' % folder); continue
        bad1 = df[df.res_ratio_1d > MAX_RES]; bad2 = df[df.ALLAE_over_ALLIE.astype(float) > MAX_HG]
        print('  %s: %d cases | res read for %d, fallback in %d, discard(1D) %d | '
              'hourglass read for %d, discard(HG) %d'
              % (os.path.basename(folder.rstrip('/')), len(df), df.res_ratio_1d.notna().sum(),
                 int((df.fallback == True).sum()), len(bad1), df.ALLAE_over_ALLIE.notna().sum(), len(bad2)))
        if len(bad1):
            print('    worst 1D:', ', '.join('%d (%.3g%s)' % (r.case, r.res_ratio_1d, ', fallback' if r.fallback else '')
                                         for r in bad1.sort_values('res_ratio_1d', ascending=False).head(10).itertuples()))
        out = os.path.join(folder, 'gate_scan.csv')
        df.to_csv(out, index=False)
        print('    -> %s' % out)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--scan', nargs='*', default=[])
    a = ap.parse_args()
    part_a()
    if a.scan:
        part_b(a.scan)
    print('\n%s' % ('ALL SYNTHETIC TESTS PASSED' if fails == 0 else '%d TEST(S) FAILED' % fails))
    sys.exit(1 if fails else 0)

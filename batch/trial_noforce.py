#!/usr/bin/env python3
"""
trial_noforce.py  (7 Oct 2026) -- TRIAL ONLY, does not touch batch results
=========================================================================
Re-solves FEniCS for existing case folders with warping_core_trial.py:
  --force-constraints off : strains only (alpha/beta resultant constraints
                            removed), i.e. no force targets
  --rotation-constraint   : 'product' (current u1*u2) or 'rotation'
                            ((X1-C1) u2 - (X2-C2) u1)
Each mode writes into <case>/trial_<mode>/ (VTU, section_X0, errors.csv,
section_frame_check.csv, overlays) and a table of all cases/modes vs the
original batch result goes to <cases-dir>/trial_summary.csv.

Uses the same inputs as batch_driver.py: 1D strains at Z_MID pulled back with
R_local, shear-strain correction per config, FEniCS at E = FENICS_E with the
stress outputs rescaled by E/FENICS_E.

  cd /mnt/c/Users/macLab/Desktop/FEniCSfiles/batch
  BATCH_SECTION_SIDE=1 python3 trial_noforce.py \
      --cases-dir /mnt/c/Users/macLab/Documents/AbaqusProjects/batch_test/pilot_table_1x1
"""
import os, glob, shutil, argparse
import numpy as np
import pandas as pd

import config
from config import tprint
import warping_core_trial as wc
from dolfin import assemble, dx, tr
from compare import _compare_inprocess

STRESS_FIELDS = ('S11', 'S22', 'S33', 'S12', 'S13', 'S23', 'Mises', 'inv_1')


def rescale(sub, base, k):
    import pyvista as pv
    for suf in ('_ref.vtu', '_deformed.vtu'):
        p = os.path.join(sub, base + suf)
        if os.path.isfile(p):
            m = pv.read(p)
            for name in STRESS_FIELDS:
                if name in m.point_data.keys():
                    m.point_data[name] = m.point_data[name] * k
            m.save(p)
    p = os.path.join(sub, base + '_section_X0.csv')
    if os.path.isfile(p):
        df = pd.read_csv(p)
        for name in STRESS_FIELDS:
            if name in df.columns:
                df[name] = df[name] * k
        df.to_csv(p, index=False, float_format='%.6e')


def inputs(case_dir, pullback=True):
    b = pd.read_csv(os.path.join(case_dir, 'beam_whole.csv'))
    b.columns = b.columns.str.strip()
    r = b.iloc[(b['Z'] - config.Z_MID).abs().idxmin()]
    sv = np.array([float(r['SE3']), float(r['SE2']), float(r['SE1'])])
    if getattr(config, 'CORRECT_SHEAR_STRAINS', False):
        K = config.SHEAR_K_FACTOR * config.G * config.A
        sv = np.array([float(r['SF3']) / K, float(r['SF2']) / K, float(r['SE1'])])
    sk = np.array([float(r['SK1']), float(r['SK2']), float(r['SK3'])])
    R = wc.rotation_matrix_from_ur(float(r['UR1']), float(r['UR2']), float(r['UR3']))
    n_t = np.array([float(r['SF3']), float(r['SF2']), float(r['SF1'])])
    m_t = np.array([float(r['SM1']), float(r['SM2']), float(r['SM3'])])
    if pullback:                       # OLD batch_driver behaviour (double rotation, pre 2026-10-08)
        return R, R.T @ sv, R.T @ sk, n_t, m_t
    return R, sv, sk, n_t, m_t         # Abaqus SE/SK are already beam-local


def scalars(d):
    out = {}
    e = os.path.join(d, 'errors.csv')
    if os.path.isfile(e):
        df = pd.read_csv(e).set_index('component')
        for c in ('Mises', 'inv_1', 'S33', 'S13', 'S23'):
            if c in df.index:
                out['relL2_' + c] = df.loc[c, 'relL2']
                if 'relL2_Szz' in df.columns:
                    out['relL2_Szz_' + c] = df.loc[c, 'relL2_Szz']
        for c in ('N_resultant', 'M_resultant'):
            if c in df.index:
                out['relL2_' + c] = df.loc[c, 'relL2']
    s = os.path.join(d, 'section_frame_check.csv')
    if os.path.isfile(s):
        out.update(pd.read_csv(s).iloc[0].to_dict())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cases-dir', required=True)
    ap.add_argument('--cases', default=None, help='comma-separated ids (default: all in summary.csv)')
    ap.add_argument('--strain-frame', choices=['pullback', 'local'], default='local',
                    help="pullback: R^T applied to SE/SK (current batch); local: SE/SK used as-is")
    ap.add_argument('--modes', default='off:product,off:rotation',
                    help="comma list of <force on|off|MASK>:<rotation product|rotation>; "
                         "MASK = 6 digits n1n2n3m1m2m3, 1 = keep that resultant constraint, "
                         "e.g. 000111:rotation keeps only the moment constraints")
    a = ap.parse_args()
    cdir = os.path.abspath(a.cases_dir)
    ids = ([int(x) for x in a.cases.split(',')] if a.cases
           else list(pd.read_csv(os.path.join(cdir, 'summary.csv'))['case']))
    k = config.E / config.FENICS_E
    rows = []
    for c in ids:
        rows.append(dict(case=c, mode='original (forces on, product)',
                         **scalars(os.path.join(cdir, 'case_%04d' % c))))
    for mode in a.modes.split(','):
        fc, rc = mode.split(':')
        mask = fc if (len(fc) == 6 and set(fc) <= set('01')) else None
        if mask:
            tag = 'trial_mask%s_%s' % (mask, rc)
        else:
            tag = 'trial_%s_%s' % ('force' if fc == 'on' else 'noforce', rc)
        if a.strain_frame == 'local':
            tag += '_localstrain'
        tprint('=== mode %s ===' % tag)
        mesh, _ = wc.load_mesh(config.MESH_XMF, workdir=cdir)
        P = wc.build_problem(mesh, E=config.FENICS_E, nu_val=config.NU, L_scale=1.0,
                             force_constraints=(fc == 'on' or mask is not None),
                             rotation_constraint=rc, constraint_mask=mask)
        for c in ids:
            case_dir = os.path.join(cdir, 'case_%04d' % c)
            base = 'case_%04d' % c
            sub = os.path.join(case_dir, tag)
            os.makedirs(sub, exist_ok=True)
            try:
                R, sv, sk, n_t, m_t = inputs(case_dir, pullback=(a.strain_frame == 'pullback'))
                wc.assign_targets(P, n_t / k, m_t / k)
                wc.forward_solve(P, sv, sk, warm=False, n_steps=10)
                Phi = k * float(assemble(((P['lmbda_param']/2)*tr(P['E_strain'])**2
                                          + P['mu_param']*tr(P['E_strain']*P['E_strain']))*dx))
                pp = wc.postprocess(P, sub, base, R_local=R, R_out=np.eye(3),
                                    do_plots=False, do_vtu=True)
                parts = P['u'].split(deepcopy=True)
                ab = [parts[i].vector().get_local()[0] for i in range(7, 13)]
                rescale(sub, base, k)
                for f in ('solid_energy.csv',):
                    if os.path.isfile(os.path.join(case_dir, f)):
                        shutil.copy2(os.path.join(case_dir, f), os.path.join(sub, f))
                n_out = k * np.asarray(pp['n0_out']); m_out = k * np.asarray(pp['m0_out'])
                _compare_inprocess(sub, os.path.join(case_dir, 'solid_mid.csv'),
                                   os.path.join(sub, base + '_ref.vtu'),
                                   n_out, m_out, R @ n_t, R @ m_t, Phi, R=R)
                row = dict(case=c, mode=tag, **scalars(sub))
                row.update({'alpha1': ab[0], 'alpha2': ab[1], 'alpha3': ab[2],
                            'beta1': ab[3], 'beta2': ab[4], 'beta3': ab[5]})
                # how far the strain-driven resultants are from the 1D ones (F0/M0 units)
                row['dN_over_F0'] = float(np.linalg.norm(n_out - R @ n_t)) / config.FORCE_SCALE
                row['dM_over_M0'] = float(np.linalg.norm(m_out - R @ m_t)) / config.MOMENT_SCALE
                rows.append(row)
                tprint('[%04d] %s done' % (c, tag))
            except Exception as ex:
                tprint('[%04d] %s FAILED: %s' % (c, tag, ex))
                rows.append(dict(case=c, mode=tag, error=str(ex)))
    df = pd.DataFrame(rows).sort_values(['case', 'mode'])
    out = os.path.join(cdir, 'trial_summary.csv')
    df.to_csv(out, index=False, float_format='%.5g')
    cols = [x for x in ['case', 'mode', 'transverse_ratio_fen', 'transverse_ratio_3d', 'relL2_Mises',
                        'relL2_inv_1', 'relL2_Szz_S13', 'relL2_Szz_S23', 'dN_over_F0', 'dM_over_M0', 'error']
            if x in df.columns]
    pd.set_option('display.width', 250)
    print(df[cols].to_string(index=False))
    print('\nWritten', out)


if __name__ == '__main__':
    main()

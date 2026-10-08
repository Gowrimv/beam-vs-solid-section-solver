"""
compare.py
==========
The in-process FEniCS-vs-Abaqus-3D comparison stage that used to be a
separate compare_section.py subprocess: stress invariants/errors.csv,
the six-component contour plots, the stress centerline overlay, the NEW
strain centerline overlay (strain_overlay.png -- LE33/LE13/LE23, Abaqus
solid vs FEniCS, along the X=0 section line), and force_error.csv.
Needs pyvista (to read the FEniCS *_ref.vtu) and matplotlib -- the one
other module in the split (besides fenics_solve.py) that isn't safely
importable without the full simulation stack installed.
"""
import os
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyvista as pv
from scipy.spatial import cKDTree

import config
from config import tprint
from strain_checks import invariants


def _compare_inprocess(case_dir, solid_csv, vtu_path, n_out, m_out,
                       abq_n=None, abq_m=None, phi_fenics=None, R=None):
    abq_df = pd.read_csv(solid_csv); abq_df.columns = abq_df.columns.str.strip()
    fen = pv.read(vtu_path)

    STRESS_COMPS = ['S11','S22','S33','S12','S13','S23']

    # ── DIRECT node matching, not interpolation ─────────────────────────
    # slice_50.xmf (FEniCS mesh) was sliced directly from the 3D solid, so
    # its 121 nodes are the SAME physical nodes as solid_mid.csv, at the
    # same (X,Y). Match each Abaqus node to its nearest FEniCS node
    # instead of interpolating — with coincident meshes this is exact.
    abq_xy = abq_df[['X','Y']].values.astype(float)
    fen_xy = np.array(fen.points)[:, :2]

    tree = cKDTree(fen_xy)
    dist, idx = tree.query(abq_xy)
    max_dist = float(dist.max())
    if max_dist > 1e-4:
        tprint('  WARNING: node match max distance %.3e — meshes may not '
               'coincide as expected' % max_dist)

    def af(col): return abq_df[col].values.astype(float) if col in abq_df.columns \
                  else np.zeros(len(abq_df))
    fen_s = {c: np.asarray(fen.point_data[c], float)[idx]
             for c in STRESS_COMPS if c in fen.point_data}
    abq_s = {c: af(c) for c in STRESS_COMPS}
    fen_inv = invariants(fen_s); abq_inv = invariants(abq_s)

    def l2(a): return float(np.sqrt(np.mean(a**2)))

    # COMMON STRESS SCALE. relL2 divides each component by its OWN norm,
    # which collapses for components that are near zero over the section --
    # S13/S23 span a factor of 28 across cases for that reason, and J3,
    # being CUBIC in stress, carries roughly 3x the relative error of any
    # degree-1 quantity purely from error propagation through a product.
    # sigma_ref is one scale per case (the 3D von Mises norm), so
    # relL2_resultant is comparable ACROSS components. J3 uses sigma_ref**3
    # to stay dimensionally consistent with its stress^3 units.
    sigma_ref = l2(abq_inv['Mises'])
    if not np.isfinite(sigma_ref) or sigma_ref <= 1e-20:
        sigma_ref = np.nan
    # Second per-case common scale: the 3D axial stress norm. relL2_Szz is
    # ||dS_ij|| / ||S33_3D|| with ONE denominator for every component, so
    # shear errors are judged against the governing bending stress (not
    # against their own near-zero norm). Added 7 Oct 2026.
    szz_ref = l2(abq_s['S33'])
    if not np.isfinite(szz_ref) or szz_ref <= 1e-20:
        szz_ref = np.nan

    rows = []
    for comp in STRESS_COMPS:
        diff = fen_s[comp]-abq_s[comp]; ref = l2(abq_s[comp])
        rows.append({'component':comp, 'L2':l2(diff), 'Linf':float(np.max(np.abs(diff))),
                     'relL2': l2(diff)/ref if ref>1e-20 else np.nan,
                     'relL2_resultant': l2(diff)/sigma_ref,
                     'relL2_Szz': l2(diff)/szz_ref,
                     'L2_over_scale': l2(diff)/config.STRESS_SCALE})
    for comp in ['Mises','inv_1']:   # stress-unit invariants
        diff = fen_inv[comp]-abq_inv[comp]; ref = l2(abq_inv[comp])
        rows.append({'component':comp, 'L2':l2(diff), 'Linf':float(np.max(np.abs(diff))),
                     'relL2': l2(diff)/ref if ref>1e-20 else np.nan,
                     'relL2_resultant': l2(diff)/sigma_ref,
                     'relL2_Szz': l2(diff)/szz_ref,
                     'L2_over_scale': l2(diff)/config.STRESS_SCALE})
    # J3 = det(deviatoric stress) — stress^3 units, normalized by
    # config.J3_SCALE. No NaN masking needed (unlike the old Lode angle):
    # a determinant is defined everywhere, including at hydrostatic
    # states (J2~0) where cos(3*theta) was undefined.
    fa = fen_inv['J3']; aa = abq_inv['J3']
    diff = fa-aa; ref = l2(aa)
    rows.append({'component':'J3', 'L2':l2(diff),
                 'Linf':float(np.max(np.abs(diff))),
                 'relL2': l2(diff)/ref if ref>1e-20 else np.nan,
                 # sigma_ref**3: J3 = det(dev(sigma)) has stress^3 units, so
                 # this puts it on the same footing as the degree-1 entries
                 # and divides out the ~3x cubic amplification.
                 'relL2_resultant': l2(diff)/(sigma_ref**3),
                 'relL2_Szz': l2(diff)/(szz_ref**3),
                 'L2_over_scale': l2(diff)/config.J3_SCALE})

    if abq_n is not None and abq_m is not None:
        n_diff = np.asarray(n_out) - np.asarray(abq_n)
        m_diff = np.asarray(m_out) - np.asarray(abq_m)
        n_ref  = float(np.linalg.norm(abq_n))
        m_ref  = float(np.linalg.norm(abq_m))
        rows.append({'component':'N_resultant',
                     'L2': float(np.linalg.norm(n_diff)),
                     'Linf': float(np.max(np.abs(n_diff))),
                     'relL2': float(np.linalg.norm(n_diff))/n_ref if n_ref>1e-20 else np.nan,
                     'L2_over_scale': float(np.linalg.norm(n_diff))/config.FORCE_SCALE})
        rows.append({'component':'M_resultant',
                     'L2': float(np.linalg.norm(m_diff)),
                     'Linf': float(np.max(np.abs(m_diff))),
                     'relL2': float(np.linalg.norm(m_diff))/m_ref if m_ref>1e-20 else np.nan,
                     'L2_over_scale': float(np.linalg.norm(m_diff))/config.MOMENT_SCALE})

    energy_path = os.path.join(case_dir, 'solid_energy.csv')
    if phi_fenics is not None and not np.isnan(phi_fenics) and os.path.isfile(energy_path):
        else_per_len = float(pd.read_csv(energy_path)['U_per_length_ELSE'].iloc[0])
        e_diff = phi_fenics - else_per_len
        rows.append({'component':'Energy',
                     'L2': abs(e_diff), 'Linf': abs(e_diff),
                     'relL2': abs(e_diff)/abs(else_per_len) if abs(else_per_len)>1e-20 else np.nan,
                     'L2_over_scale': abs(e_diff)/config.ENERGY_SCALE})

    pd.DataFrame(rows).to_csv(os.path.join(case_dir,'errors.csv'),
                              index=False, float_format='%.6e')

    # ── Section-frame regime check (added 7 Oct 2026) ───────────────────
    # Beam theory assumes the in-plane stresses S11, S22, S12 are small
    # compared with the axial stress, IN THE SECTION'S OWN (rotated) FRAME.
    # Pull the global stresses back with R (the beam rotation at Z_MID,
    # local -> global) and report RMS(sqrt(S11^2+S22^2+2 S12^2)) / RMS(S33)
    # for the 3D solid and for FEniCS. R=None -> global axes (flagged).
    try:
        Rm = np.eye(3) if R is None else np.asarray(R, float)
        def _tr(s):
            T = np.array([[s['S11'], s['S12'], s['S13']],
                          [s['S12'], s['S22'], s['S23']],
                          [s['S13'], s['S23'], s['S33']]])      # (3,3,n)
            loc = np.einsum('ji,jkn,kl->iln', Rm, T, Rm)          # R^T T R
            inplane = np.sqrt(loc[0,0]**2 + loc[1,1]**2 + 2*loc[0,1]**2)
            den = l2(loc[2,2])
            return l2(inplane)/den if den > 1e-20 else np.nan
        pd.DataFrame([{
            'transverse_ratio_3d': _tr(abq_s),
            'transverse_ratio_fen': _tr(fen_s) if all(c in fen_s for c in STRESS_COMPS) else np.nan,
            'section_frame': 'identity' if R is None else 'beam_UR_at_Z_MID',
        }]).to_csv(os.path.join(case_dir, 'section_frame_check.csv'),
                   index=False, float_format='%.6e')
    except Exception as ex:
        tprint('  WARNING: section-frame check failed: %s' % ex)

    # contour plots
    x,y = abq_df['X'].values, abq_df['Y'].values
    plot_comps = ['S11','S22','S33','S12','S13','S23','Mises','inv_1','J3']   # all components (7 Oct 2026)
    fig,axes = plt.subplots(len(plot_comps),3,
                            figsize=(13,3.8*len(plot_comps)),squeeze=False)
    for r,comp in enumerate(plot_comps):
        fa = fen_inv.get(comp, fen_s.get(comp, np.zeros(len(x))))
        aa = abq_inv.get(comp, abq_s.get(comp, np.zeros(len(x))))
        if comp == 'J3':
            scale, label = config.J3_SCALE, 'J3 / \u03c3\u2080\u00b3'
        elif comp in config.STRESS_UNIT_COMPS:
            scale, label = config.STRESS_SCALE, ('%s / \u03c3\u2080' % comp)
        else:
            scale, label = 1.0, comp
        fa = np.asarray(fa) / scale
        aa = np.asarray(aa) / scale
        diff = fa-aa
        vmin,vmax = min(fa.min(),aa.min()), max(fa.max(),aa.max())
        if vmin==vmax: vmin-=1e-12; vmax+=1e-12
        for col,(data,ttl) in enumerate([(aa,'Abaqus'),(fa,'FEniCS'),(diff,'Diff')]):
            ax = axes[r,col]
            lv = np.linspace(-max(abs(diff.min()),abs(diff.max()))+1e-12,
                              max(abs(diff.min()),abs(diff.max()))+1e-12,24) \
                 if col==2 else np.linspace(vmin,vmax,24)
            sc = ax.tricontourf(x,y,data,levels=lv,
                                cmap='coolwarm' if col==2 else 'viridis',
                                extend='both')
            ax.set_aspect('equal')
            ax.set_title('%s %s'%(label,ttl),fontsize=8)
            fig.colorbar(sc,ax=ax,shrink=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(case_dir,'contours.png'),dpi=110)
    plt.close(fig)

    # ── stress overlay along Y at section centreline (X≈0) ──────────────
    cl_comps = ['S11','S22','S33','S12','S13','S23','Mises','inv_1']   # all components (7 Oct 2026)
    cl_files  = glob.glob(os.path.join(case_dir,'*_section_X0.csv'))

    x_abq_uniq = np.unique(np.round(abq_df['X'].values, 8))
    x_cl_abq   = x_abq_uniq[np.argmin(np.abs(x_abq_uniq))]
    abq_cl     = abq_df[np.abs(abq_df['X'].values - x_cl_abq) < 1e-6].sort_values('Y')

    fen_cl = None
    if cl_files:
        fen_cl = pd.read_csv(cl_files[0])
        fen_cl.columns = fen_cl.columns.str.strip()

    fig, axes = plt.subplots(2, 4, figsize=(18, 9), squeeze=False)
    axes = axes.flatten()
    for ax, comp in zip(axes, cl_comps):
        plotted = False
        if fen_cl is not None and comp in fen_cl.columns and 'x2' in fen_cl.columns:
            ax.plot(fen_cl[comp].values/config.STRESS_SCALE, fen_cl['x2'].values,
                    '-', lw=1.8, color='#2471a3', label='FEniCS (warping)')
            plotted = True
        if comp in abq_cl.columns and len(abq_cl) > 0:
            ax.plot(abq_cl[comp].values/config.STRESS_SCALE, abq_cl['Y'].values,
                    '--', lw=1.8, color='#c0392b', label='Abaqus 3D')
            plotted = True
        ax.axvline(0, color='gray', lw=0.6, ls=':')
        ax.set_title('%s / \u03c3\u2080' % comp, fontsize=10)
        ax.set_xlabel('Stress / \u03c3\u2080  (dimensionless)', fontsize=9)
        ax.set_ylabel('Y (section coord)', fontsize=9)
        ax.grid(alpha=0.3)
        if plotted: ax.legend(fontsize=8)
    fig.suptitle('Normalized stress along section centreline (X=0)  —  '
                 'FEniCS vs Abaqus 3D, global axes  (\u03c3\u2080 = E\u00b7c/L)',
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(case_dir,'centerline_overlay.png'), dpi=130)
    plt.close(fig)

    # ── same overlay in the SECTION frame (added 7 Oct 2026) ────────────
    # Both curves pulled back with R (beam rotation at Z_MID): sigma_loc =
    # R^T sigma R. Beam theory: S11, S22, S12 ~ 0 here, S33 = bending,
    # S13/S23 = transverse shear. Shows FEniCS in-plane artifacts directly.
    if R is not None and fen_cl is not None and all(
            c in fen_cl.columns for c in STRESS_COMPS) and len(abq_cl) > 0:
        try:
            Rm = np.asarray(R, float)
            def _pull(df):
                T = np.array([[df['S11'].values, df['S12'].values, df['S13'].values],
                              [df['S12'].values, df['S22'].values, df['S23'].values],
                              [df['S13'].values, df['S23'].values, df['S33'].values]])
                Lc = np.einsum('ji,jkn,kl->iln', Rm, T, Rm)
                return {'S11': Lc[0,0], 'S22': Lc[1,1], 'S33': Lc[2,2],
                        'S12': Lc[0,1], 'S13': Lc[0,2], 'S23': Lc[1,2]}
            fl = _pull(fen_cl); al = _pull(abq_cl)
            fig, axes = plt.subplots(2, 3, figsize=(15, 9), squeeze=False)
            for ax, comp in zip(axes.flatten(), ['S11','S22','S33','S12','S13','S23']):
                ax.plot(fl[comp]/config.STRESS_SCALE, fen_cl['x2'].values, '-',
                        lw=1.8, color='#2471a3', label='FEniCS (warping)')
                ax.plot(al[comp]/config.STRESS_SCALE, abq_cl['Y'].values, '--',
                        lw=1.8, color='#c0392b', label='Abaqus 3D')
                ax.axvline(0, color='gray', lw=0.6, ls=':')
                ax.set_title('%s / \u03c3\u2080 (section frame)' % comp, fontsize=10)
                ax.set_xlabel('Stress / \u03c3\u2080', fontsize=9)
                ax.set_ylabel('Y (section coord)', fontsize=9)
                ax.grid(alpha=0.3); ax.legend(fontsize=8)
            fig.suptitle('Stress along section centreline (X=0) in the SECTION '
                         'frame (R\u1d40\u03c3R, beam rotation at Z_MID)', fontsize=11)
            fig.tight_layout()
            fig.savefig(os.path.join(case_dir, 'centerline_overlay_section.png'), dpi=130)
            plt.close(fig)
        except Exception as ex:
            tprint('  WARNING: section-frame overlay failed: %s' % ex)

    # ── strain overlay along Y at section centreline (X≈0) ──────────────
    # Same idea as the stress overlay above, but for strain: Abaqus's own
    # solid strain (LE33/LE13/LE23 -- already present in solid_mid.csv/
    # abq_df; principal_strain_from_solid_centerline() reads them at a
    # single centerline node, but here the WHOLE X≈0 column is plotted, not just
    # that one node) vs FEniCS's own strain field along the same line, IF
    # the FEniCS *_section_X0.csv export has it. The column-name
    # convention on the FEniCS side isn't guaranteed here (it depends on
    # what warping_core.postprocess() actually writes -- this file only
    # knows the STRESS names it uses, S33/S13/S23/Mises, from the
    # existing overlay above), so this tries a few plausible candidate
    # names and skips the FEniCS curve (with a one-time note, not a
    # crash) if none match -- the Abaqus curve is still plotted either
    # way. If postprocess() writes strain under some other name, add it
    # to fen_strain_col_candidates below.
    strain_cl_comps = ['LE11', 'LE22', 'LE33', 'LE12', 'LE13', 'LE23']   # all (7 Oct 2026)
    fen_strain_col_candidates = [
        tuple('LE' + s for s in ('11','22','33','12','13','23')),   # same convention as Abaqus
        tuple('E' + s for s in ('11','22','33','12','13','23')),    # warping_core names (Green-Lagrange)
        tuple('EE' + s for s in ('11','22','33','12','13','23')),
    ]
    fen_strain_cols = None
    if fen_cl is not None:
        for candidate in fen_strain_col_candidates:
            if all(c in fen_cl.columns for c in candidate):
                fen_strain_cols = dict(zip(strain_cl_comps, candidate))
                break
        if fen_strain_cols is None:
            tprint('  NOTE: no recognized strain columns in %s (tried %s) -- '
                   'strain_overlay.png will show Abaqus-only data. If '
                   'warping_core.postprocess() writes strain under a '
                   'different name, add it to fen_strain_col_candidates '
                   'in compare.py.'
                   % (cl_files[0] if cl_files else '<no file>',
                      fen_strain_col_candidates))

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), squeeze=False)
    axes = axes.flatten()
    for ax, comp in zip(axes, strain_cl_comps):
        plotted = False
        if fen_strain_cols is not None and 'x2' in fen_cl.columns:
            fen_col = fen_strain_cols[comp]
            ax.plot(fen_cl[fen_col].values, fen_cl['x2'].values,
                    '-', lw=1.8, color='#2471a3', label='FEniCS (warping)')
            plotted = True
        if comp in abq_cl.columns and len(abq_cl) > 0:
            # Abaqus shear strains (LE12/13/23) are ENGINEERING strains;
            # warping_core's E12/13/23 are tensor (Green-Lagrange) strains.
            # Halve the Abaqus shear so both curves are tensor components.
            _f = 0.5 if comp in ('LE12', 'LE13', 'LE23') else 1.0
            ax.plot(_f * abq_cl[comp].values, abq_cl['Y'].values,
                    '--', lw=1.8, color='#c0392b',
                    label='Abaqus 3D' + (' (x 1/2, tensor)' if _f != 1.0 else ''))
            plotted = True
        ax.axvline(0, color='gray', lw=0.6, ls=':')
        ax.set_title(comp + (' (tensor shear)' if comp in ('LE12', 'LE13', 'LE23') else ''), fontsize=10)
        ax.set_xlabel('Strain (dimensionless)', fontsize=9)
        ax.set_ylabel('Y (section coord)', fontsize=9)
        ax.grid(alpha=0.3)
        if plotted: ax.legend(fontsize=8)
    fig.suptitle('Strain along section centreline (X=0), global axes  —  '
                 'FEniCS (Green-Lagrange) vs Abaqus 3D (log strain)', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(case_dir,'strain_overlay.png'), dpi=130)
    plt.close(fig)

    if abq_n is not None and abq_m is not None:
        force_rows = [{
            'N1_FEniCS': n_out[0], 'N1_Abaqus': abq_n[0], 'N1_err': n_out[0]-abq_n[0],
            'N1_FEniCS_norm': n_out[0]/config.FORCE_SCALE, 'N1_Abaqus_norm': abq_n[0]/config.FORCE_SCALE,
            'N1_err_norm': (n_out[0]-abq_n[0])/config.FORCE_SCALE,
            'N2_FEniCS': n_out[1], 'N2_Abaqus': abq_n[1], 'N2_err': n_out[1]-abq_n[1],
            'N2_FEniCS_norm': n_out[1]/config.FORCE_SCALE, 'N2_Abaqus_norm': abq_n[1]/config.FORCE_SCALE,
            'N2_err_norm': (n_out[1]-abq_n[1])/config.FORCE_SCALE,
            'N3_FEniCS': n_out[2], 'N3_Abaqus': abq_n[2], 'N3_err': n_out[2]-abq_n[2],
            'N3_FEniCS_norm': n_out[2]/config.FORCE_SCALE, 'N3_Abaqus_norm': abq_n[2]/config.FORCE_SCALE,
            'N3_err_norm': (n_out[2]-abq_n[2])/config.FORCE_SCALE,
            'M1_FEniCS': m_out[0], 'M1_Abaqus': abq_m[0], 'M1_err': m_out[0]-abq_m[0],
            'M1_FEniCS_norm': m_out[0]/config.MOMENT_SCALE, 'M1_Abaqus_norm': abq_m[0]/config.MOMENT_SCALE,
            'M1_err_norm': (m_out[0]-abq_m[0])/config.MOMENT_SCALE,
            'M2_FEniCS': m_out[1], 'M2_Abaqus': abq_m[1], 'M2_err': m_out[1]-abq_m[1],
            'M2_FEniCS_norm': m_out[1]/config.MOMENT_SCALE, 'M2_Abaqus_norm': abq_m[1]/config.MOMENT_SCALE,
            'M2_err_norm': (m_out[1]-abq_m[1])/config.MOMENT_SCALE,
            'M3_FEniCS': m_out[2], 'M3_Abaqus': abq_m[2], 'M3_err': m_out[2]-abq_m[2],
            'M3_FEniCS_norm': m_out[2]/config.MOMENT_SCALE, 'M3_Abaqus_norm': abq_m[2]/config.MOMENT_SCALE,
        }]
        pd.DataFrame(force_rows).to_csv(os.path.join(case_dir,'force_error.csv'),
                                        index=False, float_format='%.8e')

"""
plotting.py
===========
The three end-of-batch summary plots (error histograms, energy scatter,
strain-screen audit). Deliberately the LIGHTEST module in this codebase
dependency-wise — matplotlib + pandas + numpy + config only, nothing
Abaqus/FEniCS/pyvista — so these can be regenerated from an existing
summary.csv without re-running anything else. See plot_results.py for
the standalone entry point that does exactly that.
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config import STRAIN_THRESHOLD, ENERGY_SCALE, tprint


def plot_histograms(df, out_dir):
    """Error plotted as rel-L2 PERCENTAGE — this is already dimensionless
    and self-normalizing (error divided by the reference value's own
    magnitude), so no separate physical scale factor is needed here.
    Physical-scale normalization is used instead where ACTUAL VALUES (not
    error ratios) are shown, e.g. plot_energy_scatter below — a raw value
    has no built-in reference to compare against, so it needs one.

    N_resultant / M_resultant are the TOTAL vector magnitude of the
    force/moment discrepancy (norm across all 3 components combined),
    not a per-component error.
    """
    titles = {
        'relL2_Mises':    'von Mises stress',
        'relL2_inv_1':    'I1  (hydrostatic stress, trace of sigma)',
        'relL2_J3':       'J3  (det. of deviatoric stress, 3rd invariant)',
        'relL2_S33':      'S33  (axial stress)',
        'relL2_S13':      'S13  (transverse shear stress)',
        'relL2_S23':      'S23  (transverse shear stress)',
        'relL2_N_resultant': 'Resultant force  (total vector magnitude)',
        'relL2_M_resultant': 'Resultant moment  (total vector magnitude)',
        'relL2_Energy':      'Strain energy per unit length',
    }
    cols = [c for c in df.columns if c.startswith('relL2_')]
    n = len(cols); ncols = 4; nrows = (n+ncols-1)//ncols
    fig,axes = plt.subplots(nrows,ncols,figsize=(5*ncols,4*nrows),squeeze=False)
    n_bins = 30   # plain equal-WIDTH bins on the error axis — halved bin
                  # WIDTH vs the previous n_bins=15 (bin width =
                  # range/n_bins, so halving the width means DOUBLING
                  # the count; if "half" was meant as half as many bins
                  # instead, i.e. coarser/wider bins, set this to 7 or 8).
    for ax,col in zip(axes.flatten(),cols):
        vals = df[col].dropna().values * 100
        if len(vals) < 2:
            ax.text(0.5,0.5,'not enough data',ha='center',va='center',
                    transform=ax.transAxes); continue
        # Standard histogram: bin edges are equal-width slices of the
        # error range (x-axis). Frequency (count per bin, y-axis) then
        # varies naturally bar-to-bar — that's what shows the distribution.
        ax.hist(vals, bins=n_bins, color='#2471a3', edgecolor='white',
               alpha=0.85)
        ax.axvline(np.median(vals),color='crimson',lw=1.5,
                   label='median %.2f%%'%np.median(vals))
        ax.set_title(titles.get(col,col),fontsize=9)
        ax.set_xlabel('rel-L2 error [%]')
        ax.set_ylabel('frequency')
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
    for ax in axes.flatten()[n:]: ax.axis('off')
    fig.suptitle('Error distribution — %d cases'%len(df),fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir,'histograms.png'),dpi=140)
    plt.close(fig); tprint('Saved histograms.png')


def plot_energy_scatter(df, out_dir):
    """This plots ACTUAL ENERGY VALUES (not an error ratio), so unlike
    the histograms above, there's no built-in reference to normalize
    against — a raw energy scatter plot would need arbitrary absolute
    units. Both axes are therefore divided by ENERGY_SCALE = U0 = E I / L^2
    (eq. 12; per-unit-length energy has force units, same
    non-dimensionalization convention as F0), giving a clearly-labeled
    NORMALIZED energy in both cases."""
    if not {'ELSE_per_length','Phi_FEniCS'}.issubset(df.columns): return
    v = df.dropna(subset=['ELSE_per_length','Phi_FEniCS'])
    if len(v)<2: return
    x_nd = v['ELSE_per_length'] / ENERGY_SCALE
    y_nd = v['Phi_FEniCS']      / ENERGY_SCALE
    fig,ax = plt.subplots(figsize=(6,5))
    ax.scatter(x_nd, y_nd, s=18, alpha=0.6, color='#117a3f')
    lim = [min(x_nd.min(),y_nd.min()), max(x_nd.max(),y_nd.max())]
    ax.plot(lim,lim,'k--',lw=0.8,label='1:1')
    ax.set_xlabel('Normalised energy, 3D Abaqus  =  (ELSE/length) / U0')
    ax.set_ylabel('Normalised energy, FEniCS  =  Phi / U0')
    ax.set_title('Strain energy at z=50: 3D vs FEniCS  (U0 = E I / L^2)')
    ax.legend(); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(os.path.join(out_dir,'energy_scatter.png'),dpi=140)
    plt.close(fig); tprint('Saved energy_scatter.png')


def plot_strain_screen_audit(df, out_dir):
    """POST-screen audit plot: analytical (pre-run, closed-form) root
    strain vs Abaqus-output (post-run, SE*/SK*-based) root strain, per
    case, against the 1:1 line and the STRAIN_THRESHOLD bound. Shows how
    conservative/accurate the pre-run screen actually is relative to real
    Abaqus output, which is exactly the question the pre-run screen can't
    answer about itself (it has to be evaluated before Abaqus runs)."""
    if not {'eps_root','eps_root_abaqus'}.issubset(df.columns):
        return
    v = df.dropna(subset=['eps_root','eps_root_abaqus'])
    if len(v) < 2:
        return
    fig, ax = plt.subplots(figsize=(6,5.5))
    over = v['eps_root_abaqus'] > STRAIN_THRESHOLD
    ax.scatter(v.loc[~over,'eps_root'], v.loc[~over,'eps_root_abaqus'],
               s=18, alpha=0.6, color='#117a3f', label='post-screen OK')
    if over.any():
        ax.scatter(v.loc[over,'eps_root'], v.loc[over,'eps_root_abaqus'],
                   s=24, alpha=0.85, color='#c0392b',
                   label='post-screen FAIL (Abaqus root strain > threshold)')
    lim = [0, max(v['eps_root'].max(), v['eps_root_abaqus'].max()) * 1.05]
    ax.plot(lim, lim, 'k--', lw=0.8, label='analytical = Abaqus (1:1)')
    ax.axhline(STRAIN_THRESHOLD, color='gray', lw=1.0, ls=':',
               label='STRAIN_THRESHOLD = %.1e' % STRAIN_THRESHOLD)
    ax.axvline(STRAIN_THRESHOLD, color='gray', lw=1.0, ls=':')
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel('eps_root  (analytical, pre-run screen)')
    ax.set_ylabel('eps_root_abaqus  (post-run, SE*/SK*-based)')
    ax.set_title('Pre-run screen vs Abaqus-output post-screen — %d cases' % len(v))
    ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir,'strain_screen_audit.png'), dpi=140)
    plt.close(fig); tprint('Saved strain_screen_audit.png')

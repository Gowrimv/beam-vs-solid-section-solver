#!/usr/bin/env python3
"""
plot_results.py
================
Standalone plot-only entry point — regenerates histograms.png,
energy_scatter.png, and strain_screen_audit.png from an EXISTING
summary.csv, without touching Abaqus, FEniCS, or pyvista at all.

This exists specifically so you can iterate on plot styling (bin count,
titles, colors, ...) against a batch that already ran, in seconds,
instead of re-running batch_driver.py (which needs Abaqus licenses and
a FEniCS solve per case) just to see an updated figure. It only imports
plotting.py, config.py and logutil.py — none of which touch Abaqus,
dolfin/FEniCS, or pyvista — plus pandas/matplotlib/numpy.

Usage:
    python3 plot_results.py
    python3 plot_results.py --out-dir /mnt/c/.../batch_test/cases
    python3 plot_results.py --summary-csv /path/to/summary.csv --plots-dir /path/to/write/plots/to
"""
import argparse
import os

import pandas as pd

from config import tprint
from plotting import plot_histograms, plot_energy_scatter, plot_strain_screen_audit


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out-dir', type=str, required=True,
                    help='Directory containing summary.csv (default: '
                         '%(default)s — the same default batch_driver.py '
                         'uses). Plots are written back into this same '
                         'directory unless --plots-dir is given.')
    ap.add_argument('--summary-csv', type=str, default=None,
                    help='Explicit path to summary.csv, if it is not '
                         'directly at <out-dir>/summary.csv.')
    ap.add_argument('--plots-dir', type=str, default=None,
                    help='Where to write histograms.png/energy_scatter.png/'
                         'strain_screen_audit.png. Defaults to --out-dir.')
    args = ap.parse_args()

    csv_path = args.summary_csv or os.path.join(args.out_dir, 'summary.csv')
    if not os.path.isfile(csv_path):
        tprint('No summary.csv found at %s — run batch_driver.py first '
               '(or pass --summary-csv / --out-dir).' % csv_path)
        return

    df = pd.read_csv(csv_path)
    plots_dir = args.plots_dir or args.out_dir
    os.makedirs(plots_dir, exist_ok=True)

    tprint('Loaded %d cases from %s' % (len(df), csv_path))
    plot_histograms(df, plots_dir)
    plot_energy_scatter(df, plots_dir)
    plot_strain_screen_audit(df, plots_dir)
    tprint('Done. Plots in: ' + plots_dir)


if __name__ == '__main__':
    main()

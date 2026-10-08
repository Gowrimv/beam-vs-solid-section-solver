"""
loads_io.py
===========
Export/replay of accepted (Vx, Vy, Mx) load triples between batch runs.

Use case: run the batch once on one SECTION_SIDE, export the loads that
were actually ACCEPTED (passed strain_screen / bending_dominance / the
post-Abaqus admission gates), then feed that exact same set of loads
into a run on a different SECTION_SIDE with --loads-csv, instead of
drawing a fresh independent random sample. Because a smaller section
sees higher strain for the same load, loads accepted on the SMALLER
section will almost always also pass the gates on a LARGER section --
so export from the smaller run and replay onto the larger one.

Depends on nothing but pandas -- no Abaqus, FEniCS, or pyvista -- so it
can be imported anywhere, same as strain_checks.py / plotting.py.
"""
import pandas as pd

from config import tprint


def export_accepted_loads(summary_csv, out_csv):
    """
    Post-run: pull the accepted (Vx, Vy, Mx) triples out of a completed
    batch's summary.csv and write them to a standalone loads file that
    a later run (e.g. on a different SECTION_SIDE) can replay via
    load_fixed_loads() / --loads-csv.

    summary.csv already only contains admitted cases, but if the
    'admit_all' column is present it is filtered on explicitly, in
    case that ever changes.
    """
    df = pd.read_csv(summary_csv)
    if 'admit_all' in df.columns:
        df = df[df['admit_all'] == True]  # noqa: E712

    out = df[['case', 'Vx', 'Vy', 'Mx']].reset_index(drop=True)
    out.to_csv(out_csv, index=False, float_format='%.10e')
    tprint('export_accepted_loads: wrote %d accepted loads to %s'
           % (len(out), out_csv))
    return out


def load_fixed_loads(loads_csv):
    """
    Load a previously-exported loads file (from export_accepted_loads)
    for replay on a new run. Returns a DataFrame with Vx, Vy, Mx columns
    in the same shape main()'s rng.uniform(...) draw produces, so it can
    be dropped straight in place of it.
    """
    df = pd.read_csv(loads_csv)
    missing = {'Vx', 'Vy', 'Mx'} - set(df.columns)
    if missing:
        raise ValueError('loads_csv %s is missing columns: %s'
                          % (loads_csv, missing))
    draws = df[['Vx', 'Vy', 'Mx']].reset_index(drop=True)
    tprint('load_fixed_loads: replaying %d fixed loads from %s'
           % (len(draws), loads_csv))
    return draws

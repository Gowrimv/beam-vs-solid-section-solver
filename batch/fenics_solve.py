"""
fenics_solve.py
================
The in-process FEniCS warping solve: builds the mesh/problem ONCE
(init_fenics(), called once per batch from main()) and re-solves it per
case (solve_and_postprocess()) under a module-level lock, since
ThreadPoolExecutor workers share this process's memory and
forward_solve() mutates the shared problem dict P["u"] in place. This
is the one module in the split that requires FEniCS/dolfin to import --
keeping it separate means strain_checks.py/abaqus_io.py/plotting.py
stay importable (e.g. for plot_results.py-style offline work) on a
machine without dolfin installed.

NOTE: solve_and_postprocess() would need per-worker copies of P instead
of a shared lock if this were ever switched from ThreadPoolExecutor to
ProcessPoolExecutor.
"""
import sys
import threading

import numpy as np

import config
from config import tprint

# ── FEniCS imports (WSL-native) ────────────────────────────────────────────
sys.path.insert(0, '/mnt/c/Users/macLab/Desktop/FEniCSfiles')
from warping_core import (load_mesh, build_problem, forward_solve,
                          postprocess, rotation_matrix_from_ur, assign_targets)  # noqa: F401 -- rotation_matrix_from_ur is re-exported for batch_driver.py
from dolfin import assemble, dx, tr


_fenics_lock = threading.Lock()
_mesh = None
_P    = None


def fenics_ready():
    """True once init_fenics() has built the shared problem. batch_driver.py
    checks this (rather than importing _P by name, which would freeze a
    stale None reference at import time) before attempting the fast
    FEniCS-in-process path for a case."""
    return _P is not None

def init_fenics(mesh_xmf, workdir):
    global _mesh, _P
    tprint("Loading mesh and building FEniCS problem (once)...")
    _mesh, _ = load_mesh(mesh_xmf, workdir=workdir)
    _P = build_problem(_mesh, E=config.FENICS_E, nu_val=config.NU, L_scale=1.0)
    tprint('FEniCS solved with E = %g; outputs rescaled by E/FENICS_E = %g'
           % (config.FENICS_E, config.E / config.FENICS_E))
    tprint("FEniCS ready. Mesh cells: %d" % _mesh.num_cells())
    # ── coordinate diagnostic (helps detect centroid/scale issues) ──
    coords = _mesh.coordinates()
    tprint('Mesh X: [%.6f, %.6f]  Y: [%.6f, %.6f]  centroid: (%.6f, %.6f)' % (
        coords[:,0].min(), coords[:,0].max(),
        coords[:,1].min(), coords[:,1].max(),
        coords[:,0].mean(), coords[:,1].mean()))

def solve_and_postprocess(n_t_sp, m_t_sp, sv_mat, sk_mat, case_dir, base,
                           R_local, R_out):
    """assign_targets -> forward_solve -> strain-energy assembly ->
    postprocess, all under _fenics_lock, using the module-level problem
    built by init_fenics(). Returns (Phi, n_out, m_out) — strain energy
    and the global-frame force/moment resultants, matching
    run_abaqus_comparison.py's convention exactly."""
    if _P is None:
        raise RuntimeError('FEniCS problem not built — call init_fenics() first')
    import time as _t
    t0 = _t.time()
    with _fenics_lock:
        t1 = _t.time()
        # Solve at E = FENICS_E: forces/moments scale with E, strains don't.
        k = config.E / config.FENICS_E
        assign_targets(_P, np.asarray(n_t_sp, float) / k, np.asarray(m_t_sp, float) / k)
        forward_solve(_P, sv_mat, sk_mat, warm=False, n_steps=10)
        t2 = _t.time()
        Phi = float(assemble(
            ((_P['lmbda_param']/2)*tr(_P['E_strain'])**2
             + _P['mu_param']*tr(_P['E_strain']*_P['E_strain']))*dx)
            * _P['L_scale']**2)   # matches run_abaqus_comparison.py exactly
        result_pp = postprocess(_P, case_dir, base,
                                R_local=R_local, R_out=R_out,
                                do_plots=getattr(config, 'FENICS_PLOTS', True),
                                do_vtu=True)
        t3 = _t.time()
    # ── back to the real E (stress units scale by k) ──────────────────────
    _rescale_stress_outputs(case_dir, base, k)
    tprint('[%s] FEniCS timing: waited for lock %.0fs, solve %.0fs, '
           'postprocess %.0fs' % (base, t1-t0, t2-t1, t3-t2))
    n_out = k * np.asarray(result_pp['n0_out'], float)
    m_out = k * np.asarray(result_pp['m0_out'], float)
    return k * Phi, n_out, m_out


# Stress-unit fields written by warping_core.postprocess (strains E11..E33
# and displacement U are independent of E and are left alone).
_STRESS_FIELDS = ('S11', 'S22', 'S33', 'S12', 'S13', 'S23', 'Mises', 'inv_1')

def _rescale_stress_outputs(case_dir, base, k):
    """Multiply the stress-unit fields of <base>_ref.vtu, <base>_deformed.vtu
    and <base>_section_X0.csv by k = E / FENICS_E, in place, right after
    postprocess() has (re)written them. No-op if k == 1."""
    if abs(k - 1.0) < 1e-12:
        return
    import os
    import pandas as pd
    for suffix in ('_ref.vtu', '_deformed.vtu'):
        p = os.path.join(case_dir, base + suffix)
        if not os.path.isfile(p):
            continue
        try:
            import pyvista as pv
            mesh = pv.read(p)
            for data in (mesh.point_data, mesh.cell_data):
                for name in _STRESS_FIELDS:
                    if name in data.keys():
                        data[name] = data[name] * k
            mesh.field_data['stress_scale_applied'] = [k]
            mesh.save(p)
        except Exception as ex:
            tprint('[%s] WARNING: could not rescale %s: %s' % (base, p, ex))
    p = os.path.join(case_dir, base + '_section_X0.csv')
    if os.path.isfile(p):
        try:
            df = pd.read_csv(p)
            for name in _STRESS_FIELDS:
                if name in df.columns:
                    df[name] = df[name] * k
            df.to_csv(p, index=False, float_format='%.6e')
        except Exception as ex:
            tprint('[%s] WARNING: could not rescale %s: %s' % (base, p, ex))

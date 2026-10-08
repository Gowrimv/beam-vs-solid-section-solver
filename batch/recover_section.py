"""
recover_section.py
==================
Reads beam_mid.csv (written by extract_odb.py), maps SE/SK/UR -> FEniCS
material frame, runs warping forward solve, postprocesses stress field.

Usage:
    python3 recover_section.py <case_dir> <mesh.xmf> <beam_mid.csv> [Z_mid] [E] [nu]

Outputs in <case_dir>/:
    *_ref.vtu          FEniCS nodal stress field (for l2 comparison)
    *_section_X0.csv   centerline stress profile
    *_resultants.csv   n0_out / m0_out in Abaqus frame
"""
import sys, os
import numpy as np
import pandas as pd

# Ensure warping_core.py is findable.
# By default it is expected alongside this script.
# Override with:  WARPING_CORE_DIR=/your/path  wsl python3 recover_section.py ...
_here   = os.path.dirname(os.path.abspath(__file__))
_wc_dir = os.environ.get('WARPING_CORE_DIR', _here)
for _p in [_here, _wc_dir]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

from warping_core import (load_mesh, build_problem, forward_solve,
                          postprocess, rotation_matrix_from_ur, assign_targets)

case_dir = sys.argv[1]
msh_file = sys.argv[2]
csv_file = sys.argv[3]
Z_MID    = float(sys.argv[4]) if len(sys.argv) > 4 else 50.0
E_MOD    = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
NU       = float(sys.argv[6]) if len(sys.argv) > 6 else 0.3

os.makedirs(case_dir, exist_ok=True)

# ── load beam CSV, snap to Z_MID ─────────────────────────────────────────────
df = pd.read_csv(csv_file)
df.columns = df.columns.str.strip()
idx = (df['Z'] - Z_MID).abs().idxmin()
row = df.iloc[idx]
print("Snapped to Z=%.4f (row %d)" % (float(row['Z']), idx))

# ── SE/SK -> FEniCS spatial frame (beam axis = Z, section axes = X,Y) ────────
# Abaqus B31 with beam along Z: SE1=axial(Z), SE2=shear(Y), SE3=shear(X)
#                                SK1=bend(X),  SK2=bend(Y),  SK3=twist(Z)
# FEniCS frame: sv[0]=v1(x), sv[1]=v2(y), sv[2]=v3(axial)
sv_sp = np.array([float(row['SE3']), float(row['SE2']), float(row['SE1'])])
sk_sp = np.array([float(row['SK1']), float(row['SK2']), float(row['SK3'])])

ur1, ur2, ur3 = float(row['UR1']), float(row['UR2']), float(row['UR3'])
R_local = rotation_matrix_from_ur(ur1, ur2, ur3)
# Abaqus B31 SE/SK are already components in the beam-local (director)
# basis d_i = R e_i, i.e. exactly v0 = R^T v, k0 = R^T k of Arora et al.
# (2019). Pulling back again with R_local.T was a double rotation
# (fixed 2026-10-08, see PARAMETERS_AND_CHANGES.md 2.19). R_local is
# still used on the OUTPUT side (postprocess, abq_n/abq_m).
sv_mat  = sv_sp
sk_mat  = sk_sp

# Force/moment targets: raw beam-local SF/SM, same ordering as the batch
# (fenics_solve.solve_and_postprocess). build_problem() ENFORCES
# n0 = n0_t and m0 = m0_t; without assign_targets() both default to ZERO.
n_t_sp  = np.array([float(row['SF3']), float(row['SF2']), float(row['SF1'])])
m_t_sp  = np.array([float(row['SM1']), float(row['SM2']), float(row['SM3'])])

print("sv_mat=%s  sk_mat=%s" % (sv_mat, sk_mat))

# ── build FEniCS problem ──────────────────────────────────────────────────────
mesh, _ = load_mesh(msh_file, workdir=case_dir)
P       = build_problem(mesh, E=E_MOD, nu_val=NU, L_scale=1.0)
assign_targets(P, n_t_sp, m_t_sp)
base    = os.path.basename(case_dir)

# ── forward solve ─────────────────────────────────────────────────────────────
f_vec = forward_solve(P, sv_mat, sk_mat, warm=False, n_steps=10)
n0, m0 = f_vec[:3], f_vec[3:]

# ── postprocess -> _ref.vtu + centerline CSV ──────────────────────────────────
# R_MAP: FEniCS material (x1,x2 in section, x3=axial) -> Abaqus (X,Y,Z)
R_MAP  = np.eye(3)   # section axes already aligned; axial=Z in both
result = postprocess(P, case_dir, base,
                     R_local=R_local, R_out=R_MAP,
                     do_plots=True, do_vtu=True)

n_out = result['n0_out']
m_out = result['m0_out']
print("Resultants (Abaqus frame): N=%s  M=%s" % (n_out, m_out))

pd.DataFrame([{
    'Z_mid': Z_MID, 'E': E_MOD, 'nu': NU,
    'SE1': row['SE1'], 'SE2': row['SE2'], 'SE3': row['SE3'],
    'SK1': row['SK1'], 'SK2': row['SK2'], 'SK3': row['SK3'],
    'UR1': ur1, 'UR2': ur2, 'UR3': ur3,
    'n1': n_out[0], 'n2': n_out[1], 'n3': n_out[2],
    'm1': m_out[0], 'm2': m_out[1], 'm3': m_out[2],
}]).to_csv(os.path.join(case_dir, base+'_resultants.csv'),
           index=False, float_format='%.8e')
print("Done:", case_dir)
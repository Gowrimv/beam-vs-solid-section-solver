"""
warping_core.py
===============
HCB cross-sectional warping solver (Arora, Kumar, Steinmann 2019),
restructured function-by-function so the same forward solve serves:
  - single forward runs   (forward_main.py)
  - inverse Newton solves (inverse_main.py)

Design rules:
  * build_problem() runs ONCE. Strains enter through v0_c/k0_c Constants
    updated with .assign() -> no form recompilation between solves.
  * u persists inside the problem handle P -> warm starting is automatic.
    Only solve_with_load_stepping() resets u.
  * The SOLVE works purely in the FEniCS local material frame. R_local
    (deformed/spatial rotation) is applied only in postprocess(), never
    fed into build_problem() or the SNES residual/Jacobian.

Units convention (matches your current script with E_scale = 1):
  * E_mod is the PHYSICAL modulus (e.g. 2.1e5 MPa) set in build_problem.
  * Mesh coordinates are non-dimensional; L_scale converts to physical.
  * n0 carries L_scale**2, m0 carries L_scale**3 hidden as
    (x_hat_phys = x*L_scale) * dx * L_scale**2.

FRAME NOTE (read this before touching sv/sk/target in a driver script):
  v0_c/k0_c are the rod's LOCAL, frame-invariant strain measures. They
  are fed straight into F via R_id = Identity(3) in build_problem(), and
  build_problem() takes no rotation argument at all -- nothing in the
  SOLVE (F_rod, E_strain, the SNES residual/Jacobian) is ever rotated.
  compute_n0_m0() likewise always returns the FEniCS material-frame
  resultants -- this is the frame forward_solve/inverse_solve's Newton
  loop actually iterates in.

  postprocess()'s R_local and R_out are the ONE place in this file that
  rotate anything, and they are applied strictly on the way OUT, after
  the physics is done, to a plain numpy/UFL copy of the already-computed
  quantities -- never fed back into the solve. R_local is applied FIRST
  (this station's deformed/spatial rotation), R_out SECOND (e.g. R_MAP,
  axis convention), composed as R_total = R_out @ R_local:
    - n0, m0             one-sided  (R_total @ v)        -- vectors
    - sigma, E            two-sided  (R_total @ T @ R_total.T) -- tensors
    - displacement (VTU)  one-sided  (R_total @ u)
  result["n0_out"]/result["m0_out"] are the FULLY rotated resultants --
  use these directly for Abaqus comparison. Do NOT also rotate
  result["n0"]/result["m0"] (material frame) by R_MAP in the driver
  script -- that duplicates what postprocess() already did and, if
  R_local was also passed, is simply wrong (missing the R_local factor
  postprocess() already applied). Pass both R_local and R_out=None for
  pure spatial/deformed frame with no axis relabeling.

  inverse_solve()'s `target` argument must ALREADY be in that same
  local frame -- inverse_solve does not rotate it. If the
  target is in FEniCS axis order already ( the Abaqus->FEniCS
  R_MAP relabeling by hand), rotate it into local frame in the
  DRIVER SCRIPT, before calling inverse_solve, using
  target_to_material() below -- it removes ONLY R_local (UR), not
  R_MAP. Do not add rotation logic inside inverse_solve or forward_solve
  themselves -- that would duplicate what postprocess()'s R_out already
  does correctly on the way out, and the two are not guaranteed to
  cancel if either changes later.
"""

import os
import numpy as np
import pandas as pd
import meshio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dolfin import *
from dolfin import PETScOptions
from fenics import *

set_log_level(LogLevel.WARNING)
from collections import defaultdict

# ============================================================================
# MESH
# ============================================================================
def load_mesh(mesh_file, workdir="."):
    """meshio read -> xdmf -> dolfin Mesh, for ANY format meshio can read.

    meshio.read() dispatches on extension, so .vtk / .vtu / .xdmf / .msh /
    .med / .inp all work. This function no longer assumes Gmsh:
      - physical/cell tags are optional (default 0 if the file has none)
      - boundary (line) cells are optional -- if the file has none (e.g. a
        ParaView slice), the exterior edges are computed from the surface
        mesh and marked 0, so the returned facet XDMF is always valid.

    Handles triangle AND quad surface meshes. Returns (mesh, facet_xdmf_path).
    """
    m = meshio.read(mesh_file)

    surface_types  = ("triangle", "quad")
    boundary_types = ("line", "line3")

    found_surface  = next((c.type for c in m.cells if c.type in surface_types), None)
    found_boundary = next((c.type for c in m.cells if c.type in boundary_types), None)

    if found_surface is None:
        raise ValueError(
            f"No triangle or quad cells found in {mesh_file}. "
            f"Cell types present: {[c.type for c in m.cells]}"
        )

    print(f"Mesh cell type: {found_surface}  "
          f"boundary type: {found_boundary or '(none in file -> will synthesize)'}")

    # --- tags: look for a physical-id array under several common keys ---------
    def get_tags(cell_type):
        n = len(m.get_cells_type(cell_type))
        for key in ("gmsh:physical", "medit:ref", "cell_tags", "CellEntityIds"):
            if key in m.cell_data_dict and cell_type in m.cell_data_dict[key]:
                return np.asarray(m.cell_data_dict[key][cell_type], dtype=int)
        return np.zeros(n, dtype=int)          # no tags in this format -> all 0

    def make_mesh(cell_type, cells, tags, prune_z=True):
        pts = m.points[:, :2] if prune_z else m.points
        return meshio.Mesh(points=pts,
                           cells={cell_type: cells},
                           cell_data={"name_to_read": [tags]})

    # --- boundary: use file's line cells, or compute exterior edges ----------
    if found_boundary is not None:
        bnd_cells = m.get_cells_type(found_boundary)
        bnd_tags  = get_tags(found_boundary)
        bnd_type  = found_boundary
    else:
        surf = m.get_cells_type(found_surface)
        edge_count  = defaultdict(int)
        edge_sample = {}
        for cell in surf:                       # works for tri (3) and quad (4)
            k = len(cell)
            for i in range(k):
                a, b = int(cell[i]), int(cell[(i + 1) % k])
                key = (a, b) if a < b else (b, a)
                edge_count[key] += 1
                edge_sample.setdefault(key, (a, b))
        bnd_cells = np.array([edge_sample[k] for k, c in edge_count.items() if c == 1],
                             dtype=int)         # edges owned by exactly one cell
        bnd_tags  = np.zeros(len(bnd_cells), dtype=int)
        bnd_type  = "line"
        print(f"  synthesized {len(bnd_cells)} boundary edges (all marked 0)")

    # --- write the two XDMF files dolfin needs -------------------------------
    base       = os.path.splitext(os.path.basename(mesh_file))[0]
    mesh_xdmf  = os.path.join(workdir, f"{base}_mesh.xdmf")
    facet_xdmf = os.path.join(workdir, f"{base}_facet_mesh.xdmf")

    meshio.write(mesh_xdmf,
                 make_mesh(found_surface, m.get_cells_type(found_surface),
                           get_tags(found_surface)))
    meshio.write(facet_xdmf,
                 make_mesh(bnd_type, bnd_cells, bnd_tags))

    mesh = Mesh()
    with XDMFFile(mesh_xdmf) as f:
        f.read(mesh)

    print(f"Cells: {mesh.num_cells()}, Vertices: {mesh.num_vertices()}, "
          f"Cell type: {found_surface}")
    return mesh, facet_xdmf


# ============================================================================
# MESH
# ============================================================================
# def load_mesh(msh_file, workdir="."):
#     """meshio read -> xdmf -> dolfin Mesh.

#     Works for triangle AND quad surface meshes -- cell type is auto-detected
#     from the .msh file rather than hardcoded.

#     The two intermediate .xdmf files (mesh.xdmf, facet_mesh.xdmf) are a
#     FEniCS 2019 (legacy dolfin) requirement: dolfin cannot read .msh
#     directly and must load via XDMF. They are written to workdir and can
#     be ignored after the run.

#     Returns (mesh, facet_xdmf_path).
#     """
#     mesh_from_file = meshio.read(msh_file)

#     # auto-detect 2D surface cell type (triangle or quad)
#     surface_types = {"triangle", "quad"}
#     boundary_types = {"line", "line3"}
#     found_surface = None
#     found_boundary = None
#     for cell_block in mesh_from_file.cells:
#         if cell_block.type in surface_types and found_surface is None:
#             found_surface = cell_block.type
#         if cell_block.type in boundary_types and found_boundary is None:
#             found_boundary = cell_block.type

#     if found_surface is None:
#         raise ValueError(
#             f"No triangle or quad cells found in {msh_file}. "
#             f"Cell types present: {[c.type for c in mesh_from_file.cells]}"
#         )
#     if found_boundary is None:
#         raise ValueError(
#             f"No line/boundary cells found in {msh_file}. "
#             f"Cell types present: {[c.type for c in mesh_from_file.cells]}"
#         )

#     print(f"Mesh cell type: {found_surface}  boundary type: {found_boundary}")

#     def create_mesh(m, cell_type, prune_z=False):
#         cells     = m.get_cells_type(cell_type)
#         cell_data = m.get_cell_data("gmsh:physical", cell_type)
#         points    = m.points[:, :2] if prune_z else m.points
#         return meshio.Mesh(points=points, cells={cell_type: cells},
#                            cell_data={"name_to_read": [cell_data]})

#     msh_basename = os.path.splitext(os.path.basename(msh_file))[0]
#     facet_xdmf = os.path.join(workdir, f"{msh_basename}_facet_mesh.xdmf")
#     mesh_xdmf  = os.path.join(workdir, f"{msh_basename}_mesh.xdmf")
#     meshio.write(facet_xdmf, create_mesh(mesh_from_file, found_boundary,
#                                           prune_z=True))
#     meshio.write(mesh_xdmf,  create_mesh(mesh_from_file, found_surface,
#                                           prune_z=True))

#     mesh = Mesh()
#     with XDMFFile(mesh_xdmf) as f:
#         f.read(mesh)

#     print(f"Cells: {mesh.num_cells()}, Vertices: {mesh.num_vertices()}, "
#           f"Cell type: {found_surface}")
#     return mesh, facet_xdmf


# ============================================================================
# KINEMATICS / MATERIAL (module-level helpers, same as your script)
# ============================================================================
def deformation_gradient_rod(v0, k0, x_field, R):
    e1 = Constant((1, 0, 0))
    e2 = Constant((0, 1, 0))
    e3 = Constant((0, 0, 1))
    term1 = outer(v0, e3)
    term2 = outer(cross(k0, x_field), e3)
    grad_x = grad(x_field)
    term3 = outer(grad_x[:, 0], e1) + outer(grad_x[:, 1], e2)
    return R * (term1 + term2 + term3)


def strain(F):
    return 0.5 * (F.T * F - Identity(3))


class InitialGuess(UserExpression):
    def eval(self, values, x):
        values[0] = x[0]
        values[1] = x[1]
        values[2] = 0
        values[3:15] = 0.0
    def value_shape(self):
        return (15,)
    
def assign_targets(P, n0_t, m0_t):
    L = P["L_scale"]
    P["n0_t"].assign(Constant(tuple(float(c) for c in n0_t)))
    P["m0_t"].assign(Constant(tuple(float(c) for c in m0_t)))

# ============================================================================
# BUILD PROBLEM  (runs ONCE)
# ============================================================================
def build_problem(mesh, E=2.1e5, nu_val=0.3, L_scale=1.0,
                  snes_tol=1e-8, quad_degree=8,
                  force_constraints=True, rotation_constraint="product",
                  constraint_mask=None):
    # constraint_mask (6 ints, order n1 n2 n3 m1 m2 m3): which resultant
    # constraints to keep when force_constraints=True (default: all 6).
    # TRIAL COPY (7 Oct 2026) of warping_core.build_problem with two switches:
    #   force_constraints=False -> the alpha/beta resultant constraints are
    #       removed (alpha, beta pinned to 0), so the section is driven by
    #       the prescribed strains (sv, sk) ONLY.
    #   rotation_constraint="rotation" -> third m_con component is the
    #       in-plane rigid-rotation constraint (X1-C1) u2 - (X2-C2) u1
    #       instead of u1*u2.
    """Build function space, Constants, variational form, solver, and all
    post-processing UFL forms ONCE. Returns handle dict P."""

    parameters["form_compiler"]["quadrature_degree"] = quad_degree

    # -- material (physical units, like your E_scale=1 version) ------------
    E_mod = Constant(E)
    nu = Constant(nu_val)
    rho0 = Constant(1)
    mu_param = E_mod / (2 * (1 + nu))
    lmbda_param = E_mod * nu / ((1 + nu) * (1 - 2 * nu))


    # -- function space (12 Reals: lambda, mu, alpha, beta) ----------------
    P2_vec = VectorElement("Lagrange", mesh.ufl_cell(), 2, dim=3)
    R_elem = FiniteElement("Real", mesh.ufl_cell(), 0)
    mixed_element = MixedElement([P2_vec] + [R_elem] * 12)
    V = FunctionSpace(mesh, mixed_element)

    u = Function(V)
    v_test = TestFunction(V)
    u_parts = split(u)
    u_sol    = u_parts[0]
    lambda_u = as_vector([u_parts[1], u_parts[2], u_parts[3]])
    mu_u     = as_vector([u_parts[4], u_parts[5], u_parts[6]])
    alpha    = as_vector([u_parts[7], u_parts[8], u_parts[9]])
    beta     = as_vector([u_parts[10], u_parts[11], u_parts[12]])

  

    

    print(f"Total DOFs: {V.dim()}")

 

    # -- geometry -----------------------------------------------------------
    V0 = FunctionSpace(mesh, "DG", 0)
    one_fn = Function(V0)
    one_fn.vector()[:] = 1.0
    total_area = assemble(one_fn * dx)
    x_coord = SpatialCoordinate(mesh)
    cx = assemble(x_coord[0] * dx) / total_area
    cy = assemble(x_coord[1] * dx) / total_area
    print(f"Area: {total_area:.8e}, centroid: ({cx:.8e}, {cy:.8e})")
    C = Constant((cx, cy, 0.0))

    # -- strain Constants: THE ONLY THING THAT CHANGES BETWEEN SOLVES -------
    v0_c = Constant((0.0, 0.0, 1.0))
    k0_c = Constant((0.0, 0.0, 0.0))

    # -- variational form (built once, compiled once) -----------------------
    # R_id = Identity(3): the solve ALWAYS happens in the local/material
    # frame. R_local is never used here -- see FRAME NOTE at top of file.
    R_id = Identity(3)
    X = SpatialCoordinate(mesh)
    u_rel = u_sol - C

    F_rod = deformation_gradient_rod(v0_c, k0_c, u_rel, R_id)
    E_strain = strain(F_rod)
    W = (lmbda_param / 2) * tr(E_strain) ** 2 + mu_param * tr(E_strain * E_strain)
    # constraint vector m 
    # m_con = as_vector([u_rel[1] * u_rel[2],
    #                    u_rel[2] * u_rel[0],
    #                    (X[0] - C[0]) * u_rel[1] - (X[1] - C[1]) * u_rel[0]])
    if rotation_constraint == "rotation":
        m_con = as_vector([u_rel[1] * u_rel[2],
                           u_rel[2] * u_rel[0],
                           (X[0] - C[0]) * u_rel[1] - (X[1] - C[1]) * u_rel[0]])
    else:
        m_con = as_vector([u_rel[1] * u_rel[2],
                           u_rel[2] * u_rel[0],
                           u_rel[0] * u_rel[1]])
    

    e3_const = Constant((0, 0, 1))
    S_pp = lmbda_param * tr(E_strain) * Identity(3) + 2 * mu_param * E_strain
    P_pp = F_rod * S_pp
    Pe3_pp = P_pp * e3_const
    x_hat_phys = u_rel * L_scale
    xPe3_pp = cross(x_hat_phys, Pe3_pp)

    n0_t, m0_t = Constant((0.,0.,0.)), Constant((0.,0.,0.))
    A = Constant(total_area)
    n_ref = Constant(E * total_area * L_scale**2)          # scaling
    m_ref = Constant(E * total_area * L_scale**3)

    n0_form = Pe3_pp * L_scale**2
    m0_form = xPe3_pp * L_scale**2

    Pi_total = (W + rho0 * (inner(lambda_u, u_rel) + inner(mu_u, m_con))) * dx
    if force_constraints and constraint_mask is not None:
        _mk = [int(v) for v in constraint_mask]
        _terms = 0
        for _i in range(3):
            _terms = _terms + (alpha[_i] * (n0_form[_i] - n0_t[_i]/A) / n_ref if _mk[_i]
                               else 0.5 * alpha[_i]**2)
            _terms = _terms + (beta[_i] * (m0_form[_i] - m0_t[_i]/A) / m_ref if _mk[3+_i]
                               else 0.5 * beta[_i]**2)
        Pi_total += _terms * dx
    elif force_constraints:
        Pi_total += (inner(alpha, (n0_form - n0_t/A) / n_ref)
                + inner(beta,  (m0_form - m0_t/A) / m_ref)) * dx
    else:
        # keep the alpha/beta DOFs but pin them to 0 (otherwise their rows
        # in the Jacobian are empty and LU fails)
        Pi_total += 0.5 * (inner(alpha, alpha) + inner(beta, beta)) * dx
   



    F_residual = derivative(Pi_total, u, v_test)
    J_jacobian = derivative(F_residual, u)

    problem = NonlinearVariationalProblem(F_residual, u, [], J_jacobian)

    PETScOptions.clear()
    PETScOptions.set("snes_type", "newtonls")
    PETScOptions.set("snes_linesearch_type", "bt")
    PETScOptions.set("ksp_type", "preonly")
    PETScOptions.set("pc_type", "lu")
    PETScOptions.set("snes_max_it", 200)
    PETScOptions.set("snes_rtol", snes_tol)
    PETScOptions.set("snes_atol", snes_tol)
    PETScOptions.set("snes_stol", snes_tol)

    solver = NonlinearVariationalSolver(problem)
    solver.parameters["nonlinear_solver"] = "snes"
    solver.parameters["snes_solver"]["error_on_nonconvergence"] = True
    solver.parameters["snes_solver"]["report"] = False

    # -- post-processing forms (built once, against the SAME F_rod/E_strain
    #    used in the solve above -- no second deformation-gradient build) ---
    # S_pp = lmbda_param * tr(E_strain) * Identity(3) + 2 * mu_param * E_strain  # 2nd PK
    # P_pp = F_rod * S_pp                     # 1st PK
    # e3_const = Constant((0, 0, 1))
    # Pe3_pp = P_pp * e3_const                # traction on e3 face
    

    J_pp = det(F_rod)
    sigma_pp = (1 / J_pp) * F_rod * S_pp * F_rod.T   # Cauchy stress

    # -- output-rotation Constant: ONE object, reused across every
    #    postprocess() call via .assign(). A fresh Constant(...) built
    #    inside postprocess() each time would be a new UFL terminal ->
    #    FFC cache miss -> recompiled projection kernels every call.
    #    Only the Constant itself needs to persist here -- the
    #    sigma_out/E_out expressions built FROM it can be (and are)
    #    rebuilt fresh in postprocess() every call with no cost, since
    #    UFL's form signature is keyed on terminal identity, not on the
    #    Python object identity of composite expressions built from them.
    #    Default value is Identity -- unrotated (FEniCS frame) output.
    #    Holds the COMPOSED rotation (R_out @ R_local) once postprocess()
    #    assigns into it -- see postprocess() docstring. ---
    R_out_c = Constant(np.eye(3))

    P = {
        "mesh": mesh, "V": V, "u": u,
        "v0_c": v0_c, "k0_c": k0_c,
        "solver": solver,
        "C": C, "cx": cx, "cy": cy, "total_area": total_area,
        "L_scale": L_scale,
        "E_mod": E_mod, "nu": nu,
        "mu_param": mu_param, "lmbda_param": lmbda_param,
        # UFL handles for post-processing (all track v0_c/k0_c/u automatically)
        "F_rod": F_rod, "E_strain": E_strain,
        "S_pp": S_pp, "P_pp": P_pp, "Pe3_pp": Pe3_pp, "xPe3_pp": xPe3_pp,
        "J_pp": J_pp, "sigma_pp": sigma_pp,
        "R_out_c": R_out_c,
        "u_rel": u_rel, "X": X,
        "n_solves": [0],
        "n0_t": n0_t, "m0_t": m0_t,
    }
    return P


# ============================================================================
# SOLVE PRIMITIVES
# ============================================================================
def set_initial_guess(P):
    """u <- identity map (x1, x2, 0, zeros). ONLY call for cold starts."""
    P["u"].interpolate(InitialGuess())


def _assign_strains(P, sv, sk):
    L = P["L_scale"]
    P["v0_c"].assign(Constant((float(sv[0]), float(sv[1]), float(sv[2]) + 1.0)))
    P["k0_c"].assign(Constant((float(sk[0]) * L, float(sk[1]) * L,
                               float(sk[2]) * L)))


def solve_at(P, sv, sk):
    """One SNES solve at (sv, sk) from whatever u currently holds (warm).
    Returns (niter, converged)."""
    _assign_strains(P, sv, sk)
    P["n_solves"][0] += 1
    return P["solver"].solve()


from contextlib import contextmanager

@contextmanager
def _state_guard(P):
    """Snapshot u, v0_c, k0_c on entry; restore them on exit (even on
    exception). Used by any routine that perturbs P's state temporarily
    (FD sweeps, stiffness probes) and must leave P exactly as found."""
    u = P["u"]
    u_saved  = Function(u.function_space())
    u_saved.assign(u)
    v0_saved = np.array(P["v0_c"].values())
    k0_saved = np.array(P["k0_c"].values())
    try:
        yield u_saved
    finally:
        u.assign(u_saved)
        P["v0_c"].assign(Constant(tuple(v0_saved)))
        P["k0_c"].assign(Constant(tuple(k0_saved)))


def solve_with_load_stepping(P, sv, sk, n_steps=10, verbose=False):
    """Cold solve: reset u, ramp strains 0 -> target linearly.
    Raises RuntimeError on failure."""
    sv = np.asarray(sv, dtype=float)
    sk = np.asarray(sk, dtype=float)
    set_initial_guess(P)
    for step in range(1, n_steps + 1):
        a = step / n_steps
        niter, conv = solve_at(P, a * sv, a * sk)
        if verbose:
            print(f"  load step {step}/{n_steps}: {niter} it, converged={conv}")
        if not conv:
            raise RuntimeError(f"load stepping failed at step {step}/{n_steps}")


# ============================================================================
# RESULTANTS  (Arora Eq. 26)
# ============================================================================
def compute_n0_m0(P):
    """n0 = int P e3 dA, m0 = int x_hat x P e3 dA on the CURRENT solution.
    Physical units (E already physical; L_scale factors included)."""
    L = P["L_scale"]
    n0 = np.array([assemble(P["Pe3_pp"][i] * dx) for i in range(3)]) * L ** 2
    m0 = np.array([assemble(P["xPe3_pp"][i] * dx) for i in range(3)]) * L ** 2
    return n0, m0


def forward_solve(P, sv, sk, warm=True, n_steps=10):
    """Full forward map (sv, sk) -> 6-vector [n0, m0]. sv/sk are the LOCAL
    (material-frame) rod strains -- do not rotate them before calling this.
    warm=True: solve directly from current u; falls back to load stepping
    if the warm solve fails."""
    if warm:
        niter, conv = solve_at(P, sv, sk)
        if not conv:
            solve_with_load_stepping(P, sv, sk, n_steps=n_steps)
    else:
        solve_with_load_stepping(P, sv, sk, n_steps=n_steps)
    n0, m0 = compute_n0_m0(P)
    return np.concatenate([n0, m0])


# ============================================================================
# STIFFNESS MAP -- work-conjugate resultant for each active strain component
# (FEniCS material frame throughout -- no R_local/R_MAP rotation applied here)
# ============================================================================
_STIFFNESS_MAP = {
    # label : (resultant_type, resultant_idx, human_name)
    "v1": ("n0", 0, "K_sv1 = dn0_1/dsv1"),
    "v2": ("n0", 1, "K_sv2 = dn0_2/dsv2"),
    "v3": ("n0", 2, "K_sv3 = dn0_3/dsv3"),
    "k1": ("m0", 0, "K_sk1 = dm0_1/dsk1"),
    "k2": ("m0", 1, "K_sk2 = dm0_2/dsk2"),
    "k3": ("m0", 2, "K_sk3 = dm0_3/dsk3"),
}


# ============================================================================
# LINEARIZED WARPING  (d(x_hat)/dq at an arbitrary base strain state)
# ----------------------------------------------------------------------------
# Call with sv_base=sk_base=zeros for scheme (ii) -- reference-state
# linearization matching Figs 5/6/8 of Arora et al. 2019.
# Call from inside a sweep loop with the current converged sv/sk for a
# per-step finite-strain linearization (closer to scheme iii).
#
# Saves and fully restores u, v0_c, k0_c -- zero side effect on caller.
# ============================================================================
def linearize_warping(P, active_type, active_component, sv_base, sk_base,
                      eps=1e-4, warm=False, n_steps=10):
    """
    d(x_hat)/d(active_component) via central FD around (sv_base, sk_base).

    active_type      : 'v' or 'k'
    active_component : 0, 1, or 2 (index into sv or sk)v
    sv_base, sk_base : base strain state for the FD (zeros = reference state)
    eps              : perturbation size (default 1e-4)
    warm             : if True, try warm solve first, fall back to load stepping
    n_steps          : load steps for cold solve

    Returns dxhat : (N, 3) float array, d(x_hat)/dq at every mesh node.
    """
    mesh = P["mesh"]
    u    = P["u"]

    sv_base = np.asarray(sv_base, dtype=float)
    sk_base = np.asarray(sk_base, dtype=float)

    def _solve_pert(pert):
        sv = sv_base.copy()
        sk = sk_base.copy()
        if active_type == 'v':
            sv[active_component] += pert
        else:
            sk[active_component] += pert
        if warm:
            niter, conv = solve_at(P, sv, sk)
            if not conv:
                solve_with_load_stepping(P, sv, sk, n_steps=n_steps)
        else:
            solve_with_load_stepping(P, sv, sk, n_steps=n_steps)
        return u.split(deepcopy=True)[0].compute_vertex_values(mesh).reshape(3, -1).T

    with _state_guard(P):
        x_plus  = _solve_pert(+eps)
        x_minus = _solve_pert(-eps)
        return (x_plus - x_minus) / (2.0 * eps)


# ============================================================================
# STIFFNESS FD  (single step -- called inside the sweep loop)
# ----------------------------------------------------------------------------
# Matches minimisation_2.py compute_K66_fd_at_step exactly:
#   - perturbs alpha by eps_fd (NOT the strain directly)
#   - forward difference, both solves warm-started from u_base
#   - K = (f(alpha+eps) - f(alpha)) / (eps_fd * q_target)
# Keeps u/v0_c/k0_c fully restored after the call.
# ============================================================================
def compute_stiffness_at_step(P, sv_step, sk_step, sv_target, sk_target,
                               active_label, eps_fd=1e-6):
    """
    Finite-difference tangent stiffness K = df_i/dq_i at the current
    load step (sv_step, sk_step), using alpha-based perturbation.

    Both FD solves are warm-started from u_base (the converged solution
    at this step), matching minimisation_2.py's approach exactly.

    active_label : one of 'v1'..'k3' -- determines which resultant
                   component is differentiated (FEniCS material frame).
    sv_target, sk_target : the full target strains (used to compute
                   alpha = sv_step / sv_target or sk_step / sk_target,
                   and q_target for the FD denominator).

    Returns K (float) or np.nan on failure.
    """
    if active_label not in _STIFFNESS_MAP:
        raise ValueError(f"active_label must be one of {list(_STIFFNESS_MAP)}")

    resultant_type, resultant_idx, stiffness_name = _STIFFNESS_MAP[active_label]
    active_type      = active_label[0]
    active_component = int(active_label[1]) - 1

    sv_target = np.asarray(sv_target, dtype=float)
    sk_target = np.asarray(sk_target, dtype=float)
    q_target  = sk_target[active_component] if active_type == 'k' \
                else sv_target[active_component]

    if abs(q_target) < 1e-4:
        return np.nan   # denominator too small

    u = P["u"]

    # infer current alpha from the step strains
    # (use the largest-magnitude component to avoid division by near-zero)
    _all = np.concatenate([sv_target, sk_target])
    _idx = int(np.argmax(np.abs(_all)))
    _q_t = _all[_idx]
    _q_s = np.concatenate([sv_step, sk_step])[_idx]
    a = _q_s / _q_t if abs(_q_t) > 1e-14 else 1.0

    try:
        with _state_guard(P) as u_base:
            # f_minus: resultant at alpha, warm from u_base
            u.assign(u_base)
            solve_at(P, sv_step, sk_step)
            n0_m, m0_m = compute_n0_m0(P)
            f_minus = n0_m if resultant_type == "n0" else m0_m

            # f_plus: resultant at alpha + eps_fd, warm from u_base
            sv_p = (a + eps_fd) * sv_target
            sk_p = (a + eps_fd) * sk_target
            u.assign(u_base)
            niter_p, conv_p = solve_at(P, sv_p, sk_p)
            if not conv_p:
                raise RuntimeError("stiffness FD perturbed solve failed")
            n0_p, m0_p = compute_n0_m0(P)
            f_plus = n0_p if resultant_type == "n0" else m0_p

            K = (f_plus[resultant_idx] - f_minus[resultant_idx]) \
                / (eps_fd * q_target)
            return K

    except RuntimeError:
        return np.nan


# ============================================================================
# INVERSE SOLVE  (Newton on strains, FD Jacobian)
# ============================================================================
def inverse_solve(P, target, sv0, sk0,
                  tol_rel=1e-8, tol_abs=1e-8, max_iter=25,
                  fd_rel=1e-6, fd_abs=1e-9, verbose=True):
    """Find (sv, sk) such that forward_solve == target (6-vector [n;m]).

    IMPORTANT: this function does NOT rotate `target`. It must already be
    expressed in the same LOCAL/material frame as sv/sk (i.e. the frame
    forward_solve's output n0/m0 comes back in, before any R_local
    rotation). If your target is in FEniCS axis order (Abaqus->FEniCS
    R_MAP already applied by hand), convert it in the driver script with
    target_to_material() BEFORE calling inverse_solve -- do not add
    that rotation here.

    Returns (x, f, history). x[:3]=sv, x[3:]=sk.
    """
    target = np.asarray(target, dtype=float)
    x = np.concatenate([np.asarray(sv0, float), np.asarray(sk0, float)])
    scale = max(np.linalg.norm(target), 1.0)

    # first solve cold (load stepping); all later solves warm
    f = forward_solve(P, x[:3], x[3:], warm=False)

    history = []
    for it in range(max_iter):
        r = f - target
        rnorm = np.linalg.norm(r)
        history.append({"iter": it, "x": x.copy(), "f": f.copy(),
                        "rnorm": rnorm})
        if verbose:
            print(f"\n--- Newton iter {it} ---")
            print(f"  sv = {x[:3]}")
            print(f"  sk = {x[3:]}")
            print(f"  f  = {f}")
            print(f"  |r| = {rnorm:.6e} (rel {rnorm/scale:.3e})")

        if rnorm < tol_abs or rnorm / scale < tol_rel:
            print(f"\nCONVERGED in {it} iterations, "
                  f"{P['n_solves'][0]} SNES solves total")
            return x, f, history

        # FD Jacobian: 6 warm-started perturbed solves
        J = np.zeros((6, 6))
        for p in range(6):
            eps = max(fd_abs, fd_rel * abs(x[p]))
            xp = x.copy()
            xp[p] += eps
            fp = forward_solve(P, xp[:3], xp[3:], warm=True)
            J[:, p] = (fp - f) / eps
        # restore base solution in u before stepping
        f = forward_solve(P, x[:3], x[3:], warm=True)

        asym = np.linalg.norm(J - J.T) / max(np.linalg.norm(J), 1e-30)
        if verbose:
            print(f"  Jacobian asymmetry: {asym:.3e}")

        try:
            dxs = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            dxs = np.linalg.lstsq(J, -r, rcond=None)[0]
            print("  (singular Jacobian, lstsq used)")

        # damped update with backtracking on |r|
        alpha_ls = 1.0
        accepted = False
        for _ in range(6):
            x_new = x + alpha_ls * dxs
            f_new = forward_solve(P, x_new[:3], x_new[3:], warm=True)
            if np.linalg.norm(f_new - target) < rnorm:
                accepted = True
                break
            alpha_ls *= 0.5
        if not accepted:
            print("  line search failed to reduce |r|; taking full step")
            x_new = x + dxs
            f_new = forward_solve(P, x_new[:3], x_new[3:], warm=True)

        x, f = x_new, f_new

    print(f"\nNOT converged after {max_iter} iterations")
    return x, f, history


def linear_beam_initial_guess(P, target):
    """Linear beam theory guess for (sv, sk) from target [n1,n2,n3,m1,m2,m3].

    Same frame requirement as inverse_solve: `target` must already be in
    the LOCAL/material frame. Convert with target_to_material() in the
    driver script first if needed.
    """
    L = P["L_scale"]
    x_coord = SpatialCoordinate(P["mesh"])
    A = P["total_area"] * L ** 2
    I11 = assemble((x_coord[1] - P["cy"]) ** 2 * dx) * L ** 4  # about e1
    I22 = assemble((x_coord[0] - P["cx"]) ** 2 * dx) * L ** 4  # about e2
    Jt = I11 + I22

    E = float(P["E_mod"])
    mu = float(P["mu_param"])

    n1, n2, n3, m1, m2, m3 = target
    sv = np.array([n1 / (mu * A), n2 / (mu * A), n3 / (E * A)])
    sk = np.array([m1 / (E * I11), m2 / (E * I22), m3 / (mu * Jt)])
    return sv, sk


# ============================================================================
# FRAMES
# ----------------------------------------------------------------------------
# The core solver works entirely in the FEniCS material frame (ei, s=0, R=I).
# No Abaqus frame conversion here. Driver scripts handle that.
#
# R_local (from UR) lives only in driver scripts:
#   input:  sv_mat = R_local.T @ sv_spatial   (pull-back strains)
#   output: v = R_local @ sv_mat              (push-forward, inverse_solve)
#
# postprocess accepts optional R_out for rotating outputs (stress/disp/VTU)
# to any desired frame. Pass R_MAP from the Abaqus driver to get Abaqus-frame
# stress and displacement in the VTU. Pass None for pure FEniCS frame output.
# ============================================================================

def rotation_matrix_from_ur(ur1, ur2, ur3):
    """R from an Abaqus rotation VECTOR (axis-angle), via Rodrigues."""
    theta = np.array([ur1, ur2, ur3], float)
    phi = np.linalg.norm(theta)
    if phi < 1e-12:
        return np.eye(3)
    k = theta / phi
    K = np.array([[    0, -k[2],  k[1]],
                  [ k[2],     0, -k[0]],
                  [-k[1],  k[0],    0]])
    return np.eye(3) + np.sin(phi)*K + (1 - np.cos(phi))*(K @ K)

def target_to_material(target_spatial, R_local):
    """Pull back 6-vector [n;m] from spatial to FEniCS material frame.
    Use in driver scripts before inverse_solve.
    """
    R = np.asarray(R_local, dtype=float)
    t = np.asarray(target_spatial, dtype=float)
    return np.concatenate([R.T @ t[:3], R.T @ t[3:]])


# ============================================================================
# POST-PROCESSING
# ----------------------------------------------------------------------------
# All physics quantities computed in FEniCS material frame.
# R_out (optional): rotation matrix applied to stress/displacement OUTPUTS only
# (e.g. pass R_MAP from the Abaqus driver to get Abaqus-frame VTU/plots/CSV).
# No physics is affected -- only how results are expressed in outputs.
# ============================================================================
def postprocess(P, OUTDIR, base_name, R_local=None, R_out=None, do_plots=True,
                do_vtu=True, n_line_pts=100):
    """Post-processing. Physics assembled in FEniCS material frame; this is
    the ONE place both rotations get applied, composed, to produce the final
    physically-correct output -- driver scripts do not rotate anything
    themselves afterward.

    R_local : optional (3,3) rotation, this station's deformed/spatial
              rotation (from UR). Applied FIRST.
    R_out   : optional (3,3) rotation, axis convention (e.g. R_MAP,
              FEniCS->Abaqus). Applied SECOND, on top of R_local.
    Composed as R_total = R_out @ R_local (both None/identity -> no
    rotation, pure FEniCS material frame). Applied to:
              - n0, m0      one-sided  (R_total @ v)      -- vectors
              - sigma, E    two-sided  (R_total @ T @ R_total.T) -- tensors
              - displacement, VTU mesh points  one-sided (R_total @ v)
    The material-frame n0/m0 are always also printed/saved (result["n0"],
    result["m0"]) -- that's the frame the SNES solve and inverse_solve
    Newton loop actually work in, useful as a sanity check. The fully
    rotated resultants are result["n0_out"]/result["m0_out"] -- use these
    directly for Abaqus comparison; no further rotation needed downstream.

    Sections:
      1. resultants   n0, m0            (material frame + fully rotated)
      2. decomposition  Pe3 term-by-term
      3. fields         sigma, E -> CG1  (fully rotated)
      4. line query     centerline CSV
      5. plots          stress/strain centerline
      6. vtu            displacement + stress + final deformed mesh (rotated)
      7. checks         Lagrange multipliers, displacement range

    Returns dict with n0, m0, n0_out, m0_out, decomposition, sig_cg, strain_cg.
    """
    os.makedirs(OUTDIR, exist_ok=True)
    mesh = P["mesh"]
    L    = P["L_scale"]
    u    = P["u"]

    # output rotation: assign the COMPOSED matrix into the ONE persistent
    # Constant on P rather than building a new Constant(...) here -- see
    # build_problem's R_out_c comment. Only R_out_c itself needs to
    # persist; sigma_out/E_out below are rebuilt fresh from it every call
    # at no extra cost.
    R_local_np = np.asarray(R_local, dtype=float) if R_local is not None else np.eye(3)
    R_out_np_in = np.asarray(R_out, dtype=float) if R_out is not None else np.eye(3)
    R_total_np = R_out_np_in @ R_local_np   # R_local first, R_out (R_MAP) second
    frame_tag  = "rotated" if (R_local is not None or R_out is not None) else "FEniCS"
    P["R_out_c"].assign(Constant(R_total_np))
    print(f"R_total_np =\n{P['R_out_c'].values()}")
    # u.split(deepcopy=True) is needed ONLY for point-evaluation / VTU output
    # below (Point()/compute_vertex_values need standalone Functions). It is
    # NOT used for physics: F, E, S, sigma, n0, m0 all come straight from the
    # persistent UFL forms in P, built once in build_problem against the
    # SAME u_rel that drove the SNES solve. This guarantees postprocess()
    # reports exactly the state forward_solve/inverse_solve converged to --
    # no second, possibly-inconsistent deformation-gradient build.
    u_split   = u.split(deepcopy=True)
    x_sol     = u_split[0]

    # rotated outputs: stress and strain rotate two-sided (rank-2 tensors,
    # R @ T @ R.T); resultants rotate one-sided (vectors, R @ v).
    # Displacement (section 6) is one-sided too. Rebuilding this expression
    # every call is cheap -- R_out_c and sigma_pp/E_strain are the same
    # persistent objects each time, so the form signature is unchanged.
    R_out_c   = P["R_out_c"]
    sigma_out = R_out_c * P["sigma_pp"] * R_out_c.T
    E_out     = R_out_c * P["E_strain"] * R_out_c.T

    # ========================================================================
    # 1. RESULTANTS -- computed in FEniCS material frame (Arora Eq. 26, same
    #    ground truth compute_n0_m0() gives forward_solve/inverse_solve), then
    #    rotated one-sided into the fully composed (R_local then R_out) frame.
    # ========================================================================
    n0, m0 = compute_n0_m0(P)
    n0_out = R_total_np @ n0
    m0_out = R_total_np @ m0

    print("\n" + "=" * 60)
    print("RESULTANTS (FEniCS material frame)")
    print(f"  n0 = [{n0[0]: .6e}, {n0[1]: .6e}, {n0[2]: .6e}]")
    print(f"  m0 = [{m0[0]: .6e}, {m0[1]: .6e}, {m0[2]: .6e}]")
    if R_local is not None or R_out is not None:
        print(f"RESULTANTS ({frame_tag} frame)")
        print(f"  n0 = [{n0_out[0]: .6e}, {n0_out[1]: .6e}, {n0_out[2]: .6e}]")
        print(f"  m0 = [{m0_out[0]: .6e}, {m0_out[1]: .6e}, {m0_out[2]: .6e}]")

    pd.DataFrame([{"n1": n0[0], "n2": n0[1], "n3": n0[2],
                   "m1": m0[0], "m2": m0[1], "m3": m0[2],
                   "n1_out": n0_out[0], "n2_out": n0_out[1], "n3_out": n0_out[2],
                   "m1_out": m0_out[0], "m2_out": m0_out[1], "m3_out": m0_out[2]}]).to_csv(
        os.path.join(OUTDIR, f"{base_name}_resultants.csv"), index=False)

    # ========================================================================
    # 2. DECOMPOSITION
    # ========================================================================
    # F = P["F_rod"]
    # S = P["S_pp"]
    # dec = {
    #     "shear_F22_S23_ordinary": float(assemble(F[1,1]*S[1,2]*dx)) * L**2,
    #     "shear_F23_S33_leakage":  float(assemble(F[1,2]*S[2,2]*dx)) * L**2,
    #     "axial_F33_S33_dominant": float(assemble(F[2,2]*S[2,2]*dx)) * L**2,
    #     "axial_F32_S23_w2_leak":  float(assemble(F[2,1]*S[1,2]*dx)) * L**2,
    #     "axial_F31_S13_w1_leak":  float(assemble(F[2,0]*S[0,2]*dx)) * L**2,
    # }
    # print(f"\n--- Pe3 decomposition ({frame_tag} frame) ---")
    # for k, v in dec.items():
    #     print(f"  {k}: {v: .6e}")
    # pd.DataFrame([dec]).to_csv(
    #     os.path.join(OUTDIR, f"{base_name}_decomposition.csv"), index=False)

    # ========================================================================
    # 3. FIELDS -- CG1 projection of sigma_out (rotated if R_out given)
    # ========================================================================
    V_cg = FunctionSpace(mesh, "CG", 1)

    sig_cg = {f"S{i+1}{j+1}": project(sigma_out[i, j], V_cg)
              for i in range(3) for j in range(i, 3)}
    strain_cg = {f"E{i+1}{j+1}": project(E_out[i, j], V_cg)
                 for i in range(3) for j in range(i, 3)}

    Mises = sqrt(0.5 * ((sigma_out[0,0]-sigma_out[1,1])**2 +
                        (sigma_out[1,1]-sigma_out[2,2])**2 +
                        (sigma_out[2,2]-sigma_out[0,0])**2) +
                 3*(sigma_out[0,1]**2 + sigma_out[1,2]**2 + sigma_out[0,2]**2))
    sig_cg["Mises"] = project(Mises, V_cg)
    sig_cg["inv_1"] = project(tr(sigma_out), V_cg)

    # ========================================================================
    # 4. LINE QUERY -- centerline x1=0
    # ========================================================================
    X_mesh  = mesh.coordinates()
    all_fields = {**sig_cg, **strain_cg}
    # Quad-mesh-safe line query: DOLFIN 2019.1.0 point-location
    # (Point()/Function.__call__) only supports simplex cells, so instead
    # of sampling n_line_pts arbitrary points we use the mesh VERTICES that
    # already lie on x1=0 -- exact for a structured, centered grid where
    # that line coincides with a vertex column (small tolerance for
    # scale()/translate() roundoff). Falls back to the old Point()-based
    # sampling (works on simplex meshes, degrades to NaN on quads) if no
    # such vertex column is found.
    on_line = np.where(np.abs(X_mesh[:, 0]) < 1e-9)[0]
    if on_line.size > 0:
        order = on_line[np.argsort(X_mesh[on_line, 1])]
        x2_vals = X_mesh[order, 1]
        line = {}
        for k, fn in all_fields.items():
            vv = fn.compute_vertex_values(mesh)
            line[k] = vv[order]
    else:
        x2_vals = np.linspace(X_mesh[:,1].min(), X_mesh[:,1].max(), n_line_pts)
        line = {k: np.full(n_line_pts, np.nan) for k in all_fields}
        for i, x2 in enumerate(x2_vals):
            pt = Point(0.0, x2)
            for k, fn in all_fields.items():
                try:
                    line[k][i] = float(fn(pt))
                except RuntimeError:
                    pass
    pd.DataFrame({"x1": 0.0, "x2": x2_vals * L, **line}).to_csv(
        os.path.join(OUTDIR, f"{base_name}_section_X0.csv"),
        index=False, float_format="%.6e")

    # ========================================================================
    # 5. PLOTS
    # ========================================================================
    if do_plots:
        skeys   = ["S11","S12","S13","S22","S23","S33",
                   "Mises","inv_1"]
        slabels = {"S11":"S11","S12":"S12","S13":"S13",
                   "S22":"S22","S23":"S23","S33":"S33",
                   "Mises":"von Mises","inv_1":"First Invariant"}
        scolors = {"S11":"#1f77b4","S12":"#d62728","S13":"#9467bd",
                   "S22":"#8c564b","S23":"#e377c2","S33":"#7f7f7f",
                   "Mises":"#2ca02c","inv_1":"#ff7f0e"}
        fig, axes = plt.subplots(2, 4, figsize=(24, 8))
        fig.suptitle(f"Stress ({frame_tag} frame)")
        for ax, k in zip(axes.flatten(), skeys):
            ax.plot(line[k], x2_vals*L, color=scolors[k])
            ax.set_title(slabels[k])
            ax.axvline(0, color="k", lw=0.6)
            ax.grid(True)
        axes[0,0].set_ylabel(r"$X_2$")
        plt.tight_layout()
        fig.savefig(os.path.join(OUTDIR, f"{base_name}_Stress_centerline.png"),
                    dpi=150)
        plt.close(fig)

        ekeys = ["E11","E12","E13","E22","E23","E33"]
        fig2, axes2 = plt.subplots(2, 3, figsize=(18, 8))
        fig2.suptitle(f"Green-Lagrange strain ({frame_tag} frame)")
        for ax, k in zip(axes2.flatten(), ekeys):
            ax.plot(line[k], x2_vals*L)
            ax.set_title(k)
            ax.axvline(0, color="k", lw=0.6)
            ax.grid(True)
        plt.tight_layout()
        fig2.savefig(os.path.join(OUTDIR, f"{base_name}_Strain_centerline.png"),
                     dpi=150)
        plt.close(fig2)

    # ========================================================================
    # 6/7 shared: evaluate x_sol at every mesh node ONCE (was done twice --
    # identically -- in the old VTU and CHECKS sections).
    # ========================================================================
    num_nodes = X_mesh.shape[0]
    # Quad-mesh-safe node evaluation: compute_vertex_values() reads the
    # P2 field through each cell's own local basis (vertex coords are a
    # subset of the P2 nodes), not via point-location, so it works on
    # quads where Point()/Function.__call__ does not.
    _x_vv = x_sol.compute_vertex_values(mesh)   # flat, length 3*num_nodes
    x_arr = _x_vv.reshape(3, num_nodes).T        # -> (num_nodes, 3)
    P0    = np.column_stack([X_mesh[:,0]*L, X_mesh[:,1]*L,
                              np.zeros(num_nodes)])
    U_fen = x_arr*L - P0            # displacement, FEniCS frame
    U_out = (R_total_np @ U_fen.T).T  # rotated to output frame

    # Reference and current geometry live in the SAME local frame as U --
    # rotating U alone and leaving P0 unrotated mixes frames the moment
    # R_total != Identity. P0 must get the same rotation as U.
    P0_out = (R_total_np @ P0.T).T

    # ========================================================================
    # 6. VTU -- reference mesh + displacement field, and the final
    #    (deformed) mesh directly, both rotated by R_out if given.
    # ========================================================================
    if do_vtu:
        nodal = {k: fn.compute_vertex_values(mesh)
                for k, fn in {**sig_cg, **strain_cg}.items()}
        _cell_name = {"triangle": "triangle",
                       "quadrilateral": "quad"}[mesh.ufl_cell().cellname()]
        cells = [(_cell_name, mesh.cells())]

        P0_ref = (R_out_np_in @ P0.T).T           # reference: R_MAP only (no R_local)
        x_def  = (R_total_np @ (x_arr * L).T).T    # deformed: full R_total on current coords

        # REFERENCE config -> compare THIS against the Abaqus z=50 slice. No warp.
        meshio.write(os.path.join(OUTDIR, f"{base_name}_ref.vtu"),
                    meshio.Mesh(points=P0_ref, cells=cells,
                                point_data={"U": U_out, **nodal}))

        # DEFORMED config -> matches Abaqus deformed positions.
        meshio.write(os.path.join(OUTDIR, f"{base_name}_deformed.vtu"),
                    meshio.Mesh(points=x_def, cells=cells,
                                point_data={"U": U_out, **nodal}))
    # ========================================================================
    # 7. CHECKS
    # ========================================================================
    lam_vals = np.array([u_split[i].vector().get_local()[0]
                         for i in [1,2,3]])
    mu_vals  = np.array([u_split[i].vector().get_local()[0]
                         for i in [4,5,6]])
    print(f"\nLagrange multipliers: lambda={lam_vals}  mu={mu_vals}")

    print(f"x1 diff: {np.max(np.abs(x_arr[:,0]-X_mesh[:,0])):.4e}")
    print(f"x2 diff: {np.max(np.abs(x_arr[:,1]-X_mesh[:,1])):.4e}")
    print(f"x3 range: {np.max(np.abs(x_arr[:,2])):.4e}")

    # f_cauchy_ref_area = np.array([assemble(sigma[i,2]*dx) for i in range(3)]) * L**2
    # print(f"Cauchy forces (FEniCS): {f_cauchy_ref_area}")
    # print(f"vs n0 (1stPK):          {n0}")

    return {"n0": n0, "m0": m0, "n0_out": n0_out, "m0_out": m0_out,
            #  "decomposition": dec,
            "sig_cg": sig_cg, "strain_cg": strain_cg}
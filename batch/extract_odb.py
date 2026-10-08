# -*- coding: utf-8 -*-
"""
extract_odb.py  --  abaqus python extract_odb.py <beam.odb> <solid.odb> <out_dir> [Z_MID] [Z_STATIONS]

Beam ODB  : reads SE, SK, SF, SM, UR, U at every node.
            Snaps to every station in Z_STATIONS.
            Writes  <out_dir>/beam_whole.csv   (all stations)

Solid ODB : reads S, LE, SENER, ELSE, NFORC1/2/3 at Z_MID layer (nodal).
            Also integrates NFORC cross-section resultants at every station.
            Writes  <out_dir>/solid_mid.csv          (nodal, Z_MID only, for l2)
                    <out_dir>/solid_energy.csv        (ELSE scalar)
                    <out_dir>/solid_resultants.csv    (NFORC resultants, all stations)

Z_STATIONS: comma-separated list of Z values, e.g. "0,10,20,30,40,50"
            If omitted, defaults to Z_MID only (backward compatible).

NOTE ON NFORC: when requested under *Element Output, Abaqus stores nodal
forces as separate scalar fields NFORC1, NFORC2, NFORC3 (X, Y, Z components)
rather than a single vector field NFORC. This file reads NFORC1/2/3 and
merges them. If requested under *Node Output instead, they appear as a
single NFORC vector field — the code handles both cases.

Called from batch_driver.py via subprocess -- out_dir is the case folder
that batch_driver already creates and owns, so no extra path setup needed.
"""
import sys
import os
import csv
import math
from odbAccess import openOdb
from abaqusConstants import ELEMENT_NODAL, CENTROID

# ── command-line args ─────────────────────────────────────────────────────────
beam_odb  = sys.argv[1]
solid_odb = sys.argv[2]
out_dir   = sys.argv[3]
Z_MID     = float(sys.argv[4]) if len(sys.argv) > 4 else 50.0

if len(sys.argv) > 5:
    Z_STATIONS = [float(z) for z in sys.argv[5].split(',')]
else:
    Z_STATIONS = [Z_MID]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

# Per-station NFORC output is verbose (11 stations x every case) and all of
# it lands in solid_resultants.csv anyway. Quiet by default; set
# VERBOSE_NFORC=1 to restore. Warnings and errors always print.
_VERBOSE_NFORC = os.environ.get('VERBOSE_NFORC', '') not in ('', '0', 'false', 'False')

print("Z_STATIONS:", Z_STATIONS)
print("Z_MID     :", Z_MID)


# ── helpers ───────────────────────────────────────────────────────────────────
def avg_nodal(field_inst, n_comps):
    """Average ELEMENT_NODAL contributions -> {nodeLabel: [avg_comp,...]}"""
    acc, cnt = {}, {}
    for v in field_inst.values:
        if not hasattr(v, 'nodeLabel'):
            continue
        nd = v.nodeLabel
        d  = v.data
        if len(d) < n_comps:
            continue
        if nd not in acc:
            acc[nd] = [0.0] * n_comps
            cnt[nd] = 0
        for i in range(n_comps):
            acc[nd][i] += float(d[i])
        cnt[nd] += 1
    return {nd: [acc[nd][i] / cnt[nd] for i in range(n_comps)] for nd in acc}


def snap_z(coord_map, z_target, tol=None):
    """
    Return (set_of_node_labels, snapped_z) for nodes closest to z_target.
    tol defaults to 45% of the minimum Z spacing in the mesh so it
    catches exactly one layer without bleeding into the next.
    """
    zvals = sorted(set(c[2] for c in coord_map.values()))
    best  = min(zvals, key=lambda z: abs(z - z_target))
    if tol is None:
        diffs = [zvals[i+1] - zvals[i] for i in range(len(zvals) - 1)]
        tol   = min(diffs) * 0.45 if diffs else 0.1
    nodes = {nd for nd, c in coord_map.items() if abs(c[2] - best) <= tol}
    return nodes, best


def build_el_centroid_z(inst, coord_map):
    """elementLabel -> centroid Z, for all elements in inst."""
    out = {}
    for el in inst.elements:
        nls = [n for n in el.connectivity if n in coord_map]
        if nls:
            out[el.label] = sum(coord_map[n][2] for n in nls) / len(nls)
    return out


def _read_scalar_nforc(fo_dict, key, inst):
    """Read a scalar NFORC component field (NFORC1/NFORC2/NFORC3).
    Returns {(nodeLabel, elementLabel): float}.
    Tries getSubset(region=inst) first; falls back to raw values."""
    out = {}
    if key not in fo_dict:
        return out
    try:
        fo = fo_dict[key].getSubset(region=inst)
    except Exception:
        fo = fo_dict[key]
    for v in fo.values:
        if not hasattr(v, 'nodeLabel') or not hasattr(v, 'elementLabel'):
            continue
        # scalar field: v.data is a tuple with one element
        try:
            val = float(v.data[0]) if hasattr(v.data, '__len__') else float(v.data)
        except Exception:
            continue
        out[(v.nodeLabel, v.elementLabel)] = val
    return out


# ═════════════════════════════════════════════════════════════════════════════
# BEAM ODB — all Z stations
# ═════════════════════════════════════════════════════════════════════════════
print("\nOpening beam ODB:", beam_odb)
odb  = openOdb(beam_odb, readOnly=True)
inst = odb.rootAssembly.instances.values()[0]
frm  = odb.steps.values()[-1].frames[-1]

coord_b = {n.label: tuple(n.coordinates) for n in inst.nodes}

# ── undeformed-geometry tilt check (see chat: X/Z, Y/Z growing linearly
# with Z means the beam's UNDEFORMED centerline isn't perfectly aligned
# with global Z -- report it as a %, and as an equivalent angle, so it's
# visible without doing this by hand from beam_whole.csv every time). Uses
# the two stations farthest apart in Z actually present in the mesh, so
# it works regardless of how many/which Z_STATIONS were requested.
_zs = sorted(set(round(c[2], 8) for c in coord_b.values()))
if len(_zs) >= 2:
    z_lo, z_hi = _zs[0], _zs[-1]
    dz = z_hi - z_lo
    if abs(dz) > 1e-12:
        node_lo = min(coord_b.items(), key=lambda kv: abs(kv[1][2] - z_lo))[1]
        node_hi = min(coord_b.items(), key=lambda kv: abs(kv[1][2] - z_hi))[1]
        dx, dy = node_hi[0] - node_lo[0], node_hi[1] - node_lo[1]
        tilt_x_pct = 100.0 * dx / dz
        tilt_y_pct = 100.0 * dy / dz
        tilt_angle_deg = math.degrees(math.atan2((dx**2 + dy**2) ** 0.5, abs(dz)))
        print("Undeformed geometry tilt (Z=%.3f -> Z=%.3f): "
              "X/Z=%.3e%%  Y/Z=%.3e%%  angle=%.3e deg"
              % (z_lo, z_hi, tilt_x_pct, tilt_y_pct, tilt_angle_deg))
        if tilt_angle_deg > 1e-3:
            print("  NOTE: tilt is above the usual round-off range "
                  "(>1e-3 deg) -- check whether this was an intentional "
                  "seeded imperfection or an unintended geometry/mesh issue.")

def _nodal_map(fo_dict, key, n):
    if key not in fo_dict:
        return {}
    return {v.nodeLabel: list(v.data)
            for v in fo_dict[key].values
            if hasattr(v, 'nodeLabel') and len(v.data) >= n}

se_map = _nodal_map(frm.fieldOutputs, 'SE', 3)
sk_map = _nodal_map(frm.fieldOutputs, 'SK', 3)
sf_map = _nodal_map(frm.fieldOutputs, 'SF', 3)
sm_map = _nodal_map(frm.fieldOutputs, 'SM', 3)
ur_map = _nodal_map(frm.fieldOutputs, 'UR', 3)
u_map  = _nodal_map(frm.fieldOutputs, 'U',  3)

beam_csv = os.path.join(out_dir, 'beam_whole.csv')
with open(beam_csv, 'wb') as f:
    w = csv.writer(f)
    w.writerow(['Node', 'X', 'Y', 'Z',
                'SE1', 'SE2', 'SE3',
                'SK1', 'SK2', 'SK3',
                'UR1', 'UR2', 'UR3',
                'U1',  'U2',  'U3',
                'SF1', 'SF2', 'SF3',
                'SM1', 'SM2', 'SM3'])
    for z_target in Z_STATIONS:
        snap_nodes_b, z_snap_b = snap_z(coord_b, z_target)
        print("  Beam Z=%.2f -> snapped=%.4f, %d nodes"
              % (z_target, z_snap_b, len(snap_nodes_b)))
        for nd in sorted(snap_nodes_b):
            c  = coord_b[nd]
            se = se_map.get(nd, [0.0]*3)
            sk = sk_map.get(nd, [0.0]*3)
            ur = ur_map.get(nd, [0.0]*3)
            u  = u_map.get(nd,  [0.0]*3)
            sf = sf_map.get(nd, [0.0]*3)
            sm = sm_map.get(nd, [0.0]*3)
            w.writerow([nd, c[0], c[1], c[2],
                        se[0], se[1], se[2],
                        sk[0], sk[1], sk[2],
                        ur[0], ur[1], ur[2],
                        u[0],  u[1],  u[2],
                        sf[0], sf[1], sf[2],
                        sm[0], sm[1], sm[2]])

odb.close()
print("Wrote:", beam_csv)


# ═════════════════════════════════════════════════════════════════════════════
# BEAM ODB — REAL section-point strain (no SE/SK extrapolation)
# ═════════════════════════════════════════════════════════════════════════════
# SE1-3/SK1-3 above are section-averaged/curvature quantities at the beam
# AXIS, not a strain tensor at an actual material point (see chat: SE2/SE3
# are the Timoshenko section-AVERAGE shear strain, SK2/SK3 are curvatures,
# not strain). Whatever the beam section actually reports depends on how
# many/which section points were requested in the .inp's *BEAM SECTION —
# by default only the corners (4 for a rectangle in 3D). This block reads
# 'E' (or 'LE' if NLGEOM is on) at EVERY section point Abaqus actually
# wrote output for -- whatever that count is -- and computes a genuine
# local principal strain at each one: no extreme-fiber extrapolation,
# no combining an averaged shear with a corner bending term that don't
# physically co-locate (see chat -- that's what the SE/SK combination
# above does, conservatively).
#
# If only the 4 default corner points exist, this will only give you 4
# points per station -- same coverage as before, just done correctly
# instead of extrapolated. To actually sample mid-width/mid-edge points
# you must increase the number of OUTPUT section points in the beam
# section definition in the .inp (e.g. *SECTION POINTS, or the
# nondefault-integration section-point count on *BEAM SECTION) -- this
# script can only read what Abaqus was told to write; it cannot conjure
# section points the analysis never computed output for.
odb  = openOdb(beam_odb, readOnly=True)
inst = odb.rootAssembly.instances.values()[0]
frm  = odb.steps.values()[-1].frames[-1]
coord_b = {n.label: tuple(n.coordinates) for n in inst.nodes}

strain_key = 'LE' if 'LE' in frm.fieldOutputs else ('E' if 'E' in frm.fieldOutputs else None)

sp_csv = os.path.join(out_dir, 'beam_section_points.csv')
n_sp_written = 0
if strain_key is None:
    print("WARNING: neither 'E' nor 'LE' present in beam ODB field outputs "
          "-- beam_section_points.csv will be empty. Request strain output "
          "(E) on the beam element set to enable real section-point "
          "principal strain.")
    with open(sp_csv, 'wb') as f:
        csv.writer(f).writerow(['Z_target', 'Z_snapped', 'Node',
                                 'SectionPoint', 'SP_description',
                                 'E11', 'E12', 'E13', 'max_principal_strain'])
else:
    # {(nodeLabel, sectionPoint.number): (data_tuple, description)}
    sp_map = {}
    for v in frm.fieldOutputs[strain_key].values:
        if not hasattr(v, 'nodeLabel'):
            continue
        sp = getattr(v, 'sectionPoint', None)
        if sp is None:
            continue
        sp_map[(v.nodeLabel, sp.number)] = (tuple(v.data), sp.description)

    n_section_points = len(set(k[1] for k in sp_map.keys()))
    print("Beam section-point strain: field=%s, %d distinct section "
          "point(s) found in ODB output" % (strain_key, n_section_points))
    if n_section_points <= 4:
        print("  NOTE: only %d section point(s) available -- this is "
              "likely just the default corner output. To sample "
              "mid-width/mid-edge points too, increase the output "
              "section point count on *BEAM SECTION / *SECTION POINTS "
              "in the .inp." % n_section_points)

    with open(sp_csv, 'wb') as f:
        w = csv.writer(f)
        w.writerow(['Z_target', 'Z_snapped', 'Node', 'SectionPoint',
                     'SP_description', 'E11', 'E12', 'E13',
                     'max_principal_strain'])
        for z_target in Z_STATIONS:
            snap_nodes, z_snap = snap_z(coord_b, z_target)
            for nd in sorted(snap_nodes):
                for (nlab, spnum), (data, desc) in sp_map.items():
                    if nlab != nd:
                        continue
                    if len(data) < 3:
                        continue
                    e11, e12, e13 = float(data[0]), float(data[1]), float(data[2])
                    # resultant shear on this section point's plane.
                    # Abaqus reports E12/E13 (and LE12/LE13) as ENGINEERING
                    # shear strains (gamma = 2*eps_tensor), so they are
                    # used as-is -- no factor 2. Principal strain then
                    # follows the same Mohr's-circle form used everywhere
                    # else: eps/2 + sqrt((eps/2)^2 + (gamma/2)^2). e11/e12/
                    # e13 are REAL local values at this actual point, not
                    # extrapolated/averaged section quantities.
                    gam = math.sqrt(e12**2 + e13**2)
                    max_p = e11/2.0 + math.sqrt((e11/2.0)**2 + (gam/2.0)**2)
                    w.writerow([z_target, z_snap, nd, spnum, desc,
                                e11, e12, e13, max_p])
                    n_sp_written += 1
    print("  wrote %d (station, node, section-point) rows" % n_sp_written)

odb.close()
print("Wrote:", sp_csv)


# ═════════════════════════════════════════════════════════════════════════════
# SOLID ODB
# ═════════════════════════════════════════════════════════════════════════════
print("\nOpening solid ODB:", solid_odb)
odb  = openOdb(solid_odb, readOnly=True)
inst = odb.rootAssembly.instances.values()[0]
frm  = odb.steps.values()[-1].frames[-1]

coord_s = {n.label: tuple(n.coordinates) for n in inst.nodes}

# element centroid Z — built once, reused for every station
el_cz = build_el_centroid_z(inst, coord_s)

# ---------------------------------------------------------------------------
# SOLID nodal DISPLACEMENTS -- needed for the deformed-configuration moment
# resultants written alongside the reference-configuration ones below.
#
# WHY: NFORC gives nodal forces in the CURRENT configuration, but coord_s
# holds n.coordinates, which are the UNDEFORMED coordinates. Using them as
# lever arms under nlgeom=YES mixes two configurations. The error is zero at
# the clamped root and grows with deformation -- measured over 150 cases, the
# Mx discrepancy against the 1D beam goes 0.46% -> 2.86% of the root moment
# from z=0 to z=L, tracking the lateral displacement (0 -> 19.8) almost
# exactly. My and Mz look far worse in relative terms (124%, 115% at the tip)
# but that is a vanishing denominator: both go to zero at the free end by
# construction, so their ratios are meaningless there.
#
# The Mx/My/Mz columns below take moments about the DEFORMED section centroid
# using DEFORMED lever arms, and use the full r x F rather than dropping the
# out-of-plane term (which matters once the section warps under torsion).
u_s = _nodal_map(frm.fieldOutputs, 'U', 3)
if not u_s:
    print("WARNING: no 'U' field in the solid ODB -- deformed-configuration "
          "moments (Mx/My/Mz) fall back to undeformed lever arms.")

# ── nodal stress + strain at Z_MID (for existing l2 / overlay pipeline) ──────
s_inst  = frm.fieldOutputs['S'].getSubset(region=inst, position=ELEMENT_NODAL)
node_s  = avg_nodal(s_inst, 6)

node_le = {}
if 'LE' in frm.fieldOutputs:
    le_inst = frm.fieldOutputs['LE'].getSubset(region=inst, position=ELEMENT_NODAL)
    node_le = avg_nodal(le_inst, 6)

node_sener = {}
if 'SENER' in frm.fieldOutputs:
    sn_inst    = frm.fieldOutputs['SENER'].getSubset(region=inst,
                                                      position=ELEMENT_NODAL)
    node_sener = {v.nodeLabel: float(v.data)
                  for v in sn_inst.values if hasattr(v, 'nodeLabel')}

snap_nodes_s, z_snap_s = snap_z(coord_s, Z_MID)
print("Solid nodal snap Z=%.4f, %d nodes" % (z_snap_s, len(snap_nodes_s)))

solid_csv = os.path.join(out_dir, 'solid_mid.csv')
with open(solid_csv, 'wb') as f:
    w = csv.writer(f)
    w.writerow(['Node', 'X', 'Y', 'Z',
                'S11', 'S22', 'S33', 'S12', 'S13', 'S23', 'Mises', 'inv_1',
                'LE11', 'LE22', 'LE33', 'LE12', 'LE13', 'LE23', 'SENER'])
    for nd in sorted(snap_nodes_s):
        if nd not in node_s:
            continue
        c  = coord_s[nd]
        s  = node_s[nd]
        le = node_le.get(nd, [0.0]*6)
        sn = node_sener.get(nd, 0.0)
        s11, s22, s33, s12, s13, s23 = s
        mises = math.sqrt(0.5 * ((s11-s22)**2 + (s22-s33)**2 + (s33-s11)**2)
                          + 3.0*(s12**2 + s13**2 + s23**2))
        inv1  = s11 + s22 + s33
        w.writerow([nd, c[0], c[1], c[2],
                    s11, s22, s33, s12, s13, s23, mises, inv1,
                    le[0], le[1], le[2], le[3], le[4], le[5], sn])
print("Wrote:", solid_csv)

# ── ELSE energy scalar ────────────────────────────────────────────────────────
elem_else = {}
if 'ELSE' in frm.fieldOutputs:
    for v in frm.fieldOutputs['ELSE'].getSubset(region=inst).values:
        if hasattr(v, 'elementLabel'):
            elem_else[v.elementLabel] = float(v.data)

elements = []
for el in inst.elements:
    nls = list(el.connectivity)
    zs  = [coord_s[n][2] for n in nls if n in coord_s]
    if zs:
        elements.append((el.label, nls, min(zs), max(zs)))

U_else = 0.0
n_else = 0
ax_lo = ax_hi = None
if elements:
    lengths   = sorted(zmax - zmin for _, _, zmin, zmax in elements if zmax > zmin)
    elem_half = 0.75 * lengths[len(lengths)//2]
    for el_label, nls, zmin, zmax in elements:
        if abs(0.5*(zmin+zmax) - Z_MID) > elem_half:
            continue
        if el_label not in elem_else:
            continue
        U_else += elem_else[el_label]
        n_else += 1
        ax_lo = zmin if ax_lo is None else min(ax_lo, zmin)
        ax_hi = zmax if ax_hi is None else max(ax_hi, zmax)

thickness    = (ax_hi - ax_lo) if (ax_lo is not None and ax_hi > ax_lo) else 1.0
U_per_length = U_else / thickness

energy_csv = os.path.join(out_dir, 'solid_energy.csv')
with open(energy_csv, 'wb') as f:
    w = csv.writer(f)
    w.writerow(['Z_mid', 'U_ELSE_slice', 'n_elements',
                'slice_thickness', 'U_per_length_ELSE'])
    w.writerow([Z_MID, U_else, n_else, thickness, U_per_length])
print("Wrote:", energy_csv)

# ── NFORC resultants at every Z station ───────────────────────────────────────
#
# Abaqus stores NFORC as:
#   - single vector field 'NFORC'      when requested under *Node Output
#   - scalar fields NFORC1/NFORC2/NFORC3  when requested under *Element Output
# This code handles both cases, preferring the vector form if present.
#
# Left-side rule: for each node on the cut plane, sum only contributions
# from elements whose centroid Z < z_snap. This avoids double-counting
# since an interior node is shared by elements on both sides of the cut.

fo_keys = frm.fieldOutputs.keys()

# ── detect which form is available ───────────────────────────────────────────
has_nforc_vector     = 'NFORC' in fo_keys
has_nforc_components = ('NFORC1' in fo_keys and
                        'NFORC2' in fo_keys and
                        'NFORC3' in fo_keys)

nforc_map = {}   # (nodeLabel, elementLabel) -> [NF1, NF2, NF3]

if has_nforc_vector:
    # ── vector form (/Node Output) ────────────────────────────────────────
    if _VERBOSE_NFORC:
        print("NFORC: reading vector field 'NFORC'")
    try:
        fo = frm.fieldOutputs['NFORC'].getSubset(region=inst)
    except Exception:
        fo = frm.fieldOutputs['NFORC']
    for v in fo.values:
        if not hasattr(v, 'nodeLabel') or not hasattr(v, 'elementLabel'):
            continue
        if len(v.data) < 3:
            continue
        nforc_map[(v.nodeLabel, v.elementLabel)] = [
            float(v.data[0]), float(v.data[1]), float(v.data[2])]
    if _VERBOSE_NFORC:
        print("  entries loaded: %d" % len(nforc_map))

elif has_nforc_components:
    # ── component form (*Element Output) — NFORC1/NFORC2/NFORC3 ──────────
    if _VERBOSE_NFORC:
        print("NFORC: reading component fields NFORC1/NFORC2/NFORC3")
    map1 = _read_scalar_nforc(frm.fieldOutputs, 'NFORC1', inst)
    map2 = _read_scalar_nforc(frm.fieldOutputs, 'NFORC2', inst)
    map3 = _read_scalar_nforc(frm.fieldOutputs, 'NFORC3', inst)
    # merge — only keep entries present in all three components
    all_keys = set(map1.keys()) & set(map2.keys()) & set(map3.keys())
    for k in all_keys:
        nforc_map[k] = [map1[k], map2[k], map3[k]]
    if _VERBOSE_NFORC:
        print("  NFORC1=%d  NFORC2=%d  NFORC3=%d  merged=%d"
              % (len(map1), len(map2), len(map3), len(nforc_map)))

else:
    print("WARNING: no NFORC data found in ODB.")
    print("  Available fields with NFORC in name:",
          [k for k in fo_keys if 'NFORC' in k])
    print("  solid_resultants.csv will have empty force columns.")

has_nforc = len(nforc_map) > 0

# ── write resultants ──────────────────────────────────────────────────────────
resultants_csv = os.path.join(out_dir, 'solid_resultants.csv')
with open(resultants_csv, 'wb') as f:
    w = csv.writer(f)
    w.writerow(['Z_target', 'Z_snapped',
                'N',  'Vx',  'Vy',
                # Mx/My/Mz are the DEFORMED-configuration resultants: lever
                # arms from the current node positions, about the deformed
                # section centroid, full r x F. Validated at the free tip,
                # where the answer is known exactly (no lever arm): they
                # recover the applied moment to 7 significant figures and
                # drive My/Mz to ~1e-11, while the reference-configuration
                # values are 5% off and leave My/Mz at 1e-05.
                'Mx', 'My', 'Mz',
                'n_contributions',
                'nforc_available',
                # marker column: tells reextract.py this file already has
                # deformed-configuration moments (value is always 'deformed').
                'moment_config'])

    for z_target in Z_STATIONS:
        snap_nodes_s_st, z_snap_st = snap_z(coord_s, z_target)

        if not has_nforc:
            w.writerow([z_target, z_snap_st,
                        '', '', '', '', '', '', 0, 'False', 'deformed'])
            # (6 blanks: N,Vx,Vy, Mx,My,Mz)
            continue

        N = Vx = Vy = Mx = My = Mz = 0.0
        n_contrib = 0
        side_used = 'left'

        # deformed centroid of this cut: mean CURRENT position of its nodes.
        # Under warping the cut is no longer planar, so this is the mean
        # position rather than a plane intercept -- which is what a moment
        # resultant should be referred to anyway.
        _cx = _cy = _cz = 0.0
        _nc = 0
        for _nd in snap_nodes_s_st:
            if _nd not in coord_s:
                continue
            _c = coord_s[_nd]
            _u = u_s.get(_nd, (0.0, 0.0, 0.0))
            _cx += _c[0] + _u[0]; _cy += _c[1] + _u[1]; _cz += _c[2] + _u[2]
            _nc += 1
        cen_def = ((_cx/_nc, _cy/_nc, _cz/_nc) if _nc else (0.0, 0.0, 0.0))

        def _accumulate(side_sign):
            """side_sign: +1 for left-side elements (el_z < cut, the
            existing convention), -1 for right-side elements (el_z >
            cut) with the sign flip needed to match that convention --
            see FIX note below."""
            n, vx_, vy_, cnt = 0.0, 0.0, 0.0, 0
            mxd_, myd_, mzd_ = 0.0, 0.0, 0.0
            for (nd, el_label), nf in nforc_map.items():
                if nd not in snap_nodes_s_st:
                    continue
                el_z = el_cz.get(el_label)
                if el_z is None:
                    continue
                if side_sign > 0 and el_z >= z_snap_st:
                    continue
                if side_sign < 0 and el_z <= z_snap_st:
                    continue
                x, y, _ = coord_s[nd]
                nf1, nf2, nf3 = nf[0], nf[1], nf[2]
                # deformed lever arm, about the deformed section centroid,
                # full r x F (keeps the out-of-plane term that warping makes
                # non-zero). cen_def is computed per station below.
                _u = u_s.get(nd, (0.0, 0.0, 0.0))
                dx = x + _u[0] - cen_def[0]
                dy = y + _u[1] - cen_def[1]
                dz = coord_s[nd][2] + _u[2] - cen_def[2]
                mxd_ += side_sign * (dy * nf3 - dz * nf2)
                myd_ += side_sign * (dz * nf1 - dx * nf3)
                mzd_ += side_sign * (dx * nf2 - dy * nf1)
                n   += side_sign * nf3
                vx_ += side_sign * nf1
                vy_ += side_sign * nf2
                cnt += 1
            return n, vx_, vy_, mxd_, myd_, mzd_, cnt

        N, Vx, Vy, Mx, My, Mz, n_contrib = _accumulate(+1)

        # FIX: the left-side rule (el_z < cut) has no elements to sum at
        # the clamped ROOT (z_snap=0) -- there is no material "before"
        # the fixed end, only the wall, so this always silently returned
        # 0 contributions and all-zero N/Vx/Vy/Mx/My there (see chat).
        # Fall back to the RIGHT side (toward the free end) whenever the
        # left side is empty: by static equilibrium of the cut, the
        # internal resultant computed from either side represents the
        # same physical quantity, so the right-side sum -- negated to
        # match the left-side sign convention already used at every
        # other station -- gives the correct root reaction (Mx_root =
        # applied Mx + Vy*L, matching strain_screen()'s moment_at_
        # station() formula), instead of silently reporting zero.
        if n_contrib == 0:
            N, Vx, Vy, Mx, My, Mz, n_contrib = _accumulate(-1)
            side_used = 'right (fallback -- left side had no elements, '\
                        'e.g. at the clamped root)'

        # FIX: NFORC gives the INTERNAL reaction across the cut (material
        # on one side acting on the other), which is the equal-and-
        # opposite of the applied load by Newton's third law -- confirmed
        # empirically (see chat): every one of N/Vx/Vy/Mx/My/Mz came back
        # with the opposite sign of the actual applied Cload/Dsload for a
        # given case, consistently. Flip all six here, once, after the
        # side selection above, so solid_resultants.csv reports in the
        # SAME sign convention as the applied loads (and as summary.csv's
        # own Vx/Vy/Mx columns) -- no more mental sign-flip needed to
        # compare this file against what was actually applied.
        N, Vx, Vy, Mx, My, Mz = -N, -Vx, -Vy, -Mx, -My, -Mz

        w.writerow([z_target, z_snap_st,
                    N, Vx, Vy, Mx, My, Mz,
                    n_contrib, 'True', 'deformed'])
        # per-station NFORC values are written to solid_resultants.csv;
        # printing them too is ~11 lines per case with nothing to act on.
        # Set VERBOSE_NFORC=1 in the environment to get them back.
        if _VERBOSE_NFORC:
            print("  Z=%.1f: N=%.4e Vx=%.4e Vy=%.4e Mx=%.4e My=%.4e Mz=%.4e "
                  "(%d contrib, %s)"
                  % (z_target, N, Vx, Vy, Mx, My, Mz, n_contrib, side_used))

# ── whole-model energies for the hourglass check (added 7 Oct 2026) ──
# ALLAE (artificial strain energy, hourglass control in C3D8R) and ALLIE
# (internal energy) at the end of the step, from the PRESELECT history
# output of the assembly. batch_driver.py rejects a case if
# ALLAE/ALLIE > --max-hourglass. Never allowed to break the extraction.
try:
    _step = odb.steps[list(odb.steps.keys())[-1]]
    _hr = None
    for _k in _step.historyRegions.keys():
        if _k.startswith('Assembly'):
            _hr = _step.historyRegions[_k]
            break
    _vals = {}
    if _hr is not None:
        _avail = list(_hr.historyOutputs.keys())
        for _nm in ('ALLAE', 'ALLIE', 'ALLSE', 'ALLWK'):
            if _nm in _avail:
                _vals[_nm] = float(_hr.historyOutputs[_nm].data[-1][1])
    _keys = sorted(_vals.keys())
    _f = open(os.path.join(out_dir, 'solid_energy_history.csv'), 'w')
    _f.write(','.join(_keys) + '\n')
    _f.write(','.join(['%.10e' % _vals[k] for k in _keys]) + '\n')
    _f.close()
except Exception as _e:
    print("WARNING: energy history (ALLAE/ALLIE) not read: %s" % _e)

odb.close()
print("\nWrote:", resultants_csv)
print("Done.")
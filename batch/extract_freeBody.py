# -*- coding: utf-8 -*-
"""
extract_freeBody.py  --  abaqus viewer noGUI=extract_freeBody.py -- <solid.odb> <out_dir> [Z_MID] [Z_STATIONS]

Uses Abaqus/Viewer's free body cut tool to extract section resultants
(forces and moments) at every Z station without manual NFORC integration.
This is more accurate and simpler than integrating NFORC manually.

Run from WSL via:
    abaqus viewer noGUI=extract_freeBody.py -- solid.odb out_dir 50.0 "0,10,20,30,40,50"

Note the -- before the arguments: abaqus viewer requires this to separate
its own arguments from the script's arguments.

Writes: <out_dir>/solid_freeBody.csv
"""
import sys
import os
import csv

# ── parse args (everything after the -- separator) ───────────────────────────
# sys.argv in abaqus viewer noGUI looks like:
#   [scriptname, arg1, arg2, ...]
# The -- is consumed by the abaqus launcher, not passed through.
solid_odb = sys.argv[-4] if len(sys.argv) >= 5 else sys.argv[1]
out_dir   = sys.argv[-3] if len(sys.argv) >= 5 else sys.argv[2]
Z_MID     = float(sys.argv[-2]) if len(sys.argv) >= 5 else 50.0
Z_STATIONS_STR = sys.argv[-1] if len(sys.argv) >= 5 else str(Z_MID)
Z_STATIONS = [float(z) for z in Z_STATIONS_STR.split(',')]

print("solid_odb  :", solid_odb)
print("out_dir    :", out_dir)
print("Z_MID      :", Z_MID)
print("Z_STATIONS :", Z_STATIONS)

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

# ── Abaqus/Viewer imports (only available inside abaqus viewer) ───────────────
from abaqus import session
from abaqusConstants import WHOLE_MODEL

# ── open ODB ──────────────────────────────────────────────────────────────────
print("\nOpening ODB:", solid_odb)
odb = session.openOdb(name=solid_odb)

# use the last step, last frame (same convention as extract_odb.py)
step     = odb.steps.values()[-1]
frame    = step.frames[-1]
instance = odb.rootAssembly.instances.values()[0]

print("Step     :", step.name)
print("Frame    :", frame.frameValue)
print("Instance :", instance.name)

# ── create a viewport (required even in noGUI mode for free body cuts) ────────
vp = session.Viewport(name='vp_freeBody', origin=(0, 0), width=200, height=150)
vp.setValues(displayedObject=odb)

# display the last frame
vp.odbDisplay.setFrame(step=step.name, frame=len(step.frames) - 1)

# ── helper: get node coordinates from ODB ─────────────────────────────────────
coord_s = {n.label: tuple(n.coordinates) for n in instance.nodes}

def snap_z(coord_map, z_target):
    """Return the actual Z value in the mesh closest to z_target."""
    zvals = sorted(set(round(c[2], 8) for c in coord_map.values()))
    return min(zvals, key=lambda z: abs(z - z_target))

# ── free body cut at each Z station ──────────────────────────────────────────
#
# The free body cut works by:
#   1. Defining a cutting plane (normal = Z-axis, passing through z_target)
#   2. Asking Abaqus/Viewer to integrate the stress field over that plane
#   3. Reading back the resultant force and moment vectors
#
# The result is equivalent to NFORC integration but done internally by
# Abaqus using its own Gauss-point data, so it is more accurate.

results = []

for z_target in Z_STATIONS:
    z_snap = snap_z(coord_s, z_target)
    print("\nProcessing Z=%.4f (snapped=%.4f)" % (z_target, z_snap))

    try:
        # ── define the free body cut ─────────────────────────────────────
        # A free body cut needs:
        #   - a name
        #   - the cutting plane: defined by a point on the plane and
        #     the normal direction (Z-axis for a Z=const cut)
        #   - the region: whole model or a named element set
        #
        # modelLinearisation is the Abaqus term for the free body cut
        # result object. It gives forces/moments in the global frame.

        cut_name = 'cut_Z%.4f' % z_snap

        # point on the cutting plane (any point with Z = z_snap)
        # use the section centroid X=0, Y=0 if known; otherwise any point
        cut_point = (0.0, 0.0, z_snap)

        # normal to the cutting plane = beam axis = Z
        cut_normal = (0.0, 0.0, 1.0)

        # create the free body cut
        # session.freeBodies is the container; addFromOdb creates one
        # from the current ODB display state
        free_body = session.FreeBodyFromOdb(
            name=cut_name,
            odb=odb,
            step=step.name,
            frame=len(step.frames) - 1,
            sectionPoint=cut_point,
            normal=cut_normal,
            # summationPoint: where moments are computed about
            # use section centroid (0,0,z_snap) so Mx/My are
            # bending moments about the centroid, not a random point
            summationPoint=cut_point,
            # componentResolution: global frame
            componentResolution=WHOLE_MODEL,
        )

        # ── read the resultants ───────────────────────────────────────────
        # Abaqus returns force and moment as tuples (F1, F2, F3) and
        # (M1, M2, M3) in the global coordinate system.
        # For beam axis = Z:
        #   F3 = N   (axial)
        #   F1 = Vx  (shear X)
        #   F2 = Vy  (shear Y)
        #   M1 = Mx  (bending about X)
        #   M2 = My  (bending about Y)
        #   M3 = Mz  (torque about Z)

        force  = free_body.force    # (F1, F2, F3)
        moment = free_body.moment   # (M1, M2, M3)

        N  = float(force[2])    # axial
        Vx = float(force[0])    # shear X
        Vy = float(force[1])    # shear Y
        Mx = float(moment[0])   # bending about X
        My = float(moment[1])   # bending about Y
        Mz = float(moment[2])   # torque about Z

        print("  N=%.4e  Vx=%.4e  Vy=%.4e  Mx=%.4e  My=%.4e  Mz=%.4e"
              % (N, Vx, Vy, Mx, My, Mz))

        results.append({
            'Z_target':       z_target,
            'Z_snapped':      z_snap,
            'N':              N,
            'Vx':             Vx,
            'Vy':             Vy,
            'Mx':             Mx,
            'My':             My,
            'Mz':             Mz,
            'status':         'ok',
        })

        # clean up the free body object to avoid name collisions on next
        # iteration — Abaqus keeps them in session.freeBodies
        try:
            del session.freeBodies[cut_name]
        except Exception:
            pass

    except Exception as e:
        print("  ERROR at Z=%.4f: %s" % (z_target, e))
        results.append({
            'Z_target':  z_target,
            'Z_snapped': z_snap,
            'N': '', 'Vx': '', 'Vy': '',
            'Mx': '', 'My': '', 'Mz': '',
            'status': 'error: %s' % e,
        })

# ── write CSV ─────────────────────────────────────────────────────────────────
out_csv = os.path.join(out_dir, 'solid_freeBody.csv')
fields  = ['Z_target', 'Z_snapped', 'N', 'Vx', 'Vy', 'Mx', 'My', 'Mz', 'status']

with open(out_csv, 'wb') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for r in results:
        w.writerow(r)

print("\nWrote:", out_csv)

odb.close()
print("Done.")
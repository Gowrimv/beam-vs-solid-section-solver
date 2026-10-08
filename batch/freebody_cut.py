# -*- coding: utf-8 -*-
"""
freebody_cut.py -- STANDALONE.  Not part of the batch pipeline.

    abaqus python freebody_cut.py <solid.odb> [options]

Computes the section force/moment resultant at one or more Z cuts, TWO
independent ways, and prints them side by side so they can be checked
against each other.

ROUTE A -- free-body cut (the real summation)
    Take every NFORC nodal force on the cut that belongs to an element on
    one side of it, and sum  r x f  with r measured in the DEFORMED
    configuration from the DEFORMED centroid of that cut:

        F = sum f_i
        M = sum (x_i + u_i - c_def) x f_i

    This is the only route that sees how the traction is actually
    distributed over the face, so it is the one that stays valid if the
    loading ever stops being a single end load.

ROUTE B -- equilibrium transport (no NFORC at all)
    With no loads applied between the tip and the cut, the internal force
    resultant is constant and the moment simply transports:

        F(z) = F_applied
        M(z) = M_applied + (r_tip_def - r_z_def) x F_applied

    Three numbers per station (the deformed centroid) instead of ~400
    nodal force vectors, and it is exact at the discrete level too,
    because NFORC satisfies nodal equilibrium exactly.

WHY BOTH.  Route B is cheaper and immune to how the traction distributes,
but it makes the moment a DERIVED quantity -- it can only tell you about
the deformed geometry, never about the stress field.  Route A is a genuine
independent measurement.  Run them together once: if they agree, route B
is safe to use for the rest of the sweep.

CAVEAT on route B's anchor.  Fx/Fy are applied as a uniform *Dsload TRVEC
over the tip face, so their resultant acts through that face's DEFORMED
area centroid; this script uses the mean deformed position of the tip-cut
nodes as the estimate.  Mx is a *Cload on a distributing-coupling
reference node -- a pure couple, so its point of application does not
matter at all.  The A-vs-B agreement at the stations FARTHEST from the tip
is the test of the centroid estimate, since that is where the lever arm,
and therefore the sensitivity to it, is largest.

SIGN CONVENTION.  NFORC reports the internal reaction across the cut,
which is equal and opposite to the applied load.  Route A flips all six
components once at the end so both routes report in the same convention as
the applied Cload/Dsload.

OPTIONS
    --stations 0,10,20,...   Z cuts (default: 11 evenly spaced over the span)
    --inp <case_3d.inp>      read Fx, Fy, Mx from the input deck
    --fx --fy --mx           or give them directly (override --inp)
    --area A                 tip-face area, used to turn the TRVEC traction
                             read by --inp into a force (default: config.A)
    --csv <path>             also write the table
"""
from __future__ import print_function
import sys
import os
import csv

from odbAccess import openOdb

# Tip-face area for turning the *Dsload TRVEC traction into a force.
# Follows config.SECTION_SIDE (area = side**2); --area overrides it.
# config.py is PARSED AS TEXT, not imported: this script runs under Abaqus
# Python, which may be 2.7, and config.py uses Python-3-only syntax.
def _area_from_config():
    import re
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.py')
    try:
        txt = open(p).read()
    except Exception:
        return None
    m = re.search(r'^\s*SECTION_SIDE\s*=\s*([0-9.eE+-]+)', txt, re.M)
    return float(m.group(1)) ** 2 if m else None

_DEFAULT_AREA = _area_from_config()
if _DEFAULT_AREA is None:
    print('WARNING: could not read SECTION_SIDE from config.py -- '
          'using area = 1.0; pass --area if the section is not 1 x 1.')
    _DEFAULT_AREA = 1.0


# ---------------------------------------------------------------- arguments
def parse_args(argv):
    o = {'odb': None, 'stations': None, 'inp': None,
         'fx': None, 'fy': None, 'mx': None, 'my': 0.0, 'mz': 0.0, 'csv': None,
         'area': _DEFAULT_AREA}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--stations':
            i += 1; o['stations'] = [float(t) for t in argv[i].split(',')]
        elif a in ('--inp', '--csv'):
            i += 1; o[a[2:]] = argv[i]
        elif a == '--area':
            i += 1; o['area'] = float(argv[i])
        elif a in ('--fx', '--fy', '--mx', '--my', '--mz'):
            i += 1; o[a[2:]] = float(argv[i])
        elif o['odb'] is None:
            o['odb'] = a
        else:
            raise SystemExit('unexpected argument: %s' % a)
        i += 1
    if o['odb'] is None:
        raise SystemExit(__doc__)
    return o


def loads_from_inp(path, area=1.0):
    """Fx/Fy from *Dsload TRVEC (traction magnitude x unit vector x tip-face
    AREA), moments from *Cload dof 4/5/6.  Returns (fx, fy, mx, my, mz).
    area must be the tip-face area (--area); default 1.0 = 1x1 section."""
    f = [0.0, 0.0, 0.0]
    m = [0.0, 0.0, 0.0]
    lines = open(path).read().splitlines()
    for i, ln in enumerate(lines):
        head = ln.strip().lower()
        if i + 1 >= len(lines):
            break
        nxt = [t.strip() for t in lines[i + 1].split(',')]
        if head.startswith('*dsload') and len(nxt) >= 6 and nxt[1].upper() == 'TRVEC':
            mag = float(nxt[2])
            vec = [float(nxt[3]), float(nxt[4]), float(nxt[5])]
            nrm = (vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2) ** 0.5 or 1.0
            for k in range(3):
                f[k] += mag * area * vec[k] / nrm
        elif head.startswith('*cload') and len(nxt) >= 3:
            dof, val = int(nxt[1]), float(nxt[2])
            if dof in (1, 2, 3):
                f[dof - 1] += val
            elif dof in (4, 5, 6):
                m[dof - 4] += val
    return f[0], f[1], m[0], m[1], m[2]


# ------------------------------------------------------------------ helpers
def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def nodal_vec(fo_dict, key, n):
    if key not in fo_dict:
        return {}
    out = {}
    for v in fo_dict[key].values:
        if hasattr(v, 'nodeLabel') and len(v.data) >= n:
            out[v.nodeLabel] = [float(c) for c in v.data[:n]]
    return out


def read_nforc(fo_dict, inst):
    """{(node, element): [f1,f2,f3]} from either a NFORC vector field or
    the NFORC1/NFORC2/NFORC3 scalar fields."""
    if 'NFORC' in fo_dict:
        out = {}
        for v in fo_dict['NFORC'].values:
            if hasattr(v, 'nodeLabel') and hasattr(v, 'elementLabel'):
                out[(v.nodeLabel, v.elementLabel)] = [float(c) for c in v.data[:3]]
        return out
    comp = []
    for k in ('NFORC1', 'NFORC2', 'NFORC3'):
        d = {}
        if k in fo_dict:
            try:
                fo = fo_dict[k].getSubset(region=inst)
            except Exception:
                fo = fo_dict[k]
            for v in fo.values:
                if hasattr(v, 'nodeLabel') and hasattr(v, 'elementLabel'):
                    dd = v.data
                    d[(v.nodeLabel, v.elementLabel)] = \
                        float(dd[0]) if hasattr(dd, '__len__') else float(dd)
        comp.append(d)
    keys = set(comp[0]) & set(comp[1]) & set(comp[2])
    return dict((k, [comp[0][k], comp[1][k], comp[2][k]]) for k in keys)


def snap(coord, z):
    zs = sorted(set(c[2] for c in coord.values()))
    best = min(zs, key=lambda v: abs(v - z))
    diffs = [zs[i + 1] - zs[i] for i in range(len(zs) - 1)]
    tol = min(diffs) * 0.45 if diffs else 0.1
    return set(nd for nd, c in coord.items() if abs(c[2] - best) <= tol), best


# --------------------------------------------------------------------- main
def main():
    o = parse_args(sys.argv[1:])
    if o['inp'] and o['fx'] is None:
        o['fx'], o['fy'], o['mx'], o['my'], o['mz'] = loads_from_inp(o['inp'], o['area'])
    for k in ('fx', 'fy', 'mx'):
        if o[k] is None:
            raise SystemExit('need --inp or --fx/--fy/--mx')
    F_app = (o['fx'], o['fy'], 0.0)
    M_app = (o['mx'], o['my'], o['mz'])

    odb = openOdb(o['odb'], readOnly=True)
    inst = odb.rootAssembly.instances.values()[0]
    frm = odb.steps.values()[-1].frames[-1]
    coord = dict((n.label, tuple(n.coordinates)) for n in inst.nodes)
    u = nodal_vec(frm.fieldOutputs, 'U', 3)
    nforc = read_nforc(frm.fieldOutputs, inst)
    if not nforc:
        raise SystemExit('no NFORC in this odb -- request it under *Element Output')

    el_cz = {}
    for el in inst.elements:
        nl = [n for n in el.connectivity if n in coord]
        if nl:
            el_cz[el.label] = sum(coord[n][2] for n in nl) / float(len(nl))

    zmax = max(c[2] for c in coord.values())
    stations = o['stations'] or [zmax * i / 10.0 for i in range(11)]

    def deformed_centroid(nodes):
        sx = sy = sz = 0.0
        n = 0
        for nd in nodes:
            if nd not in coord:
                continue
            c, du = coord[nd], u.get(nd, [0.0, 0.0, 0.0])
            sx += c[0] + du[0]; sy += c[1] + du[1]; sz += c[2] + du[2]
            n += 1
        return (sx / n, sy / n, sz / n) if n else (0.0, 0.0, 0.0)

    tip_nodes, _ = snap(coord, zmax)
    r_tip = deformed_centroid(tip_nodes)

    print('odb            : %s' % o['odb'])
    print('applied F      : %+.6e %+.6e %+.6e' % F_app)
    print('applied M      : %+.6e %+.6e %+.6e' % M_app)
    print('tip cut centroid (deformed): %+.6f %+.6f %+.6f' % r_tip)
    print('')

    rows = []
    for z in stations:
        nodes, zs = snap(coord, z)
        c_def = deformed_centroid(nodes)

        # ---- ROUTE A: free-body NFORC cut, deformed lever arms ----------
        def accumulate(side):
            Fa = [0.0, 0.0, 0.0]
            Ma = [0.0, 0.0, 0.0]
            cnt = 0
            for (nd, el), f in nforc.items():
                if nd not in nodes:
                    continue
                ez = el_cz.get(el)
                if ez is None or (side > 0 and ez >= zs) or (side < 0 and ez <= zs):
                    continue
                c, du = coord[nd], u.get(nd, [0.0, 0.0, 0.0])
                r = (c[0] + du[0] - c_def[0],
                     c[1] + du[1] - c_def[1],
                     c[2] + du[2] - c_def[2])
                rxf = cross(r, f)
                for k in range(3):
                    Fa[k] += side * f[k]
                    Ma[k] += side * rxf[k]
                cnt += 1
            return Fa, Ma, cnt

        Fa, Ma, cnt = accumulate(+1)
        side_used = 'left'
        if cnt == 0:                      # clamped root: no material below
            Fa, Ma, cnt = accumulate(-1)
            side_used = 'right'
        Fa = [-v for v in Fa]             # internal -> applied convention
        Ma = [-v for v in Ma]

        # ---- ROUTE B: equilibrium transport from the tip ----------------
        arm = (r_tip[0] - c_def[0], r_tip[1] - c_def[1], r_tip[2] - c_def[2])
        axf = cross(arm, F_app)
        Fb = list(F_app)
        Mb = [M_app[k] + axf[k] for k in range(3)]

        rows.append((zs, Fa, Ma, Fb, Mb, cnt, side_used, c_def))

    hdr = ('%6s %4s  %-4s %14s %14s %14s %10s'
           % ('z', 'n', '', 'route A', 'route B', 'A - B', 'rel diff'))
    print(hdr)
    print('-' * len(hdr))
    for zs, Fa, Ma, Fb, Mb, cnt, side, c_def in rows:
        scale = max(abs(v) for v in Ma + Mb) or 1.0
        fscale = max(abs(v) for v in Fa + Fb) or 1.0
        first = True
        for lab, a, b, sc in (('N ', Fa[2], Fb[2], fscale),
                              ('Vx', Fa[0], Fb[0], fscale),
                              ('Vy', Fa[1], Fb[1], fscale),
                              ('Mx', Ma[0], Mb[0], scale),
                              ('My', Ma[1], Mb[1], scale),
                              ('Mz', Ma[2], Mb[2], scale)):
            stem = ('%6.1f %4d  ' % (zs, cnt)) if first else ' ' * 12
            print('%s%-4s %+14.6e %+14.6e %+14.6e %9.4f%%'
                  % (stem, lab, a, b, a - b, 100.0 * (a - b) / sc))
            first = False
        print('')

    # ---- normalisation ---------------------------------------------------
    # One scale per CASE, not per station, so stations are comparable within
    # a case AND cases are comparable with each other.  scale_M is the ROOT
    # moment the case carries, |M_applied| + |F_applied| * span: the largest
    # moment anywhere in the model, so a normalised difference reads directly
    # as "fraction of the biggest thing this case is doing".  Normalising by
    # the LOCAL |M(z)| instead would blow up wherever a component passes
    # through zero -- the same near-zero-denominator trap as relL2_J3.
    scale_F = (F_app[0] ** 2 + F_app[1] ** 2 + F_app[2] ** 2) ** 0.5
    scale_M = (M_app[0] ** 2 + M_app[1] ** 2 + M_app[2] ** 2) ** 0.5 \
              + scale_F * zmax
    scale_F = scale_F or 1.0
    scale_M = scale_M or 1.0
    print('scale_F = %.6e   scale_M = |M_app| + |F_app|*%.1f = %.6e'
          % (scale_F, zmax, scale_M))

    if o['csv']:
        with open(o['csv'], 'wb') as fh:
            w = csv.writer(fh)
            w.writerow(['Z', 'n_contrib', 'side',
                        'cenx_def', 'ceny_def', 'cenz_def',
                        'N_A', 'Vx_A', 'Vy_A', 'Mx_A', 'My_A', 'Mz_A',
                        'N_B', 'Vx_B', 'Vy_B', 'Mx_B', 'My_B', 'Mz_B',
                        'scale_F', 'scale_M',
                        'dN_norm', 'dVx_norm', 'dVy_norm',
                        'dMx_norm', 'dMy_norm', 'dMz_norm',
                        'dF_resultant_norm', 'dM_resultant_norm'])
            for zs, Fa, Ma, Fb, Mb, cnt, side, c_def in rows:
                dF = [Fa[k] - Fb[k] for k in range(3)]
                dM = [Ma[k] - Mb[k] for k in range(3)]
                dFn = sum(v * v for v in dF) ** 0.5 / scale_F
                dMn = sum(v * v for v in dM) ** 0.5 / scale_M
                w.writerow([zs, cnt, side] + list(c_def)
                           + [Fa[2], Fa[0], Fa[1], Ma[0], Ma[1], Ma[2]]
                           + [Fb[2], Fb[0], Fb[1], Mb[0], Mb[1], Mb[2]]
                           + [scale_F, scale_M,
                              dF[2] / scale_F, dF[0] / scale_F, dF[1] / scale_F,
                              dM[0] / scale_M, dM[1] / scale_M, dM[2] / scale_M,
                              dFn, dMn])
        print('wrote %s' % o['csv'])

    odb.close()


if __name__ == '__main__':
    main()

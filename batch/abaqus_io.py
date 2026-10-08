"""
abaqus_io.py
============
Abaqus/WSL job I/O: WSL<->Windows path translation, .inp editing
(write_inp), job submission with stdin starved and a hard timeout
(run_abaqus_case/_run_abaqus_job), and the leftover-output-file cleanup
that keeps Abaqus from ever popping an "Overwrite?" prompt. No FEniCS,
no pyvista -- just subprocess/os/re -- so this module has no Python
package dependencies beyond the stdlib and can be imported without
either of those installed.
"""
import os
import re
import subprocess
import threading

import config
from config import tprint


def win(wsl_path):
    """WSL /mnt/X/... -> Windows X:\\... for abaqus.bat args.
    Non-/mnt/ paths (e.g. /home/...) are returned unchanged."""
    p = os.path.abspath(wsl_path)
    m = re.match(r'^/mnt/([a-zA-Z])(/.*)?$', p)
    if m:
        drive = m.group(1).upper()
        rest  = (m.group(2) or '').replace('/', '\\')
        return drive + ':' + rest
    return p

def abq(*args):
    """Build a cmd.exe command list for any abaqus.bat invocation."""
    return ['cmd.exe', '/c', win(config.ABAQUS_BAT)] + list(args)

RE_1D_VX = re.compile(r'(free,\s*1\s*,\s*)([+-]?\d[\d.eE+-]*)')
RE_1D_VY = re.compile(r'(free,\s*2\s*,\s*)([+-]?\d[\d.eE+-]*)')
RE_1D_MX = re.compile(r'(free,\s*4\s*,\s*)([+-]?\d[\d.eE+-]*)')
RE_3D_VX = re.compile(r'(free,\s*TRVEC\s*,\s*)([+-]?\d[\d.eE+-]*)(,\s*1\.,\s*0\.,\s*0\.)')
RE_3D_VY = re.compile(r'(free,\s*TRVEC\s*,\s*)([+-]?\d[\d.eE+-]*)(,\s*0\.,\s*1\.,\s*0\.)')
RE_3D_MX = re.compile(r'(m_Set-1,\s*4\s*,\s*)([+-]?\d[\d.eE+-]*)')

_RE_SOLID_ELEMENT = re.compile(r'^(\*Element,\s*type=)([A-Za-z0-9]+)', re.I | re.M)
_RE_BRICK8 = re.compile(r'^C3D8[A-Z]*$', re.I)

def set_solid_element(txt, elem_type):
    """Replace the element type on every '*Element, type=...' line of a 3D
    template. Only 8-node bricks are accepted, because the template mesh is
    first-order (8 nodes per element); anything else would give Abaqus the
    wrong number of nodes per element."""
    if not _RE_BRICK8.match(elem_type):
        raise ValueError(
            "ELEMENT_3D=%r not allowed: the 3D template mesh has 8-node "
            "bricks, so only C3D8* types (C3D8R, C3D8, C3D8I, C3D8H, "
            "C3D8RH...) can be swapped in. Quadratic or tetrahedral "
            "elements need a new mesh (and a matching FEniCS section "
            "mesh)." % elem_type)
    txt, n = _RE_SOLID_ELEMENT.subn(
        lambda m: m.group(1) + elem_type.upper(), txt)
    if n == 0:
        raise ValueError("no '*Element, type=' line found in the 3D template")
    return txt

def write_inp(src, dst, vx, vy, mx):
    with open(src, 'r') as f: txt = f.read()
    is_1d = '_1d' in os.path.basename(dst)
    def sub(t, rx, val):
        def repl(m):
            suffix = m.group(3) if m.re.groups >= 3 else ''
            return m.group(1) + ('%.6e' % val) + suffix
        return rx.sub(repl, t, count=1)
    # 1D: *Cload on the tip node -> the value IS the force.
    # 3D: *Dsload TRVEC on the tip face -> the value is a TRACTION
    # (force per unit area), so divide by the section area to apply a
    # total shear force of vx / vy. (With config.A = 1 this is a no-op.)
    if is_1d:
        fvx, fvy = vx, vy
    else:
        fvx, fvy = vx / config.A, vy / config.A
    txt = sub(txt, RE_1D_VX if is_1d else RE_3D_VX, fvx)
    txt = sub(txt, RE_1D_VY if is_1d else RE_3D_VY, fvy)
    txt = sub(txt, RE_1D_MX if is_1d else RE_3D_MX, mx)
    if not is_1d:
        txt = _last_frame_only(txt)
        if getattr(config, 'ELEMENT_3D', None):
            txt = set_solid_element(txt, config.ELEMENT_3D)
    with open(dst, 'w') as f: f.write(txt)

# extract_odb.py only ever reads the LAST frame of the 3D ODB, but the
# template's "*Output, field" has no frequency, so Abaqus writes every
# increment (element-nodal S/E/LE/NE/NFORC... for 10k elements, ~10+ MB
# per increment). "number interval=1" writes the end of the step only.
_RE_FIELD_OUT = re.compile(r'^(\*Output,\s*field)([^\n]*)$', re.I | re.M)

def _last_frame_only(txt):
    def repl(m):
        opts = m.group(2)
        if re.search(r'frequency|number\s*interval|time\s*points', opts, re.I):
            return m.group(0)
        return m.group(1) + opts + ', number interval=1'
    return _RE_FIELD_OUT.sub(repl, txt)

def run(cmd, cwd=None, dry=False):
    """Run any subprocess with stdin starved. Without this, child
    processes inherit this terminal's stdin — if anything in the chain
    (cmd.exe, abaqus python, etc.) ever reads from it, keystrokes or
    buffered lines can get consumed and later replayed against the
    bash prompt, which looks exactly like garbage commands appearing
    on their own. DEVNULL makes any such read return EOF immediately."""
    tprint(' '.join(str(x) for x in cmd))
    if dry: return 0
    return subprocess.run(cmd, cwd=cwd, stdin=subprocess.DEVNULL).returncode

ABAQUS_TIMEOUT_1D = 300    # seconds — B31 jobs finish in seconds; 5 min = clearly hung
ABAQUS_TIMEOUT_3D = 2400   # seconds — finer 20k-element solid: hard NLGEOM cases take up to ~2030 s (was 900, which dropped 60 converged cases)

# Job-scoped output files Abaqus creates. Deleting these before submission
# guarantees no "Overwrite?" prompt ever appears.
_JOB_FILE_EXTS = ['.odb', '.dat', '.msg', '.sta', '.prt', '.lck', '.res',
                  '.mdl', '.stt', '.023', '.log', '.com', '.sim', '.SMABulk']

_SCRATCH_EXTS = ['.sim', '.stt', '.mdl', '.prt', '.res', '.023', '.com',
                 '.SMABulk', '.lck']

def clean_job_files(case_dir, job_name):
    """Delete any leftover output files for job_name so Abaqus never
    prompts to overwrite. Called before every submission."""
    for ext in _JOB_FILE_EXTS:
        p = os.path.join(case_dir, job_name + ext)
        if os.path.exists(p):
            try: os.remove(p)
            except OSError: pass

def _run_abaqus_job(job, inp_wsl, case_dir, timeout):
    """Launch one Abaqus job with stdin starved (so any unexpected
    prompt fails fast instead of hanging) and a hard timeout (so a
    genuinely stuck process gets killed rather than blocking forever).
    Returns True on clean completion, False otherwise."""
    cmd = abq('job='+job, 'input='+win(inp_wsl), 'interactive', 'ask_delete=off')
    try:
        proc = subprocess.Popen(cmd, cwd=case_dir, stdin=subprocess.DEVNULL)
    except Exception as ex:
        tprint('  %s: failed to launch (%s)' % (job, ex))
        return False
    try:
        rc = proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        tprint('  %s: TIMED OUT after %ds — killing process' % (job, timeout))
        proc.kill()
        try: proc.wait(timeout=10)
        except Exception: pass
        return False
    if rc != 0:
        tprint('  %s: exited with code %d' % (job, rc))
        return False
    return True

def job_completed(case_dir, job):
    """True only if Abaqus itself reports a clean finish. abaqus.bat run
    through cmd.exe often returns exit code 0 even when pre.exe or
    standard.exe failed (seen: 'memory allocation request failed' in
    pre.exe left a corrupt .odb and rc=0), so the exit code alone is
    not trusted. Requires the .sta file to contain
    'THE ANALYSIS HAS COMPLETED SUCCESSFULLY'."""
    sta = os.path.join(case_dir, job + '.sta')
    try:
        with open(sta, 'r', errors='ignore') as f:
            return 'COMPLETED SUCCESSFULLY' in f.read().upper()
    except OSError:
        return False

def run_abaqus_case(case_id, job1d, inp1d, job3d, inp3d, case_dir):
    """Run 1D and 3D Abaqus jobs concurrently via threads (not raw Popen
    racing on shared stdin — each job's Popen has its own DEVNULL stdin,
    so they can't interfere with each other or block on a shared console).
    Returns (ok1d, ok3d)."""
    results = {}
    def _worker(name, job, inp, timeout):
        results[name] = _run_abaqus_job(job, inp, case_dir, timeout)

    t1 = threading.Thread(target=_worker,
                          args=('1d', job1d, inp1d, ABAQUS_TIMEOUT_1D))
    t2 = threading.Thread(target=_worker,
                          args=('3d', job3d, inp3d, ABAQUS_TIMEOUT_3D))
    t1.start(); t2.start()
    t1.join(); t2.join()
    ok1 = results.get('1d', False) and job_completed(case_dir, job1d)
    ok3 = results.get('3d', False) and job_completed(case_dir, job3d)
    # Solver-internal files are not needed once the job has finished;
    # keep .odb (data), .sta/.msg/.dat/.log (diagnostics).
    for job, ok in ((job1d, ok1), (job3d, ok3)):
        if ok:
            for ext in _SCRATCH_EXTS:
                p = os.path.join(case_dir, job + ext)
                if os.path.exists(p):
                    try: os.remove(p)
                    except OSError: pass
    for name, job, rc_ok, ok in (('1d', job1d, results.get('1d'), ok1),
                                 ('3d', job3d, results.get('3d'), ok3)):
        if rc_ok and not ok:
            tprint('  %s: exit code 0 but .sta has no "COMPLETED '
                   'SUCCESSFULLY" -- treating as failed' % job)
    return ok1, ok3


# ── 1D convergence check (added 7 Oct 2026) ────────────────────────────────
# Abaqus/Standard judges force equilibrium against a time-averaged force.
# When the loads are tiny in absolute terms (e.g. ~1e-7 with E = 1), that
# average drops below Abaqus's "zero force" level and Abaqus substitutes a
# fallback of 1e-2, so the force check always passes and an UNCONVERGED
# solution is accepted (seen: residual force up to 300% of the applied
# shear). This reads the LAST equilibrium iteration in the .msg file and
# reports the residuals relative to the applied loads.
_RE_RES_F = re.compile(r'LARGEST RESIDUAL FORCE\s+([-+\d.E]+)')
_RE_RES_M = re.compile(r'LARGEST RESIDUAL MOMENT\s+([-+\d.E]+)')
_RE_TAVG_F = re.compile(r'TIME AVG\. FORCE\s+([-+\d.E]+)')

def residual_check_1d(case_dir, job, vx, vy, mx, L, tail_bytes=8000):
    """Return dict with res_force, res_moment, tavg_force, res_ratio_1d and
    tavg_fallback_1d (True if Abaqus used the 1e-2 fallback). res_ratio_1d
    = max(res_force / F_ref, res_moment / M_ref) with
    F_ref = max(|V|, |Mx|/L), M_ref = max(|Mx|, |V| L), so it is
    dimensionless and never divides by ~0. None values if unreadable."""
    out = {'res_force_1d': None, 'res_moment_1d': None, 'tavg_force_1d': None,
           'res_ratio_1d': None, 'tavg_fallback_1d': None}
    p = os.path.join(case_dir, job + '.msg')
    try:
        with open(p, 'rb') as fh:
            fh.seek(0, 2); n = fh.tell(); fh.seek(max(0, n - tail_bytes))
            t = fh.read().decode(errors='ignore')
    except OSError:
        return out
    rf = _RE_RES_F.findall(t); rm = _RE_RES_M.findall(t); ta = _RE_TAVG_F.findall(t)
    try:
        V = (vx*vx + vy*vy) ** 0.5
        F_ref = max(V, abs(mx) / L); M_ref = max(abs(mx), V * L)
        if rf:
            out['res_force_1d'] = abs(float(rf[-1]))
        if rm:
            out['res_moment_1d'] = abs(float(rm[-1]))
        if ta:
            out['tavg_force_1d'] = float(ta[-1])
            out['tavg_fallback_1d'] = abs(out['tavg_force_1d'] - 1e-2) < 1e-9
        ratios = []
        if out['res_force_1d'] is not None and F_ref > 0:
            ratios.append(out['res_force_1d'] / F_ref)
        if out['res_moment_1d'] is not None and M_ref > 0:
            ratios.append(out['res_moment_1d'] / M_ref)
        if ratios:
            out['res_ratio_1d'] = max(ratios)
    except (ValueError, TypeError):
        pass
    return out

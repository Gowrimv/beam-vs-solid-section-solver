#!/usr/bin/env python3
"""
monte_Carlo.py
==============
Standalone: estimates what fraction of random (Vx, Vy, Mx) draws pass the
pre-run strain screen, by brute-force sampling.

Uses the SAME section constants, load ranges and screen as the batch
(config.py + strain_checks.strain_screen), so it always matches whatever
SECTION_SIDE / ranges are currently set -- no private copies.

Run:
    python3 monte_Carlo.py [N]          # default N = 200000
"""
import sys
import numpy as np

import config
from strain_checks import strain_screen, bending_dominance, range_safety_check


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
    rng = np.random.default_rng(0)
    vx = rng.uniform(*config.VX_RANGE, N)
    vy = rng.uniform(*config.VY_RANGE, N)
    mx = rng.uniform(*config.MX_RANGE, N)

    eps = np.vectorize(strain_screen)(vx, vy, mx)
    passed = eps <= config.STRAIN_THRESHOLD
    bd = np.vectorize(bending_dominance)(vx[passed], vy[passed], mx[passed])

    print(f"Section side:       {config.SECTION_SIDE}  (A={config.A:.4g}, "
          f"I={config.I:.4g}, c={config.c:.4g})")
    print(f"Ranges:             Vx,Vy +-{config.VX_RANGE[1]:.3e}   "
          f"Mx +-{config.MX_RANGE[1]:.3e}")
    print(f"N samples:          {N}")
    print(f"Accepted:           {passed.sum()}  ({passed.mean():.2%})")
    print(f"Median strain:      {np.median(eps):.4e}  "
          f"(accepted only: {np.median(eps[passed]):.4e})")
    print(f"Max strain:         {eps.max():.4e}  threshold {config.STRAIN_THRESHOLD:.1e}")
    if passed.any():
        print(f"Bending dominance at Z_MID (accepted, theory): median "
              f"{np.median(bd):.3g}, below {config.MIN_BENDING_DOMINANCE_POST:g}: "
              f"{(bd < config.MIN_BENDING_DOMINANCE_POST).mean():.2%}")
    print()
    range_safety_check(verbose=True)


if __name__ == '__main__':
    main()

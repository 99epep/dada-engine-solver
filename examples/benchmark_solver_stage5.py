"""Serial exact Stage-5 replays with compilation explicitly outside timing."""
import argparse
from pathlib import Path
from benchmark_solver_stage4 import run
from benchmark_solver_acceleration import DEFAULT_MANIFEST

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--manifest', type=Path, default=DEFAULT_MANIFEST)
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--cases', nargs='+')
    p.add_argument('--no-exact-cache', action='store_true')
    p.add_argument('--legacy-replay', action='store_true')
    p.add_argument('--profile', action='store_true')
    a=p.parse_args()
    run(a.output,a.manifest,'numba',a.repeats,a.cases,a.profile,
        exact_cache=not a.no_exact_cache,warm_compile=True,shared_replay=not a.legacy_replay)

"""Serial frozen replays of the explicit Numba wall-RHS research prototype."""
import argparse
import hashlib
import json
from pathlib import Path
import time

STARTED=time.perf_counter()
from benchmark_solver_acceleration import DEFAULT_MANIFEST, run
BASE_IMPORT_SECONDS=time.perf_counter()-STARTED


def benchmark(output,manifest,cases,repeats,backend):
    if output.exists(): raise FileExistsError('Use a new measurement directory.')
    started=time.perf_counter()
    stats={}
    if backend=='numba':
        before=time.perf_counter()
        import numba
        from numba_wall_prototype import prototype_backend
        import_seconds=time.perf_counter()-before
        with prototype_backend() as stats:
            run(manifest,output,cases,repeats)
        stats.update(numba_version=numba.__version__,backend_import_seconds=import_seconds)
    else:
        run(manifest,output,cases,repeats)
    source=Path(__file__).with_name('numba_wall_prototype.py')
    stats.update(backend=backend,base_import_seconds=BASE_IMPORT_SECONDS,
        run_including_backend_import_seconds=time.perf_counter()-started,
        prototype_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        settings=dict(fastmath=False,parallel=False,disk_cache=False),
        note='First supported RHS call includes lazy compilation. Reports and convergence remain Python.')
    (output/'prototype.json').write_text(json.dumps(stats,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,default=DEFAULT_MANIFEST)
    parser.add_argument('--cases',nargs='+')
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--backend',choices=('python','numba'),default='numba')
    args=parser.parse_args()
    benchmark(args.output,args.manifest,args.cases,args.repeats,args.backend)

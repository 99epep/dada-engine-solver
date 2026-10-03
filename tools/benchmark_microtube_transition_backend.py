"""Bounded, fixed-state RHS comparison; never runs or edits a campaign.

Use a raw history record, its compatible study and a saved trajectory NPZ with
`angles` and `trajectory` arrays. An optional pre-change numerical_primitives.py
measures the former fallback-heavy kernel against identical states. No timing
assertions: scheduling, CPU and JIT startup vary between machines.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import inspect
import json
import statistics
import sys
import time

import numpy as np

from dada_solver.research.schema import compile_study, load_study
from dada_solver.wall_backend import WallBackendSettings, WallRHS, CompiledWallRHS


def reference_dispatcher(path):
    import numba
    from numba.extending import register_jitable
    name = '_dada_benchmark_reference_numeric'
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    for name, function in vars(module).items():
        if inspect.isfunction(function) and function.__module__ == module.__name__ and name != 'wall_kernel':
            register_jitable(function)
    return numba.njit(fastmath=False, parallel=False)(module.wall_kernel)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--record', type=Path, required=True, help='One raw campaign history record as JSON')
    parser.add_argument('--trajectory', type=Path, required=True)
    parser.add_argument('--reference-source', type=Path)
    parser.add_argument('--samples', type=int, default=128)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.samples < 1 or args.repeats < 1: parser.error('Samples/repeats must be positive.')
    if args.output.exists(): parser.error('Output already exists; choose a new benchmark artifact.')
    record = json.loads(args.record.read_text())
    definition = compile_study(load_study(args.study))
    start = time.perf_counter()
    wrapper = definition.adapter.build(dict(definition.fixed_parameters, **record['physical'])).build()
    construction = time.perf_counter()-start
    with np.load(args.trajectory) as trajectory:
        indices = np.linspace(0, len(trajectory['angles'])-1, args.samples, dtype=int)
        samples = [(float(trajectory['angles'][i]), trajectory['trajectory'][:, i].copy()) for i in indices]
    reference = np.array([wrapper.derivative(a, v) for a, v in samples])
    implementations = [('python', WallRHS(wrapper, WallBackendSettings('python'))),
                       ('numba_continuum', WallRHS(wrapper, WallBackendSettings('numba')))]
    if args.reference_source:
        old = WallRHS(wrapper, WallBackendSettings('numba'))
        if not isinstance(old.implementation, CompiledWallRHS):
            raise ValueError('Benchmark requires a supported compiled model.')
        old.implementation.dispatcher = reference_dispatcher(args.reference_source)
        implementations.append(('reference_kernel', old))
    results = {}
    for label, rhs in implementations:
        start = time.perf_counter(); rhs(*samples[0]); first = time.perf_counter()-start
        # Warm every exercised path before steady timing, and check identical inputs.
        actual = np.array([rhs(a, v) for a, v in samples])
        error = np.abs(actual-reference)
        scale = 1e-11+2e-11*np.abs(reference)
        np.testing.assert_allclose(actual, reference, rtol=2e-11, atol=1e-11)
        times = []; before = rhs.snapshot()
        for _ in range(args.repeats):
            start = time.perf_counter()
            for a, v in samples: rhs(a, v)
            times.append(time.perf_counter()-start)
        snapshot = rhs.snapshot()
        results[label] = dict(first_call_seconds=first, repeat_seconds=times,
            median_seconds=statistics.median(times), calls_per_repeat=len(samples),
            maximum_absolute_rhs_error=float(error.max()), maximum_scaled_rhs_error=float((error/scale).max()),
            measured_calls=snapshot['calls']-before['calls'],
            measured_fallback_calls=snapshot['fallback_calls']-before['fallback_calls'],
            backend=snapshot)
        if label == 'reference_kernel':
            # The wrapper identity is current: name the actual old kernel separately.
            results[label]['reference_kernel_sha256'] = hashlib.sha256(args.reference_source.read_bytes()).hexdigest()
    result = dict(candidate_id=record['candidate_id'], scope='same saved states; no integration or optimization',
        input_sha256={name:hashlib.sha256(getattr(args,name).read_bytes()).hexdigest()
                      for name in ('study','record','trajectory')},
        construction_seconds=construction, results=results)
    if 'reference_kernel' in results:
        result['reference_over_continuum_ratio'] = results['reference_kernel']['median_seconds']/results['numba_continuum']['median_seconds']
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({name:{key:value for key,value in item.items() if key in
        ('first_call_seconds','median_seconds','measured_fallback_calls','maximum_scaled_rhs_error')}
        for name,item in results.items()},indent=2))


if __name__ == '__main__':
    main()

"""Local mass-conserving periodic-map study; never changes production iteration.

Central differences in nine scaled tangent coordinates are repeated at several
steps. Quadratures reset each cycle and are not shooting coordinates. The exact
last-angle cache experiment is candidate-local and stores no thermodynamic data.
"""
from dataclasses import replace
import argparse
import json
from pathlib import Path
import time
import numpy as np
from benchmark_solver_acceleration import ROOT, DEFAULT_MANIFEST, digest, environment
from refine_motor_four_stage_hx9d_variable_gas import _load_basis, _build_design
from dada_solver.exchangers.wall_cycle import WallCycleNumericalSettings
from dada_solver.wall_backend import WallRHS, WallBackendSettings

from dada_solver.exchangers.wall_iteration import tangent_coordinates, MASS_INDICES as MASS
LABELS = ['m_S', 'U_S', 'm_L', 'U_L', 'm_Hi', 'U_Hi', 'm_Ho', 'U_Ho', 'E_wall_Hi', 'E_wall_Ho']


def central_jacobian(function, dimension, step):
    columns = []
    for j in range(dimension):
        delta = np.zeros(dimension)
        delta[j] = step
        columns.append((function(delta)-function(-delta))/(2*step))
    return np.column_stack(columns)


class ExactAngleProvider:
    """Observe exact repeats, optionally reuse only the most recent angle."""
    def __init__(self, forward, enabled):
        self.forward, self.enabled = forward, enabled
        self.last_angle, self.last_value = None, None
        self.calls = self.repeats = 0
        self.angles = set()

    def __getattr__(self, name):
        return getattr(self.forward, name)

    def cylinder_volumes_and_derivatives(self, angle):
        self.calls += 1
        self.angles.add(angle)
        if angle == self.last_angle:
            self.repeats += 1
            if self.enabled:
                return self.last_value
        self.last_angle = angle
        self.last_value = self.forward.cylinder_volumes_and_derivatives(angle)
        return self.last_value


def integrate(wrapper, rhs, state, settings):
    return wrapper.integrate_cycle(state, rhs=rhs,
        integration_method=settings.integration_method,
        rtol=settings.integration_relative_tolerance,
        atol=np.array(settings.integration_absolute_tolerances),
        maximum_step_angle=settings.maximum_step_angle_radians)


def angle_experiment(wrapper, anchor, settings, repeats=3):
    rows = []
    reference = None
    # Alternate order to reduce systematic warm-up/order bias.
    for repeat in range(repeats):
        for enabled in ((False, True) if repeat % 2 == 0 else (True, False)):
            provider = ExactAngleProvider(wrapper.model.kinematics, enabled)
            candidate = replace(wrapper, model=replace(wrapper.model, kinematics=provider))
            rhs = WallRHS(candidate, WallBackendSettings('numba'))
            start = time.perf_counter()
            angles, trajectory = integrate(candidate, rhs, anchor, settings)
            elapsed = time.perf_counter()-start
            if reference is None:
                reference = angles, trajectory
            np.testing.assert_array_equal(angles, reference[0])
            np.testing.assert_array_equal(trajectory, reference[1])
            rows.append(dict(repeat=repeat, cache=enabled, seconds=elapsed,
                calls=provider.calls, consecutive_exact_repeats=provider.repeats,
                distinct_angles=len(provider.angles), bitwise_equal=True))
    return rows


def spectrum(matrix, basis):
    eigenvalues, vectors = np.linalg.eig(matrix)
    order = np.argsort(-np.abs(eigenvalues))
    result = []
    for j in order:
        mode = basis @ vectors[:, j]
        # Normalize phase and amplitude for readable signed relative-state modes.
        pivot = np.argmax(np.abs(mode))
        mode = mode / mode[pivot]
        temperatures = np.r_[mode[1:8:2]-mode[0:8:2], mode[8:10]]
        result.append(dict(real=float(eigenvalues[j].real), imag=float(eigenvalues[j].imag),
            magnitude=float(abs(eigenvalues[j])),
            relative_state_mode_real=mode.real.tolist(), relative_state_mode_imag=mode.imag.tolist(),
            relative_temperature_mode_real=temperatures.real.tolist(),
            relative_temperature_mode_imag=temperatures.imag.tolist()))
    return result


def run(output, steps=(1e-3, 3e-4, 1e-4), ode_factor=1., cases=None):
    if not steps or any(not np.isfinite(h) or h <= 0 for h in steps):
        raise ValueError('Finite positive difference steps are required.')
    if not np.isfinite(ode_factor) or ode_factor <= 0:
        raise ValueError('The ODE tolerance factor must be finite and positive.')
    if cases and set(cases)-{'production_cold', 'smooth_four_bar'}:
        raise ValueError('Only the two saved tighter-periodic anchors are supported.')
    if output.exists():
        raise FileExistsError('Use a new diagnostic output directory.')
    output.mkdir(parents=True)
    manifest = json.loads(DEFAULT_MANIFEST.read_text())
    for name, expected in manifest['input_hashes'].items():
        if digest(ROOT/name) != expected:
            raise ValueError('Changed frozen input: '+name)
    raw = dict(manifest['numerical_settings'])
    raw['integration_absolute_tolerances'] = tuple(ode_factor*np.array(raw['integration_absolute_tolerances']))
    raw['integration_relative_tolerance'] *= ode_factor
    settings = WallCycleNumericalSettings(**raw)
    _, base = _load_basis()
    metadata = dict(environment=environment(), manifest_sha256=digest(DEFAULT_MANIFEST),
        source_hashes={str(p.relative_to(ROOT)): digest(p) for p in sorted((ROOT/'src').rglob('*.py'))},
        script_sha256=digest(Path(__file__)), steps=steps, ode_tolerance_factor=ode_factor,
        state_labels=LABELS, coordinate_definition='x = anchor + diag(abs(anchor)) @ basis @ z')
    (output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    for case in manifest['cases']:
        if case['name'] not in (cases or ('production_cold', 'smooth_four_bar')):
            continue
        design = _build_design(base, case['parameters'])
        if case.get('motion') == 'base':
            design = replace(design, kinematics=base.kinematics)
        wrapper = design.build()
        anchor_path = ROOT/'outputs/solver_acceleration_stage4/numba_tight'/(case['name']+'_0.json')
        anchor = np.array(json.loads(anchor_path.read_text())['final_state'])
        scales, basis = tangent_coordinates(anchor)
        rhs = WallRHS(wrapper, WallBackendSettings('numba'))
        rhs(0., np.r_[anchor, np.zeros(5)])  # Exclude JIT from cycle measurements.
        cycle_records = []
        def physical_map(state):
            start = time.perf_counter()
            try:
                angles, trajectory = integrate(wrapper, rhs, state, settings)
            except (ValueError, RuntimeError, ArithmeticError) as error:
                cycle_records.append(dict(seconds=time.perf_counter()-start,
                    exception=type(error).__name__, message=str(error),
                    relative_initial_perturbation=((state-anchor)/scales).tolist()))
                raise
            cycle_records.append(dict(seconds=time.perf_counter()-start,
                relative_mass_drift=float(np.max(np.abs(trajectory[MASS].sum(axis=0)-state[MASS].sum()))/state[MASS].sum()),
                minimum_state=float(trajectory[:10].min()), samples=len(angles)))
            return trajectory[:10, -1]
        endpoint = physical_map(anchor)
        closure = float(np.max(np.abs(endpoint-anchor)/(settings.periodic_absolute_tolerance+
            settings.periodic_relative_tolerance*np.maximum(np.abs(endpoint), np.abs(anchor)))))
        def reduced_map(z):
            return basis.T @ ((physical_map(anchor+scales*(basis@z))-anchor)/scales)
        differences = []
        previous = None
        for step in steps:
            try:
                matrix = central_jacobian(reduced_map, 9, step)
            except (ValueError, RuntimeError, ArithmeticError) as error:
                differences.append(dict(step=step, exception=type(error).__name__, message=str(error)))
                print(case['name'], 'step', step, type(error).__name__, str(error), flush=True)
                continue
            item = dict(step=step, matrix=matrix.tolist(), spectrum=spectrum(matrix, basis),
                newton_matrix_condition=float(np.linalg.cond(matrix-np.eye(9))),
                operator_norm=float(np.linalg.norm(matrix, 2)),
                change_from_previous_step=None if previous is None else float(np.linalg.norm(matrix-previous, 2)))
            differences.append(item)
            previous = matrix
            print(case['name'], 'step', step, 'dominant eigenvalue', item['spectrum'][0]['real'],
                'imag', item['spectrum'][0]['imag'], flush=True)
        angle_rows = angle_experiment(wrapper, anchor, settings)
        result = dict(case=case['name'], anchor_source=str(anchor_path.relative_to(ROOT)),
            anchor_sha256=digest(anchor_path), anchor=anchor.tolist(), scales=scales.tolist(), basis=basis.tolist(),
            baseline_normalized_periodic_error=closure, baseline_endpoint=endpoint.tolist(),
            differences=differences, cycles=cycle_records, backend=rhs.snapshot(),
            exact_angle_experiment=angle_rows)
        (output/(case['name']+'.json')).write_text(json.dumps(result, indent=2)+'\n')
        print(case['name'], 'complete;', len(cycle_records), 'map cycles;', 'closure', closure, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--steps', nargs='+', type=float, default=[1e-3, 3e-4, 1e-4])
    parser.add_argument('--ode-factor', type=float, default=1.)
    parser.add_argument('--cases', nargs='+')
    args = parser.parse_args()
    run(args.output, args.steps, args.ode_factor, args.cases)

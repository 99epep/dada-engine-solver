"""Derived report curves: one saved-state cycle, never a campaign reevaluation."""
from dataclasses import replace
import gzip
import json
import math
from pathlib import Path
import time
import numpy as np

from dada_solver.campaign.candidate import content_hash
from dada_solver.campaign.definition import runtime_identity
from .schema import compile_study

DEFAULT_PLOTS = ('positions', 'heat', 'pressures', 'temperatures')
PLOT_NAMES = ('positions', 'volumes', 'pressures', 'heat', 'temperatures', 'flows',
              'velocity', 'acceleration', 'mechanisms')
THERMAL_PLOTS = frozenset(('pressures', 'heat', 'temperatures', 'flows'))


def plot_names(value):
    if value is None: return ()
    values = value.split(',') if isinstance(value, str) else [name for group in value for name in group.split(',')]
    result = []
    for item in values:
        for name in item.split(','):
            if name == 'none': continue
            for expanded in PLOT_NAMES if name == 'all' else DEFAULT_PLOTS if name == 'default' else (name,):
                if expanded not in PLOT_NAMES: raise ValueError(f'Unknown plot {expanded!r}; choose from {PLOT_NAMES}, default, all or none.')
                if expanded not in result: result.append(expanded)
    if 'none' in values and len(values) != 1: raise ValueError('Use plots none on its own.')
    return tuple(result)


def curves(angle, unit, series, *, method, angle_domain='solver_cycle_angle'):
    return dict(angle=np.degrees(angle).tolist(), angle_unit='deg', angle_domain=angle_domain,
                unit=unit, series=[dict(label=k, values=np.asarray(v).tolist()) for k,v in series.items()], method=method)


def kinematic_plots(study, record, names):
    """Use the production operation transform exactly once for piston curves."""
    result = {}
    requested = set(names) & {'positions','volumes','velocity','acceleration'}
    if requested:
        d = compile_study(study)
        built = d.adapter.build(dict(d.fixed_parameters, **record['physical'])).build()
        model = getattr(built, 'model', built)
        angles = np.linspace(0, 2*math.pi, 361)
        volume = {}; position = {}; velocity = {}; acceleration = {}
        for side in ('small','large'):
            limits = getattr(model.machine_volumes, side+'_cylinder')
            law = model.kinematics
            v = np.array([getattr(law, side+'_cylinder_volume')(float(t)) for t in angles])
            volume[side] = v
            position[side] = (v-limits.minimum)/limits.swept
            if 'velocity' in names:
                velocity[side] = [getattr(law, side+'_cylinder_volume_derivative')(float(t))/limits.swept for t in angles]
            if 'acceleration' in names:
                fn = getattr(law, side+'_cylinder_volume_second_derivative', None)
                if fn is not None:
                    try: acceleration[side] = [fn(float(t))/limits.swept for t in angles]
                    except NotImplementedError: pass
        for name,unit,series in [('positions','fraction of swept volume',position),('velocity','1/rad',velocity),('acceleration','1/rad²',acceleration)]:
            if name in names:
                result[name] = curves(angles,unit,series,method='Current production kinematics; no thermodynamic integration.')
                if not series: result[name]['unavailable'] = 'This family does not provide second derivatives.'
        if 'volumes' in names:
            from .visualization import sample_report_volumes
            result['volumes'] = sample_report_volumes(study,record['physical'])
    if 'mechanisms' in names:
        from .visualization import sample_motion
        result['mechanisms'] = sample_motion(study,record['physical'],samples=181)
    return result


def replay_thermal(study, record, *, notify=lambda message: None):
    """One physical cycle from the saved endpoint, without periodic reconvergence.

    Curves are new derived evidence under the current runtime. Stored metrics,
    classification, retry history and candidate identity are never overwritten.
    """
    saved = record.get('final_periodic_state')
    if not saved or not saved.get('periodic_solution'):
        raise ValueError('No saved final periodic state; thermodynamic curves are unavailable.')
    definition = compile_study(study)
    design = definition.adapter.build(dict(definition.fixed_parameters, **record['physical']))
    built = design.build(); model = getattr(built,'model',built)
    state = np.asarray(saved['values'],float)
    started = last_notice = time.monotonic()
    deadline = started+definition.candidate_budget_seconds
    def check(*args):
        nonlocal last_notice
        now = time.monotonic()
        if now-last_notice>=10:
            notify(f"Thermodynamic replay {record['candidate_id'][:12]}: {now-started:.0f}s elapsed.")
            last_notice = now
        if now > deadline:
            from dada_solver.integration import IntegrationInterrupted
            raise IntegrationInterrupted('Report replay budget exhausted.')
    notify(f"Thermodynamic replay {record['candidate_id'][:12]}: one cycle from the saved periodic state (no optimization).")
    from dada_solver.exchangers.external_stream import ExternalStreamWallMachine
    if isinstance(built,ExternalStreamWallMachine):
        if state.shape != (10,): raise ValueError('Expected a saved ten-state wall cycle.')
        from dada_solver.exchangers.wall_cycle import solve_periodic_wall_machine,WallDiagnosticCycle
        from dada_solver.diagnostic_replay import replay_wall_trajectory
        replay = solve_periodic_wall_machine(built,state,maximum_cycles=1,progress_callback=check,
            settings=replace(definition.wall_numerical_settings,accelerate_walls=False),backend=definition.wall_backend)
        if replay.trajectory is None: raise ValueError(f'Report cycle unavailable: {replay.message}')
        cycle = WallDiagnosticCycle(replay.angles,replay.trajectory)
        facts = replay_wall_trajectory(built,cycle,replay.trajectory)
        points = [s.point for s in facts.samples]
        heats = {}; temperatures = {}
        for i,side in enumerate(('heat_in','heat_out')):
            heats[side+' gas → wall'] = [-s.walls[i].gas_heat_w for s in facts.samples]
            temperatures[side+' wall'] = [s.walls[i].wall_temperature_k for s in facts.samples]
        endpoint = replay.trajectory[:10,-1]
        backend = replay.backend_statistics
    else:
        if state.shape != (8,): raise ValueError('Expected a saved eight-state reservoir cycle.')
        from dada_solver.factory import build_periodic_solver,initial_valve_topology
        from dada_solver.state import ThermodynamicState
        from dada_solver.dynamics import ValveTopology
        from dada_solver.valves import ValveState
        topology = saved.get('topology')
        topology = ValveTopology(**{k:ValveState(v) for k,v in topology.items()}) if topology else initial_valve_topology()
        cycle = build_periodic_solver(design.configuration,model,check).integrator.integrate_cycle(ThermodynamicState.from_array(state),topology)
        if not cycle.completed: raise ValueError(f'Report cycle unavailable: {cycle.message}')
        points = [model.instantaneous_point(float(a),ThermodynamicState.from_array(x),top)
            for a,x,top in zip(cycle.angles,cycle.states.T,cycle.topologies)]
        heats = {'heat_in gas → reservoir':[-model.cold_heat_transfer.heat_rate(p.temperatures[2]) for p in points],
                 'heat_out gas → reservoir':[-model.hot_heat_transfer.heat_rate(p.temperatures[3]) for p in points]}
        temperatures = {}; endpoint = cycle.states[:,-1]; backend = {'actual_backend':'python'}
    # Keep actual integrator samples, including valve transitions. No curve fit.
    angles = cycle.angles
    labels = ('small gas','large gas','heat_in gas','heat_out gas')
    pressures = {name:[p.pressures[i] for p in points] for i,name in enumerate(labels)}
    temperatures.update({name:[p.temperatures[i] for p in points] for i,name in enumerate(labels)})
    flows = {name:[getattr(p.flows,name) for p in points] for name in ('small_to_cold','cold_to_large','large_to_hot','hot_to_small')}
    method = 'One-cycle replay from saved periodic endpoint under current runtime; historical verdict and metrics unchanged.'
    result = {name:curves(angles,unit,series,method=method) for name,unit,series in
              [('pressures','Pa',pressures),('temperatures','K',temperatures),('flows','kg/s',flows),('heat','W',heats)]}
    metadata = dict(method=method,backend=backend,endpoint_relative_drift=float(np.max(np.abs(endpoint-state)/np.maximum(np.abs(state),1e-30))),
                    sample_count=len(angles),thermodynamic_replay=True)
    for plot in result.values(): plot['replay'] = metadata
    from dada_solver.valves import ValveState
    result['machine_cycle'] = dict(angle_rad=np.asarray(angles).tolist(), series={
        'S_temperature_K': temperatures['small gas'],
        'L_temperature_K': temperatures['large gas'],
        'Hi_temperature_K': temperatures['heat_in gas'],
        'Ho_temperature_K': temperatures['heat_out gas'],
        'Ho_to_S_valve_open': [int(p.topology.hot_to_small is ValveState.OPEN) for p in points],
        'Hi_to_L_valve_open': [int(p.topology.cold_to_large is ValveState.OPEN) for p in points],
    })
    return result


def candidate_plots(study, record, names, *, cache_directory=None, notify=lambda message: None):
    result = kinematic_plots(study,record,tuple(name for name in names if name != 'mechanisms'))
    mechanism_source = None
    if 'mechanisms' in names:
        from .machine_render import candidate_mechanisms
        mechanism_source = candidate_mechanisms(study,record['physical'])
        result['mechanisms'] = dict(unavailable='A complete physical mechanism is required for both cylinders.')
    if not THERMAL_PLOTS.intersection(names) and mechanism_source is None: return result
    key = content_hash(dict(version=3,study_id=study.study_id,candidate_id=record['candidate_id'],
        physical=record['physical'],state=record.get('final_periodic_state'),runtime=runtime_identity()))
    path = Path(cache_directory)/(key+'.json.gz') if cache_directory is not None else None
    try:
        cached = None
        if path is not None and path.exists():
            try:
                with gzip.open(path,'rt') as stream: cached = json.load(stream)
            except (OSError,ValueError,EOFError): cached = None
        if cached is not None and cached.get('key') == key:
            thermal = cached['plots']; notify(f"Thermodynamic replay {record['candidate_id'][:12]}: using matching derived cache.")
        else:
            thermal = replay_thermal(study,record,notify=notify)
            if path is not None:
                # Cache is disposable and separate from scientific campaign history.
                import os,tempfile
                path.parent.mkdir(parents=True,exist_ok=True)
                fd,tmp = tempfile.mkstemp(dir=path.parent,suffix='.tmp')
                try:
                    with os.fdopen(fd,'wb') as stream:
                        stream.write(gzip.compress(json.dumps(dict(key=key,plots=thermal),allow_nan=False).encode(),mtime=0))
                    os.replace(tmp,path)
                finally:
                    if os.path.exists(tmp): os.unlink(tmp)
        result.update({name:thermal[name] for name in names if name in THERMAL_PLOTS})
        if mechanism_source is not None:
            try:
                from .machine_render import render_machine_webp, layout_for_configuration
                import base64
                artifacts, configuration = mechanism_source
                animation, metadata = render_machine_webp(artifacts,thermal['machine_cycle'],
                    layout=layout_for_configuration(configuration),motor_operation=configuration.motor_operation)
                result['mechanisms'] = dict(metadata, media_type='image/webp',
                    data_uri='data:image/webp;base64,'+base64.b64encode(animation).decode('ascii'),
                    size_bytes=len(animation), replay=thermal.get('temperatures',{}).get('replay'))
            except (ValueError,RuntimeError,ArithmeticError,ImportError) as error:
                message = f'{type(error).__name__}: {error}'
                notify(f"Machine animation {record['candidate_id'][:12]} unavailable: {message}")
                result['mechanisms'] = dict(unavailable=message)
    except (ValueError,RuntimeError,ArithmeticError,ImportError) as error:
        message = f'{type(error).__name__}: {error}'
        notify(f"Thermodynamic replay {record['candidate_id'][:12]} unavailable: {message}")
        for name in names:
            if name in THERMAL_PLOTS or (name=='mechanisms' and mechanism_source is not None): result[name] = dict(unavailable=message)
    return result

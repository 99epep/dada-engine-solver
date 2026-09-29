"""Best-effort snapshots of existing locals at a microtube rejection.

Only invoked after failure. No closure is evaluated again: unavailable quantities
remain null, especially when a hydraulic solve failed before producing a flow.
The traceback is inspected immediately and is never retained in result artifacts.
"""
import math


def _snapshot(error):
    """Extract the first rejected RHS trial state, not an accepted cycle sample."""
    result = dict(criterion=str(error), state_kind='rejected_trial_state',
        exchanger=None, passage=None, angle_rad=None, time_s=None,
        mass_flow_kg_s=None, reynolds=None, prandtl=None, mach=None, knudsen=None,
        p1_pa=None, p2_pa=None, temperature_k=None, inner_diameter_m=None,
        tube_length_m=None, tube_count=None)
    frames = []
    trace = error.__traceback__
    while trace is not None:
        frames.append((trace.tb_frame.f_globals.get('__name__'),
                       trace.tb_frame.f_code.co_name, trace.tb_frame.f_locals))
        trace = trace.tb_next
    model = link = None
    for module, name, values in frames:
        if module == 'dada_solver.dynamics' and name == 'instantaneous_point':
            model = values['self']
            result['angle_rad'] = values['theta']
        if module == 'dada_solver.exchangers.external_stream' and name == 'derivative':
            model = values['self'].model
            result['angle_rad'] = values['angle']
        if module == 'dada_solver.exchangers.external_stream' and name == '<genexpr>':
            if 'index' in values:
                result['exchanger'] = {2: 'H_i', 3: 'H_o'}.get(values['index'])
        if module == 'dada_solver.exchangers.hardware' and name == '_gas_flow':
            link = values['self']
            _geometry(result, link.bank)
            result.update(p1_pa=values['pin'], p2_pa=values['pout'],
                          temperature_k=values['temperature'])
            # A pre-root estimate is not the flow of the rejected hydraulic state.
            if 'diagnostics' in values:
                result['mass_flow_kg_s'] = values['flow']
                _numbers(result, values['diagnostics'])
            elif 're' in values:
                result['hydraulic_trial_reynolds'] = values['re']
        if module == 'dada_solver.exchangers.gas_film' and name == 'evaluate':
            _geometry(result, values['self'].bank)
            result.update(mass_flow_kg_s=values.get('flow'), p1_pa=values.get('p1'),
                          p2_pa=values.get('p2'), temperature_k=values['temperature'])
            if 'd' in values:
                _numbers(result, values['d'])
            # Passages are ordered by the production flow-context constructor.
            result['passage_index'] = values.get('passage_index')
    if model is not None:
        if result['angle_rad'] is not None:
            result['time_s'] = result['angle_rad']/model.angular_speed
        if link is not None:
            for candidate, exchanger, passage in (
                (model.small_cold_link, 'H_i', 'small_to_cold'),
                (model.cold_large_valve.flow_model, 'H_i', 'cold_to_large'),
                (model.large_hot_link, 'H_o', 'large_to_hot'),
                (model.hot_small_valve.flow_model, 'H_o', 'hot_to_small')):
                if candidate is link:
                    result.update(exchanger=exchanger, passage=passage)
                    result['flow_orientation'] = 'upstream_to_downstream'
                    break
    if result['exchanger'] and result.get('passage_index') in (0, 1):
        result['passage'] = {'H_i': ('small_to_cold', 'cold_to_large'),
                             'H_o': ('large_to_hot', 'hot_to_small')}[result['exchanger']][result['passage_index']]
        result['flow_orientation'] = 'signed_named_passage'
    return {key: (None if isinstance(value, float) and not math.isfinite(value) else value)
            for key, value in result.items()}


def _geometry(result, bank):
    result.update(inner_diameter_m=bank.inner_diameter_m,
                  tube_length_m=bank.tube_length_m, tube_count=bank.tube_count)


def _numbers(result, diagnostics):
    for name in ('reynolds', 'prandtl', 'mach', 'knudsen'):
        result[name] = getattr(diagnostics, name)
    result['issues'] = list(diagnostics.issues)
    result['relative_pressure_drop'] = diagnostics.compressibility_parameter
    result['correlation_id'] = diagnostics.correlation_id


def microtube_failure_snapshot(error):
    """Instrumentation must never replace the original scientific failure."""
    try:
        return _snapshot(error)
    except Exception as diagnostic_error:
        return dict(criterion=str(error), snapshot_unavailable=type(diagnostic_error).__name__)

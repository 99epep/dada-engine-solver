"""Best-effort snapshots of existing locals at a microtube rejection.

Only invoked after failure. No closure is evaluated again: unavailable quantities
remain null, especially when a hydraulic solve failed before producing a flow.
The traceback is inspected immediately and is never retained in result artifacts.
"""
import math


def microtube_failure_snapshot(error):
    """Extract the first rejected RHS trial state, not an accepted cycle sample."""
    result = dict(criterion=str(error), state_kind='rejected_trial_state',
        exchanger=None, passage=None, angle_rad=None, time_s=None,
        mass_flow_kg_s=None, reynolds=None, prandtl=None, mach=None, knudsen=None,
        p1_pa=None, p2_pa=None, temperature_k=None, inner_diameter_m=None,
        tube_length_m=None, tube_count=None, pressure_ratio=None,
        relative_pressure_drop=None, passage_index=None, flow_orientation=None,
        hydraulic_trial_reynolds=None, issues=None, correlation_id=None)
    frames = []
    trace = error.__traceback__
    while trace is not None:
        frames.append((trace.tb_frame.f_globals.get('__name__'),
                       trace.tb_frame.f_code.co_name, trace.tb_frame.f_locals))
        trace = trace.tb_next
    model = link = None
    for module, name, values in frames:
        if module == 'dada_solver.dynamics' and name == 'instantaneous_point':
            model = values.get('self')
            result['angle_rad'] = values.get('theta')
        if module == 'dada_solver.exchangers.external_stream' and name == 'derivative':
            model = getattr(values.get('self'), 'model', model)
            result['angle_rad'] = values.get('angle', result['angle_rad'])
        if module == 'dada_solver.exchangers.external_stream' and name == '<genexpr>':
            if 'index' in values:
                result['exchanger'] = {2: 'heat_in', 3: 'heat_out'}.get(values.get('index'))
        if module == 'dada_solver.exchangers.hardware' and name == '_gas_flow':
            link = values.get('self')
            _geometry(result, getattr(link, 'bank', None))
            result.update(p1_pa=values.get('pin'), p2_pa=values.get('pout'),
                          temperature_k=values.get('temperature'))
            # A pre-root estimate is not the flow of the rejected hydraulic state.
            if 'diagnostics' in values:
                result['mass_flow_kg_s'] = values.get('flow')
                _numbers(result, values.get('diagnostics'))
            elif 're' in values:
                result['hydraulic_trial_reynolds'] = values.get('re')
        if module == 'dada_solver.exchangers.gas_film' and name == 'evaluate':
            _geometry(result, getattr(values.get('self'), 'bank', None))
            result.update(mass_flow_kg_s=values.get('flow'), p1_pa=values.get('p1'),
                          p2_pa=values.get('p2'), temperature_k=values.get('temperature'))
            # Passages are ordered by the production flow-context constructor.
            result['passage_index'] = values.get('passage_index')
        if module == 'dada_solver.exchangers.gas_correlations' and name == 'require':
            _numbers(result, values.get('diagnostics'))
    if model is not None:
        if result['angle_rad'] is not None and getattr(model, 'angular_speed', 0):
            result['time_s'] = result['angle_rad']/model.angular_speed
        if link is not None:
            for candidate, exchanger, passage in (
                (getattr(model, 'small_cold_link', None), 'heat_in', 'small_to_cold'),
                (getattr(getattr(model, 'cold_large_valve', None), 'flow_model', None), 'heat_in', 'cold_to_large'),
                (getattr(model, 'large_hot_link', None), 'heat_out', 'large_to_hot'),
                (getattr(getattr(model, 'hot_small_valve', None), 'flow_model', None), 'heat_out', 'hot_to_small')):
                if candidate is link:
                    result.update(exchanger=exchanger, passage=passage)
                    result['flow_orientation'] = 'upstream_to_downstream'
                    break
    if result['exchanger'] in ('heat_in', 'heat_out') and result.get('passage_index') in (0, 1):
        result['passage'] = {'heat_in': ('small_to_cold', 'cold_to_large'),
                             'heat_out': ('large_to_hot', 'hot_to_small')}[result['exchanger']][result['passage_index']]
        result['flow_orientation'] = 'signed_named_passage'
    return {key: (None if isinstance(value, float) and not math.isfinite(value) else value)
            for key, value in result.items()}


def _geometry(result, bank):
    result.update(inner_diameter_m=getattr(bank, 'inner_diameter_m', None),
                  tube_length_m=getattr(bank, 'tube_length_m', None), tube_count=getattr(bank, 'tube_count', None))


def _numbers(result, diagnostics):
    for name in ('reynolds', 'prandtl', 'mach', 'knudsen'):
        result[name] = getattr(diagnostics, name, None)
    issues = getattr(diagnostics, 'issues', None)
    result['issues'] = list(issues) if issues is not None else None
    result['relative_pressure_drop'] = getattr(diagnostics, 'compressibility_parameter', None)
    result['correlation_id'] = getattr(diagnostics, 'correlation_id', None)

    result['pressure_ratio'] = getattr(diagnostics, 'pressure_ratio', None)

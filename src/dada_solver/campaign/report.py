"""Factual phase summaries and deterministic suggestions for human review."""
from collections import Counter
import numpy as np


def elite_records(records, size=5):
    unique = {}
    for record in records:
        objective = record.get('objective')
        if (record['status'] == 'feasible' and objective and objective['available']
                and objective['value'] is not None and np.isfinite(objective['value'])):
            unique.setdefault(record['candidate_id'], record)
    return sorted(unique.values(), key=lambda r:(r['objective']['value'], r['candidate_id']))[:size]


def summary(record):
    if record is None: return None
    return {k:record.get(k) for k in ('candidate_id','objective','physical','normalized','metrics','constraints')}


def time_statistics(records):
    durations = [r['duration_seconds'] for r in records if not r.get('cache_hit')]
    return dict(count=len(durations), total_seconds=sum(durations),
        median_seconds=float(np.median(durations)) if durations else None,
        p90_seconds=float(np.percentile(durations, 90)) if durations else None,
        maximum_seconds=max(durations) if durations else None)


def make_report(space, before, phase, *, requested_seconds, elapsed_seconds, elite_size=5):
    start = elite_records(before, 1); all_records = before+phase
    end = elite_records(all_records, elite_size)
    start_best, end_best = (start[0] if start else None), (end[0] if end else None)
    start_coords = start_best['normalized'] if start_best else space.initial_coordinates
    start_physical = start_best['physical'] if start_best else space.decode(start_coords)
    changes = []
    if end_best:
        for i, p in enumerate(space.parameters):
            changes.append(dict(name=p.name, physical_start=start_physical[p.name],
                physical_end=end_best['physical'][p.name],
                normalized_delta=end_best['normalized'][i]-start_coords[i]))
        changes.sort(key=lambda x:(-abs(x['normalized_delta']), x['name']))
    unique_phase = [r for r in phase if not r.get('cache_hit')]
    bound_pressure = []
    for i, p in enumerate(space.parameters):
        for side in ('lower','upper'):
            def near(r): return r['normalized'][i] <= .05 if side == 'lower' else r['normalized'][i] >= .95
            count = sum(near(r) for r in unique_phase)
            elite_count = sum(near(r) for r in end)
            if count >= 2 or elite_count >= 2:
                bound_pressure.append(dict(parameter=p.name, side=side, threshold=.05,
                    phase_count=count, phase_total=len(unique_phase), elite_count=elite_count, elite_total=len(end)))
    failures = Counter(r['status'] for r in unique_phase if r['status'] != 'feasible')
    reasons = Counter(r.get('reason') or r['status'] for r in unique_phase if r['status'] != 'feasible')
    violations = Counter(c['name'] for r in unique_phase for c in r.get('constraints', []) if not c['available'] or not c['satisfied'])
    gain = start_best['objective']['value']-end_best['objective']['value'] if start_best and end_best else None
    suggestions = []
    if unique_phase and sum(r['status']=='feasible' for r in unique_phase)/len(unique_phase) < .2:
        dominant = violations.most_common(1) or reasons.most_common(1)
        suggestions.append('Few feasible candidates: inspect '+str(dominant)+' before changing bounds or assumptions.')
    for item in bound_pressure:
        if len(end) >= 3 and item['elite_count']/len(end) >= .6 and gain is not None and gain > 0:
            suggestions.append(f"Objective improved and {item['elite_count']}/{len(end)} elites are near the {item['side']} bound of {item['parameter']}; consider reviewing that bound.")
    if len(end) >= 3:
        points = np.array([r['normalized'] for r in end])
        diameter = float(np.max(np.abs(points[:,None,:]-points[None,:,:])))
        if diameter < .15:
            suggestions.append(f'Elite normalized diameter is {diameter:.3f}; consider a narrower region or a robustness test.')
        values = [r['objective']['value'] for r in end]
        if diameter > .4 and max(values)-min(values) <= .05*max(abs(min(values)), 1e-12):
            suggestions.append('Distant feasible elites have objective values within 5%; preserve both regions rather than collapsing the search.')
    if len(unique_phase) >= 5 and gain is not None and gain <= 1e-3*max(abs(start_best['objective']['value']),1e-12) and end_best:
        margins = end_best['constraints']
        if margins and all(c['available'] and c['margin'] is not None and c['margin']>0 for c in margins):
            suggestions.append('At least five new evaluations gave less than 0.1% improvement and the best listed margins are positive; consider local refinement after reviewing their physical sizes.')
    return dict(requested_duration_seconds=requested_seconds, actual_duration_seconds=elapsed_seconds,
        attempted=len(phase), cache_hits=sum(bool(r.get('cache_hit')) for r in phase),
        rejected_before_integration=sum(not r['integrated'] for r in unique_phase),
        integrated=sum(r['integrated'] for r in unique_phase),
        converged=sum(r['converged'] for r in unique_phase), feasible=sum(r['status']=='feasible' for r in unique_phase),
        best_at_phase_start=summary(start_best), best_at_phase_end=summary(end_best),
        best_in_phase=[summary(r) for r in elite_records(unique_phase, elite_size)],
        objective_improvement=gain, improvement_convention='start minus end; objectives are minimized',
        parameter_change_reference='phase_start_best' if start_best else 'configured_initial_values_not_evaluated',
        parameter_changes=changes, bound_pressure=bound_pressure,
        status_counts=dict(failures), dominant_reasons=dict(reasons.most_common()),
        violated_constraints=dict(violations.most_common()),
        top_distinct_feasible=[summary(r) for r in end], evaluation_time=time_statistics(unique_phase),
        suggestions=suggestions or ['No further deterministic suggestion is supported by this phase.'])


def readable_report(report):
    def brief(best):
        return 'none' if best is None else f"{best['candidate_id']} objective={best['objective']['value']}"
    lines = ['Optimization campaign phase',
        f"Requested: {report['requested_duration_seconds']:.1f} s; actual: {report['actual_duration_seconds']:.1f} s",
        f"Attempts: {report['attempted']}; cached: {report['cache_hits']}; preflight rejected: {report['rejected_before_integration']}; integrated: {report['integrated']}; converged: {report['converged']}; feasible: {report['feasible']}",
        f"Objective improvement: {report['objective_improvement']} (minimization)",
        f"Best at start: {brief(report['best_at_phase_start'])}",
        f"Best at end: {brief(report['best_at_phase_end'])}"]
    best = report['best_at_phase_end']
    if best:
        lines.append(f"Best physical metrics: {best['metrics']}")
        lines.append('Best constraint margins:')
        for c in best['constraints']:
            lines.append(f"  {c['name']}: margin={c['margin']}; available={c['available']}; satisfied={c['satisfied']}")
    lines.append('Largest normalized parameter changes:')
    for x in report['parameter_changes'][:8]:
        lines.append(f"  {x['name']}: {x['physical_start']:.6g} -> {x['physical_end']:.6g}; normalized delta {x['normalized_delta']:+.4f}")
    lines.append('Top distinct feasible candidates:')
    lines.extend('  '+brief(r) for r in report['top_distinct_feasible'])
    lines.extend([f"Bound pressure: {report['bound_pressure']}", f"Failure reasons: {report['dominant_reasons']}",
        f"Violated/unavailable constraints: {report['violated_constraints']}", f"Evaluation timing: {report['evaluation_time']}",
        'Suggestions (no automatic changes):', *['  '+x for x in report['suggestions']]])
    return '\n'.join(lines)+'\n'

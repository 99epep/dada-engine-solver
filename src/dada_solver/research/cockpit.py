"""Deterministic presentation evidence and review commands; never alter a study."""
from collections import Counter
import re
import shlex
from dada_solver.campaign.report import elite_records


def unique_prefixes(ids, minimum=8):
    ordered = sorted(set(ids)); result = {}
    for i, candidate in enumerate(ordered):
        neighbors = ordered[max(0,i-1):i] + ordered[i+1:i+2]
        size = minimum
        while any(other.startswith(candidate[:size]) for other in neighbors): size += 1
        result[candidate] = candidate[:size]
    return result


def reason_category(record, scientific):
    """Use declared domains, never round temperatures into guessed thresholds."""
    transport_failure=(record.get('diagnostics') or {}).get('transport_failure',{})
    if transport_failure.get('category') in ('transport_temperature_below_domain','transport_temperature_above_domain','transport_temperature_nonfinite'):
        return transport_failure['category']
    reason = record.get('reason') or ''
    match = re.search(r'Transport temperature ([\deE+.\-]+) K outside declared domain', reason)
    if match:
        temperature = float(match[1]); domains = []
        for side in ('heat_in','heat_out'):
            inputs = (scientific.get('basis',{}).get(side) or {}).get('inputs',{})
            transport = (inputs.get('gas_model') or {}).get('transport',{})
            if 'minimum_temperature' in transport and 'maximum_temperature' in transport:
                domains.append((transport['minimum_temperature'], transport['maximum_temperature']))
        if domains and temperature < min(x[0] for x in domains): return 'transport_temperature_below_domain'
        if domains and temperature > max(x[1] for x in domains): return 'transport_temperature_above_domain'
        return 'transport_temperature_outside_domain'
    lower = reason.lower()
    for fragment, category in (
        ('hydrodynamic_entry_unresolved','hydrodynamic_entry_domain'),
        ('large_relative_pressure_drop','relative_pressure_drop_domain'),
        ('transition flow','hydraulic_transition_domain'),('high_mach','mach_domain'),
        ('beyond_continuum','continuum_domain'),('knudsen','continuum_domain'),
        ('thermal closure','thermal_correlation_domain'),('slip','slip_model_domain'),
        ('outside table','fluid_table_domain')):
        if fragment in lower: return category
    return record['status']


def campaign_evidence(records, scientific, source, is_campaign):
    unique = {r['candidate_id']:r for r in records}
    sampled = list(unique.values()); elites = elite_records(sampled,5)
    active = [p for p in scientific['parameters'] if 'value' not in p]
    bounds = []
    for i,p in enumerate(active):
        if p.get('kind')=='choice': continue
        for side in ('lower','upper'):
            def near(r):
                values=r.get('normalized',[])
                return len(values)>i and (values[i]<=.05 if side=='lower' else values[i]>=.95)
            count=sum(near(r) for r in sampled); elite_count=sum(near(r) for r in elites)
            bounds.append(dict(parameter=p['name'], side=side, limit=p[side], unit=p.get('unit','1'),
                count=count, total=len(sampled), percent=100*count/len(sampled) if sampled else 0.,
                elite_count=elite_count, elite_total=len(elites),
                pressed=bool(len(elites)>=3 and elite_count/len(elites)>=.6)))
    failures = Counter(reason_category(r,scientific) for r in records if r['status']!='feasible')
    violations = Counter(c['name'] for r in records for c in r.get('constraints',[]) if c.get('available') and not c.get('satisfied'))
    report = ['dada-research','report',source]
    report_command = shlex.join(report)
    suggestions = []
    pressed = [b for b in bounds if b['pressed']]
    if pressed:
        for b in pressed:
            suggestions.append(dict(signal=f"{b['elite_count']}/{b['elite_total']} elites are within 5% normalized distance of the {b['side']} bound of {b['parameter']}. Review this bound before extending exploration; no bound change is inferred.", action=dict(type='bounds',parameter=b['parameter'],side=b['side'])))
    elif is_campaign and (len(sampled)<512 or len(elites)<5):
        suggestions.append(dict(signal=f'{len(sampled)} distinct candidates and {len(elites)} retained feasible elites: below the display reference of 512 candidates / 5 elites. Continuing the same bounds is an available action, not a convergence or saturation claim.',
            command=shlex.join(['dada-research','resume',source,'--budget','30m'])))
    for category,count in sorted(failures.items(), key=lambda x:(-x[1],x[0]))[:3]:
        suggestions.append(dict(signal=f'{count} attempts classified as {category}. Inspect exact failure evidence before changing the study or model domain.',action=dict(type='filter',field='failure_category',value=category,source=source)))
    for name,count in sorted(violations.items(), key=lambda x:(-x[1],x[0]))[:3]:
        suggestions.append(dict(signal=f'{name} is violated in {count} records with available evidence. Inspect its values and margins; limits are not relaxed automatically.',action=dict(type='filter',field='violated_constraint',value=name,source=source)))
    suggestions.append(dict(signal='Regenerate this derived report from the stored journal; no integration.',command=report_command))
    return dict(source=source,is_campaign=is_campaign,
        funnel=dict(attempted=len(records),integrated=sum(bool(r.get('integrated')) for r in records),
            converged=sum(bool(r.get('converged')) for r in records),feasible=sum(r['status']=='feasible' for r in records),
            cache_hits=sum(bool(r.get('cache_hit')) for r in records),distinct_candidates=len(sampled)),
        rejection_categories=dict(sorted(failures.items(),key=lambda x:(-x[1],x[0]))),
        violated_constraints=dict(violations),bounds=bounds,suggestions=suggestions,
        bound_convention='All distinct candidate coordinates, including rejected candidates; latest attempt per ID. Near = outer 5% of normalized search interval (including log transforms). Elites = up to five distinct feasible objective-ranked candidates. Pressed = at least 3 elites and at least 60% near a bound.')

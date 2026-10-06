"""Inspectable constraint values; reporting never invents a scientific recommendation."""
import math
from .families import MECHANICAL_METRICS, available_metrics


def margin_record(name, value, limit, relation, unit, *, method=None):
    available=value is not None and (not isinstance(value,(int,float)) or math.isfinite(value))
    margin=None
    if available:
        margin=value-limit if relation=='minimum' else limit-value if relation=='maximum' else -abs(float(value)-float(limit))
    relative=margin/abs(limit) if margin is not None and isinstance(limit,(int,float)) and limit!=0 and relation!='equal' else None
    return dict(name=name,value=value if available else None,limit=limit,relation=relation,unit=unit,
        margin=margin,relative_margin=relative,satisfied=bool(available and margin>=0),available=available,
        state='unavailable' if not available else 'satisfied' if margin>=0 else 'violated',
        near_active=bool(available and relation!='equal' and relative is not None and 0<=relative<=.05),method=method)


def validate_mechanical_constraint(row, family, *, scoped=True):
    expected={'metric','relation','limit','unit'}|({'side'} if scoped else set())
    if set(row)!=expected: raise ValueError('Mechanical constraint requires metric, relation, limit, unit and study side.')
    metric=row['metric']
    if metric not in available_metrics(family): raise ValueError(f'Metric {metric!r} is unavailable for {family}.')
    unit,relation=MECHANICAL_METRICS[metric]
    if row['unit']!=unit or row['relation'] not in ('minimum','maximum','equal') or relation is not None and row['relation']!=relation:
        raise ValueError(f'Incorrect unit or relation for mechanical metric {metric}.')
    limit=row['limit']
    if isinstance(limit,bool) or not isinstance(limit,(int,float)) or not math.isfinite(limit) or limit<0:
        raise ValueError('Mechanical limit must be finite and nonnegative.')
    if metric=='zero_crossing_count' and type(limit) is not int: raise ValueError('Reversal count must be an integer.')


PHYSICAL_CONSTRAINTS={
    'minimum_motor_power':('required_power','minimum','W'),
    'minimum_cooling_power':('required_power','minimum','W'),
    'maximum_mechanical_input_power':('limit','maximum','W'),
    'maximum_pressure':('limit','maximum','Pa'),
    'maximum_temperature':('limit','maximum','K'),
    'maximum_absolute_mass_flow':('limit','maximum','kg/s'),
    'maximum_mach_number':('limit','maximum','1'),
    'valid_thermodynamic_model':(None,'equal','1'),
    'periodic_convergence':(None,'equal','1')}


def enrich_constraints(records, declarations):
    """Existing margins are authoritative; restore value/limit without a solver replay."""
    declared={r['type']:r for r in declarations}
    out=[]
    for row in records:
        if 'value' in row:
            out.append(dict(row)); continue
        definition=declared.get(row['name'],{})
        field,relation,unit=PHYSICAL_CONSTRAINTS.get(row['name'],(None,'equal','1'))
        limit=definition.get(field) if field else True
        available=row.get('available',False)
        margin=row.get('margin')
        if field and limit is not None and available and margin is not None:
            value=limit+margin if relation=='minimum' else limit-margin
            entry=margin_record(row['name'],value,limit,relation,unit)
        elif not field:
            entry=margin_record(row['name'],row['satisfied'] if available else None,True,'equal','1',method='categorical verdict')
            # Preserve categorical +/-1 margin convention, never pretend SI margin.
            entry['margin']=margin
        else:
            entry=margin_record(row['name'],None,limit,relation,unit)
        entry.update(satisfied=row['satisfied'],available=available)
        out.append(entry)
    return out


def limiting_evidence(record, scientific):
    """Presentation-only limits from stored verdicts and declared model boundaries.

    No universal mass-flow limit, inferred missing threshold, or new feasibility
    verdict is introduced. Reject-state evidence is not a periodic-cycle margin.
    """
    rows=[dict(c,scope='configured constraint / stored verdict',context=record['status'])
          for c in record.get('constraints',[])]
    def boundary(name,value,limit,context,method):
        if limit is None:
            row=dict(name=name,value=value,limit=None,relation='maximum',unit='1',margin=None,
                relative_margin=None,available=False,satisfied=False,state='unavailable',near_active=False,method=method)
        else: row=margin_record(name,value,limit,'maximum','1',method=method)
        rows.append(dict(row,scope='model-domain evidence',context=context))
    basis=scientific.get('basis',{})
    validity=(record.get('metrics') or {}).get('validity') or {}
    thresholds=basis.get('configuration',{}).get('validity',{})
    for metric in ('maximum_cp_variation','maximum_compressibility_deviation'):
        if validity.get(metric) is not None:
            boundary('thermodynamic.'+metric,validity[metric],thresholds.get(metric),'cycle model validity','stored validity diagnostic / declared threshold')
    passages=((record.get('derived') or {}).get('microtube_gas_domains') or {}).get('passages',{})
    for passage,values in passages.items():
        side='heat_in' if passage.startswith('Hi.') else 'heat_out' if passage.startswith('Ho.') else None
        if side is None: continue
        limits=basis.get(side,{}).get('inputs',{}).get('gas_model') or {}
        for metric,declared in (('mach','maximum_mach'),('compressibility_parameter','maximum_relative_pressure_drop')):
            value=(values.get('ranges',{}).get(metric) or {}).get('maximum')
            if value is not None:
                boundary('microtube.'+metric,value,limits.get(declared),passage,'sampled cycle maximum; not a continuous guarantee')
        hydraulic_mach=(values.get('hydraulic_upstream_ranges',{}).get('mach') or {}).get('maximum')
        if hydraulic_mach is not None:
            boundary('microtube.hydraulic_mach',hydraulic_mach,limits.get('maximum_mach'),passage,'sampled hydraulic upstream maximum; not a continuous guarantee')
    failure=(record.get('diagnostics') or {}).get('first_microtube_failure')
    # A retained first trial is history, not a boundary of a recovered cycle.
    # A retry can also fail for a different reason; do not attribute that final
    # outcome to the first snapshot. Partial rejection records may omit reason.
    reason=record.get('reason') or ''
    same_cause=bool(failure) and (not reason or any(
        issue.strip() in reason for issue in (failure.get('criterion') or '').split(';') if issue.strip()))
    if failure and record['status']=='invalid_exchanger' and not record.get('converged') and same_cause:
        side=failure.get('exchanger')
        limits=basis.get(side,{}).get('inputs',{}).get('gas_model') or {}
        context=f"{side or 'unknown exchanger'} / {failure.get('passage') or 'unknown passage'}; angle={failure.get('angle_rad')}; time={failure.get('time_s')}"
        for metric,declared in (('mach','maximum_mach'),('relative_pressure_drop','maximum_relative_pressure_drop')):
            if failure.get(metric) is not None:
                boundary('microtube.'+metric,failure[metric],limits.get(declared),context,'first rejected trial state')
        rows.append(dict(name=failure.get('criterion') or 'microtube rejection',value=None,limit=None,
            relation=None,unit=None,margin=None,relative_margin=None,available=True,satisfied=False,
            state='violated',near_active=False,scope='recorded rejection',context=context,
            method='exact criterion; unavailable quantitative boundaries are not inferred'))
    def priority(row):
        severity=0 if row.get('available') and not row.get('satisfied') else 1 if not row.get('available') else 2 if row.get('near_active') else 3
        relative=row.get('relative_margin')
        return severity,abs(relative) if relative is not None else math.inf,row['name'],row.get('context','')
    return sorted(rows,key=priority)

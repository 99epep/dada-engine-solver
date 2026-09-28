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
            # Preserve historical categorical +/-1 convention, never pretend SI margin.
            entry['margin']=margin
        else:
            entry=margin_record(row['name'],None,limit,relation,unit)
        entry.update(satisfied=row['satisfied'],available=available)
        out.append(entry)
    return out

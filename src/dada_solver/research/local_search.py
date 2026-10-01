"""Validation and offline evidence for local Sobol regions."""
import math
from dada_solver.campaign.parameters import ChoiceParameter
from collections import Counter
from dada_solver.campaign.report import elite_records
from .schema import keys


def validate_search(search, space):
    common=('type','seed','scramble','domain')
    local=search.get('domain')=='local_regions_v1'
    keys(search,common,'search',('radius_fraction','allocation','evaluate_centers','regions') if local else ('evaluate_initial',))
    if search['type']!='sobol' or search['domain'] not in ('fixed_global_bounds','local_regions_v1') or type(search['seed']) is not int or search['seed']<0 or type(search['scramble']) is not bool:
        raise ValueError('Use explicit bounded Sobol settings.')
    if not local:
        if type(search.get('evaluate_initial',False)) is not bool: raise ValueError('evaluate_initial must be boolean.')
        return
    if set(search)!=set(common)|{'radius_fraction','allocation','evaluate_centers','regions'}:
        raise ValueError('Local regions require radius_fraction, allocation, evaluate_centers and regions.')
    radius=search['radius_fraction']
    if isinstance(radius,bool) or not isinstance(radius,(int,float)) or not math.isfinite(radius) or not 0<radius<=1:
        raise ValueError('Local radius_fraction must be finite and in (0, 1].')
    if search['allocation']!='round_robin' or search['evaluate_centers'] is not True:
        raise ValueError('Local regions require round_robin allocation and evaluate_centers = true.')
    if not space.parameters: raise ValueError('Local refinement requires active ordered parameters.')
    if not isinstance(search['regions'],list) or not search['regions']: raise ValueError('Provide at least one local region.')
    ids=set()
    for region in search['regions']:
        keys(region,('id','source_candidate_id','source_study_id','center'),'search.regions')
        if not isinstance(region['id'],str) or not region['id'] or region['id'] in ids: raise ValueError('Region IDs must be nonempty and unique.')
        ids.add(region['id'])
        for key in ('source_candidate_id','source_study_id'):
            value=region[key]
            if not isinstance(value,str) or len(value)!=64 or any(c not in '0123456789abcdef' for c in value): raise ValueError('Region provenance requires SHA-256 IDs.')
        if not isinstance(region['center'],dict):
            raise ValueError('Region centers require named physical values.')
        for parameter in space.parameters:
            value=region['center'].get(parameter.name)
            if not isinstance(parameter,ChoiceParameter) and (isinstance(value,bool) or not isinstance(value,(int,float))):
                raise ValueError('Numeric region coordinates require finite physical numeric values.')
        space.encode(region['center'])


def region_summary(records, scientific):
    search=scientific.get('search',{})
    if search.get('domain')!='local_regions_v1': return None
    rows=[]
    for region in search['regions']:
        selected=[r for r in records if region['id'] in r.get('search_origin',{}).get('region_ids',
            [r.get('search_origin',{}).get('region_id')])]
        from .cockpit import reason_category
        failures=Counter(reason_category(r,scientific) for r in selected if r['status']!='feasible')
        best=elite_records(selected,1)
        rows.append(dict(**region,attempts=len(selected),converged=sum(bool(r.get('converged')) for r in selected),
            feasible=sum(r['status']=='feasible' for r in selected),
            best=best[0] if best else None,
            rejection_counts=dict(failures),rejection_rates={k:v/len(selected) for k,v in failures.items()}))
    return dict(radius_fraction=search['radius_fraction'],allocation='round_robin',regions=rows)

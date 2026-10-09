"""Packaged design preferences promoted to explicit synthesis constraints."""
from importlib.resources import files
from dataclasses import dataclass, asdict
import math
import tomllib

FINAL_MECHANICAL_SAMPLES = 1440
SCREEN_FILES = {'six_bar_design': 'six_bar_design.toml', 'four_bar_design': 'four_bar_design.toml'}


def mechanical_screen(name):
    if name in (None, 'none'):
        return dict(name='none', version=None, constraints=[])
    if name not in SCREEN_FILES:
        raise ValueError(f'Unknown mechanical screen {name!r}.')
    return tomllib.loads(files('dada_solver.research').joinpath('data/'+SCREEN_FILES[name]).read_text())


def overlay_profile_constraints(defaults, overrides):
    """Replace new-policy defaults by identity, never inherited requirements."""
    result = [dict(row) for row in defaults]
    for row in overrides:
        key = (row['metric'], row['relation'])
        matches = [i for i, old in enumerate(result) if (old['metric'], old['relation']) == key]
        if len(matches) > 1: raise ValueError(f'Ambiguous profile constraint: {key!r}.')
        if matches:
            if result[matches[0]]['unit'] != row['unit']:
                raise ValueError(f'Incompatible units for {key!r}.')
            result[matches[0]] = dict(row)
        else: result.append(dict(row))
    validate_constraint_intervals(result)
    return tuple(result)


def validate_constraint_intervals(rows):
    """Reject incompatible units and empty declared mechanical intervals."""
    for row in rows:
        value=row['limit']
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
            raise ValueError(f"Mechanical limit for {row['metric']} must be finite and nonnegative.")
    for metric in {r['metric'] for r in rows}:
        group = [r for r in rows if r['metric'] == metric]
        if len({r['unit'] for r in group}) != 1:
            raise ValueError(f'Incompatible units for {metric}.')
        low = max((r['limit'] for r in group if r['relation'] in ('minimum','equal')), default=-math.inf)
        high = min((r['limit'] for r in group if r['relation'] in ('maximum','equal')), default=math.inf)
        if low > high: raise ValueError(f'Contradictory mechanical constraints for {metric}.')


def strongest_constraints(*groups):
    """Consolidate commitments without weakening any parent constraint."""
    result = []
    for rows in groups:
        for row in rows:
            old = next((r for r in result if (r['metric'],r['relation']) == (row['metric'],row['relation'])), None)
            if old is None: result.append(dict(row)); continue
            if old['unit'] != row['unit']: raise ValueError(f"Incompatible units for {row['metric']}.")
            if row['relation'] == 'minimum': old['limit'] = max(old['limit'],row['limit'])
            elif row['relation'] == 'maximum': old['limit'] = min(old['limit'],row['limit'])
            elif old['limit'] != row['limit']: raise ValueError(f"Contradictory equal constraints for {row['metric']}.")
    validate_constraint_intervals(result)
    return tuple(result)


def merge_constraints(*groups):
    result=[]
    for group in groups:
        for row in group:
            if row not in result: result.append(dict(row))
    return tuple(result)


def inherit_parent_plan(plan, parent, side):
    """Inherit this parent's side only; never union independent parents/pistons."""
    from dataclasses import replace
    from .artifacts import MechanismArtifact
    artifact=MechanismArtifact.from_data(parent['mechanisms'][side])
    merge = strongest_constraints if plan.request.mechanism_family == 'four_bar' else merge_constraints
    constraints=merge(artifact.scientific['constraints'],plan.request.mechanical_constraints)
    return replace(plan, request=replace(plan.request, mechanical_constraints=constraints))


def screen_provenance(plan, parent=None, side=None, *, requested_plan=None):
    screen=mechanical_screen(plan.request.mechanical_screen)
    inherited=[]
    if parent is not None and side in parent['mechanisms']:
        old=parent['mechanisms'][side]['provenance'].get('mechanical_screen')
        if old: inherited.append(old)
    active=screen if screen['name']!='none' or not inherited else inherited[0]
    result = dict(name=active['name'], version=active['version'], requested=screen['name'],
                activation='selected' if screen['name']!='none' else 'inherited' if inherited else 'none', inherited=inherited,
                effective_constraints=list(plan.request.mechanical_constraints),
                final_mechanical_samples=FINAL_MECHANICAL_SAMPLES if plan.request.mechanism_family in ('four_bar','six_bar') and plan.request.stages!=('primary_discovery',) else None)
    if plan.request.mechanism_family == 'four_bar':
        result['envelope_frame'] = 'slider_axis' if active['name']=='four_bar_design' else 'artifact_setting_or_global'
        result['parent_requirements'] = [] if parent is None else parent['mechanisms'][side]['scientific']['constraints']
        result['inheritance_rule'] = 'Strongest requirement per metric/relation; parent requirements cannot be relaxed.'
        requested = (requested_plan or plan).request.mechanical_constraints
        result['parent_prevents_relaxation'] = [
            dict(metric=old['metric'],relation=old['relation'],requested_limit=new['limit'],parent_limit=old['limit'])
            for old in result['parent_requirements'] for new in requested
            if (old['metric'],old['relation']) == (new['metric'],new['relation']) and
            ((old['relation']=='minimum' and old['limit']>new['limit']) or
             (old['relation']=='maximum' and old['limit']<new['limit']))]
    return result


@dataclass(frozen=True)
class PrimaryMechanicalPolicy:
    """Primary-only search filters; False/None explicitly disables one limit."""
    minimum_primary_transmission_sine: float | None = .30
    minimum_e_span_over_crank: float | None = .75

    def __post_init__(self):
        for name,value in asdict(self).items():
            if value is False:
                object.__setattr__(self,name,None)
            elif value is not None and (isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0):
                raise ValueError(f'primary_mechanical.{name} must be finite, nonnegative, or false.')
        value=self.minimum_primary_transmission_sine
        if value is not None and value>1:
            raise ValueError('Primary minimum transmission sine cannot exceed one.')

    def profile(self):
        resource=tomllib.loads(files('dada_solver.research').joinpath('data/six_bar_primary_design.toml').read_text())
        return dict(name=resource['name'],version=resource['version'],thresholds=asdict(self),
                    final_mechanical_samples=resource['canonical_samples'],scope='primary_discovery_only')

    def constraints(self, metrics):
        from .margins import margin_record
        rows=[]
        for metric,field,unit in (
                ('minimum_primary_transmission_sine','minimum_primary_transmission_sine','1'),
                ('E_span_over_crank','minimum_e_span_over_crank','crank_radius')):
            limit=getattr(self,field)
            if limit is None: continue
            record=margin_record(metric,metrics[metric],limit,'minimum',unit,
                                 method='primary full-cycle design screen; 32-epsilon roundoff tolerance')
            # Roundoff only: never move a scientific value or accept a meaningful violation.
            if record['available'] and record['margin']<0 and abs(record['margin'])<=32*2.220446049250313e-16*max(1.,limit):
                record.update(satisfied=True,state='satisfied',near_active=True)
            rows.append(record)
        return rows

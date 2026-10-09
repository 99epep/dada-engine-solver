"""Packaged design preferences promoted to explicit synthesis constraints."""
from importlib.resources import files
from dataclasses import dataclass, asdict
import math
import tomllib

FINAL_MECHANICAL_SAMPLES = 1440


def mechanical_screen(name):
    if name in (None, 'none'):
        return dict(name='none', version=None, constraints=[])
    if name != 'six_bar_design':
        raise ValueError(f'Unknown mechanical screen {name!r}.')
    return tomllib.loads(files('dada_solver.research').joinpath('data/six_bar_design.toml').read_text())


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
    constraints=merge_constraints(artifact.scientific['constraints'],plan.request.mechanical_constraints)
    return replace(plan, request=replace(plan.request, mechanical_constraints=constraints))


def screen_provenance(plan, parent=None, side=None):
    screen=mechanical_screen(plan.request.mechanical_screen)
    inherited=[]
    if parent is not None and side in parent['mechanisms']:
        old=parent['mechanisms'][side]['provenance'].get('mechanical_screen')
        if old: inherited.append(old)
    active=screen if screen['name']!='none' or not inherited else inherited[0]
    return dict(name=active['name'], version=active['version'], requested=screen['name'],
                activation='selected' if screen['name']!='none' else 'inherited' if inherited else 'none', inherited=inherited,
                effective_constraints=list(plan.request.mechanical_constraints),
                final_mechanical_samples=FINAL_MECHANICAL_SAMPLES if plan.request.mechanism_family=='six_bar' and plan.request.stages!=('primary_discovery',) else None)


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

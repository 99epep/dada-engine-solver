"""Portable local-region study creation; no physics change or evaluation import."""
import copy
import math
from pathlib import Path

from .report import inspect, select_records
from .snapshots import stored_study
from .study_export import write_study


def refine(sources, selectors, radius, output, *, tidy=False, standalone=False, name=None, purpose=None):
    if isinstance(radius,bool) or not math.isfinite(radius) or not 0<radius<=1:
        raise ValueError('Refinement radius must be finite and in (0, 1].')
    output=Path(output)
    if output.suffix!='.toml': raise ValueError('Refinement output must end in .toml.')
    datasets=[inspect(source) for source in sources]
    if not datasets: raise ValueError('Provide at least one source.')
    if len({d['study_id'] for d in datasets})!=1:
        raise ValueError('Refinement centers have incompatible scientific study identities.')
    if selectors:
        chosen=select_records({'records':[r for d in datasets for r in d['records']]},selectors)
    else:
        if any(Path(d['source']).is_dir() for d in datasets):
            raise ValueError('A campaign source requires explicit --candidate selectors.')
        chosen=[d['records'][0] for d in datasets]
    chosen=list({r['candidate_id']:r for r in chosen}.values())
    with stored_study(datasets[0]) as study:
        if study.data['schema_version'] != 3: raise ValueError('Refinement requires a schema-3 study.')
        raw=copy.deepcopy(study.data)
        if not study.space.parameters: raise ValueError('Refinement requires active ordered parameters.')
        for record in chosen: study.space.encode(record['physical'])
        # Initials aid standalone evaluate; bounds/transforms and fixed scientific
        # inputs stay untouched. Region centers carry exact source float values.
        for row in raw['parameters']:
            if 'initial' in row: row['initial']=chosen[0]['physical'][row['name']]
        raw['search']=dict(type='sobol',domain='local_regions_v1',seed=raw['search']['seed'],
            scramble=raw['search']['scramble'],radius_fraction=radius,allocation='round_robin',
            evaluate_centers=True,regions=[dict(id=f'basin_{i+1}',source_candidate_id=r['candidate_id'],
                source_study_id=datasets[0]['study_id'],center=r['physical']) for i,r in enumerate(chosen)])
        raw['study']['name'] += ' — local refinement'
        raw['execution']['default_max_candidates']=512
        return write_study(study,raw,output,tidy=tidy,standalone=standalone,name=name,purpose=purpose)

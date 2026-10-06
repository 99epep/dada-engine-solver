"""Portable local-region study creation; no physics change or evaluation import."""
import copy
import math
from pathlib import Path
import tempfile

from .report import inspect, select_records
from .snapshots import stored_study
from .schema import load_study
from .study_io import dumps


def refine(sources, selectors, radius, output):
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
        basis_path=output.with_suffix('.basis.json')
        raw['sources']['machine']['path']=basis_path.name
        contents={basis_path:study.basis.source}
        for side,artifact in study.artifacts.items():
            import json
            target=output.with_name(output.stem+'.'+side+'.mechanism.json')
            raw['kinematics'][side]['artifact']=target.name
            contents[target]=json.dumps(artifact.data,indent=2)+'\n'
        contents[output]=dumps(raw)
        if any(path.exists() for path in contents): raise ValueError('Refinement study or associated input already exists.')
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for path,text in contents.items(): (root/path.name).write_text(text)
            load_study(root/output.name)
        output.parent.mkdir(parents=True,exist_ok=True)
        for path,text in contents.items():
            with path.open('x') as stream: stream.write(text)
    return output

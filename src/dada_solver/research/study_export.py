"""Portable study formatting and optional removal of historical metadata."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile

from .schema import load_study, compile_study
from .study_io import dumps


def write_study(study, raw, output, *, tidy=False, standalone=False, name=None, purpose=None):
    """Preserve scientific inputs and search windows; never evaluate a cycle."""
    output=Path(output)
    if output.suffix!='.toml': raise ValueError('Study output must end in .toml.')
    raw=copy.deepcopy(raw)
    if name is not None: raw['study']['name']=name
    if purpose is not None: raw['study']['purpose']=purpose
    basis=study.basis.source
    if standalone:
        raw['study'].pop('parent_candidate_id',None)
        for index,region in enumerate(raw['search'].get('regions',[]),1):
            region.pop('source_candidate_id',None)
            region.pop('source_study_id',None)
            region['id']=f'region_{index}'
        data=copy.deepcopy(study.basis.data)
        data['provenance']={}
        if data['warm_start'] is not None:
            data['warm_start'].pop('source_candidate_id',None)
        basis=json.dumps(data,indent=2,allow_nan=False)+'\n'
    basis_path=output.with_suffix('.basis.json')
    raw['sources']['machine']=dict(path=basis_path.name,sha256=hashlib.sha256(basis.encode()).hexdigest())
    contents={basis_path:basis}
    for side,artifact in study.artifacts.items():
        target=output.with_name(output.stem+'.'+side+'.mechanism.json')
        data=artifact.data
        if standalone:
            data=copy.deepcopy(data);data['provenance']={}
        raw['kinematics'][side]['artifact']=target.name
        contents[target]=json.dumps(data,indent=2,allow_nan=False)+'\n'
    contents[output]=dumps(raw,tidy=tidy)
    if any(path.exists() for path in contents):
        raise ValueError('Study or associated input already exists; choose a new output.')
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory)
        for path,text in contents.items(): (root/path.name).write_text(text)
        generated=load_study(root/output.name)
        definition=compile_study(generated)
        definition.adapter.build(dict(definition.fixed_parameters,
            **{p.name:p.initial for p in generated.space.parameters}))
    output.parent.mkdir(parents=True,exist_ok=True)
    created=[]
    try:
        for path,text in contents.items():
            with path.open('x') as stream:
                created.append(path);stream.write(text)
    except BaseException:
        for path in created: path.unlink()
        raise
    return output


def export_study(source, output, *, tidy=False, standalone=False, name=None, purpose=None):
    study=load_study(source)
    return write_study(study,study.data,output,tidy=tidy,standalone=standalone,name=name,purpose=purpose)

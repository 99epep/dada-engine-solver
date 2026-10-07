"""Explicit human selection of independent physical piston mechanisms."""
from pathlib import Path

from dada_solver.campaign.candidate import content_hash
from .artifacts import MechanismArtifact, MechanismLibrary
from .families import PHYSICAL_FAMILIES


def pair_mechanisms(*, small_library, small, large_library, large, output):
    """Write one reusable pair without changing geometry or evaluating a machine."""
    mechanisms={};sources={}
    libraries={}
    for side,path,family_id in (('small',small_library,small),('large',large_library,large)):
        path=Path(path)
        key=path.resolve()
        if key not in libraries:
            try:
                libraries[key]=MechanismLibrary.load(path)
            except (ValueError,KeyError,TypeError,AttributeError) as error:
                raise ValueError(f'Invalid {side.upper()} mechanism library {path}: {error}') from error
        library=libraries[key];member=library.member(family_id)
        if side not in member['mechanisms']:
            raise ValueError(f'{side.upper()} selection {family_id!r} contains no {side.upper()} artifact.')
        raw=member['mechanisms'][side];artifact=MechanismArtifact.from_data(raw)
        settings=artifact.scientific['settings'];family=settings['family']
        if family not in PHYSICAL_FAMILIES or settings.get('component') is not None:
            raise ValueError(f'{side.upper()} requires a complete physical piston mechanism, not a primary/intermediate component.')
        # Retain the selected payload, including its provenance, exactly.
        mechanisms[side]=raw
        sources[side]=dict(source_library=str(path),source_library_hash=library.to_data()['content_hash'],
            source_family_id=family_id,artifact_hash=artifact.content_hash,physical_family=family)
    identity={s:{k:v for k,v in data.items() if k in ('source_family_id','artifact_hash')} for s,data in sources.items()}
    family_id='pair/'+content_hash(identity)[:16]
    result=MechanismLibrary((dict(family_id=family_id,mechanisms=mechanisms,
        metadata=dict(families={s:v['physical_family'] for s,v in sources.items()},
            provenance=dict(operation='manual_pair_selection',**sources))),))
    result.save(output)
    return result

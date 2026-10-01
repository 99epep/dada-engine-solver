"""Prepare the requested 16D-to-18D valve-topology study without evaluation."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile

from dada_solver.research.report import inspect, select_records
from dada_solver.research.schema import load_study
from dada_solver.research.study_io import dumps

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'outputs/pedal_cell_retro/5L_2p8Hz_coupled/campaign'
GLOBAL = ROOT/'outputs/pedal_cell_retro/5L_2p8Hz_coupled_valves18d/study.toml'
OUTPUT = ROOT/'outputs/pedal_cell_retro/5L_2p8Hz_coupled_valves18d_local/study.toml'
CANDIDATES = (
    'e02fcc89d7976753ba37415ab542c7974b47120504a4ada85c3e6655ad9c627a',
    '1e31322584b8359a39f13c2e24e9dcc1b3e445653aa18424cd25ac41284a989c',
    '22107e0bbbb88b9104240c12fdd5a95d62753d262d4a6d6b5c5a1f20c60826ea',
)
PLACEMENTS = ('valve.heat_in.placement', 'valve.heat_out.placement')
TOPOLOGIES = dict(DD=('downstream','downstream'), UD=('upstream','downstream'),
                  DU=('downstream','upstream'), UU=('upstream','upstream'))


def prepare(output=OUTPUT):
    """Explicitly lift selected coordinates into the supplied global 18D study."""
    output=Path(output)
    data=inspect(SOURCE)
    selected=select_records(data,list(CANDIDATES))
    study=load_study(GLOBAL)
    raw=copy.deepcopy(study.data)
    names={p.name for p in study.space.parameters}
    if len(names)!=18 or not set(PLACEMENTS)<=names:
        raise ValueError('Expected the global 18D valve-placement parameter space.')
    regions=[]
    for record in selected:
        if set(record['physical'])!=names-set(PLACEMENTS):
            raise ValueError('Source candidate must supply exactly the 16 numerical coordinates.')
        for topology,placements in TOPOLOGIES.items():
            center=dict(record['physical'],**dict(zip(PLACEMENTS,placements)))
            study.space.encode(center)
            regions.append(dict(id=f"basin_{record['candidate_id'][:12]}_{topology}",
                source_candidate_id=record['candidate_id'],source_study_id=data['study_id'],center=center))
    raw['search']=dict(type='sobol',domain='local_regions_v1',seed=raw['search']['seed'],
        scramble=raw['search']['scramble'],radius_fraction=.20,allocation='round_robin',
        evaluate_centers=True,regions=regions)
    for row in raw['parameters']:
        if 'initial' in row: row['initial']=regions[0]['center'][row['name']]
    raw['study']['name']='Pedal Cell Retro — 5 L / 2.8 Hz — 18D twelve-region valve refinement'
    raw['execution']['default_max_candidates']=4096
    raw['sources']['machine']['path']='study.basis.json'
    files={'study.basis.json':study.basis.source}
    for side,artifact in study.artifacts.items():
        name=f'study.{side}.mechanism.json'
        raw['kinematics'][side]['artifact']=name
        files[name]=json.dumps(artifact.data,indent=2)+'\n'
    files[output.name]=dumps(raw)
    provenance=dict(source_campaign=str(SOURCE.relative_to(ROOT)),source_study_id=data['study_id'],
        target_global_study=str(GLOBAL.relative_to(ROOT)),target_global_study_id=study.study_id,
        target_global_toml_sha256=hashlib.sha256(GLOBAL.read_bytes()).hexdigest(),
        source_candidates=list(CANDIDATES),topology_order=list(TOPOLOGIES),
        transformation='Exact source numerical coordinates plus explicit valve choices; no evaluation.',
        radius_fraction=.20,choice_behavior='Binary bin centers at 0.25/0.75 retain their choice at radius 0.20.')
    files['provenance.json']=json.dumps(provenance,indent=2)+'\n'
    if any((output.parent/name).exists() for name in files):
        raise ValueError('Destination study or associated inputs already exist.')
    with tempfile.TemporaryDirectory() as directory:
        for name,text in files.items(): (Path(directory)/name).write_text(text)
        load_study(Path(directory)/output.name)
    output.parent.mkdir(parents=True,exist_ok=True)
    for name,text in files.items():
        with (output.parent/name).open('x') as stream: stream.write(text)
    return output


if __name__=='__main__':
    print(prepare())

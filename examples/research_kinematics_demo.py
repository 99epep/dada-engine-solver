"""Run bounded Research V2 acceptance cases, not a cross-family optimization.

PYTHONPATH=src python3 examples/research_kinematics_demo.py --output outputs/research_kinematics_v2
Existing destinations are refused. Ten exact evaluations at most are requested:
eight fixed studies and two Sobol points separated by a real disk-based resume.
"""
import argparse
from pathlib import Path
import json
import tomllib
from dada_solver.campaign.history import atomic_json
from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.research.cli import evaluate
from dada_solver.research.presets import initialize_v2
from dada_solver.research.schema import load_study,compile_study
from dada_solver.research.study_io import dumps
from dada_solver.research.report import compare,render_html,inspect
from dada_solver.research.visualization import sample_motion


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.output
    if root.exists(): raise ValueError('Choose a new demonstration directory; existing results are preserved.')
    root.mkdir(parents=True)
    paths=[];summary=[]
    for family in ('harmonic','slider_crank','four_bar','six_bar','free_spline','fourier_c2','ideal_piecewise','structured_c2_15p'):
        path=initialize_v2(root/family/'study.toml',family,family,champion=family=='structured_c2_15p')
        output=path.parent/'evaluation.json'
        result=evaluate(path,output,budget='180s')
        paths.append(output)
        summary.append(dict(family=family,status=result['status'],candidate_id=result['candidate_id'],
            indicated_power_w=result['metrics'].get('indicated_power_w'),efficiency=result['metrics'].get('indicated_thermal_efficiency'),
            cycles=result['periodic_cycle_count'],reason=result['reason']))
        print(json.dumps(summary[-1]),flush=True)
        render_html(inspect(output),path.parent/'report.html')
        if family in ('six_bar','structured_c2_15p'):
            atomic_json(path.parent/'motion.json',sample_motion(load_study(path),samples=73))
    path=initialize_v2(root/'active_phase'/'study.toml')
    raw=tomllib.loads(path.read_text())
    row=next(r for r in raw['parameters'] if r['name']=='kinematics.small.phase_rad')
    initial=row.pop('value');row.update(kind='continuous',initial=initial,lower=initial-.02,upper=initial+.02,transform='linear')
    path.write_text(dumps(raw))
    campaign=OptimizationCampaign(compile_study(load_study(path)),path.parent/'run')
    first=campaign.run(180,maximum_candidates=1)
    resumed=OptimizationCampaign.resume(path.parent/'run')
    second=resumed.run(180,maximum_candidates=1)
    records=resumed.history.load()
    assert [r['sequence_index'] for r in records]==[0,1]
    assert len({r['candidate_id'] for r in records})==2
    assert all(r['integrated'] and r['converged'] for r in records)
    render_html(inspect(path.parent/'run'),path.parent/'report.html')
    render_html(compare(paths+ [path.parent/'run']),root/'comparison.html')
    atomic_json(root/'summary.json',dict(purpose='Bounded integration and usability demonstration',
        fair_optimization_comparison=False,notes=[
            'The first seven cases share a thermal basis but have independently chosen reference motions, not optimized families.',
            'The structured case replays candidate 3952 with its own original hardware and warm state.',
            'Indicated gas power is not useful shaft power. No mechanical efficiency is assumed.',
            'No mechanism synthesis or large optimization is executed.'],
        evaluations=summary,sobol=dict(phases=[first['attempted'],second['attempted']],
            candidate_ids=[r['candidate_id'] for r in records],sequence_indices=[r['sequence_index'] for r in records])))
    print('Saved bounded demonstration and offline comparison to '+str(root),flush=True)

if __name__=='__main__':main()

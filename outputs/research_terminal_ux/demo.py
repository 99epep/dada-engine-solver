"""Four synthetic evaluations through the real runner; no thermodynamic integration."""
from pathlib import Path
import io
import json
import sys
import tomllib
from dada_solver.campaign.candidate import content_hash
from dada_solver.campaign.evaluator import rejected
from dada_solver.campaign.runner import OptimizationCampaign
from dada_solver.research.presets import initialize_v2
from dada_solver.research.progress import CLIProgress
from dada_solver.research.schema import load_study, compile_study
from dada_solver.research.study_io import dumps

ROOT = Path(__file__).resolve().parent


class Clock:
    value = 0.
    def __call__(self): return self.value


class SyntheticEvaluator:
    def __init__(self,clock): self.clock=clock;self.index=0
    def evaluate_with_control(self,candidate,control):
        self.clock.value+=11
        control.check()
        self.clock.value+=2
        i=self.index;self.index+=1
        status=('feasible','converged_infeasible','invalid_exchanger','feasible')[i]
        result=rejected(status,'')
        result.update(integrated=i!=2,converged=i!=2)
        if i in (0,3):
            cop,q,p,m=(1.1472,151.8,132.3,.061) if i==0 else (1.1892,286.4,240.8,.078)
            result.update(objective=dict(name='synthetic_negative_cop',value=-cop,available=True),
                metrics=dict(cooling_cop=cop,cooling_power_w=q,indicated_mechanical_input_power_w=p,
                             maximum_absolute_mass_flow_kg_s=m))
        elif i==1:
            result.update(reason='maximum_absolute_mass_flow: violated',
                metrics=dict(maximum_absolute_mass_flow_kg_s=.084),
                constraints=[dict(name='maximum_absolute_mass_flow',satisfied=False,available=True,
                                  value=.084,limit=.08,relation='maximum',margin=-.004)])
        else: result['reason']='MicrotubeDomainError: large_relative_pressure_drop'
        return result


class Capture:
    def __init__(self,stream): self.stream=stream;self.raw=io.StringIO()
    def write(self,text): self.raw.write(text);return self.stream.write(text)
    def flush(self): self.stream.flush()
    def isatty(self): return self.stream.isatty()
    def fileno(self): return self.stream.fileno()


def main():
    path=initialize_v2(ROOT/'study.toml','harmonic','harmonic')
    raw=tomllib.loads(path.read_text())
    row=next(p for p in raw['parameters'] if p['name']=='volume.swept_ratio')
    value=row.pop('value');row.update(initial=value,lower=value*.8,upper=value*1.2,kind='continuous',transform='linear')
    raw['study']['name']='Synthetic terminal UX demonstration — no physical results'
    raw['search']=dict(type='sobol',domain='local_regions_v1',seed=42,scramble=True,
        radius_fraction=.1,allocation='round_robin',evaluate_centers=True,
        regions=[dict(id=f'basin_{i+1}',source_candidate_id=content_hash(dict(synthetic_center=i)),
                      source_study_id=content_hash(dict(synthetic_source=True)),
                      center={'volume.swept_ratio':value*(1+.05*i)}) for i in range(2)])
    path.write_text(dumps(raw))
    definition=compile_study(load_study(path));clock=Clock();capture=Capture(sys.stdout)
    events=[]
    with (ROOT/'run.log').open('w') as logfile, CLIProgress(capture,scientific=definition.study.scientific) as tty, CLIProgress(logfile,scientific=definition.study.scientific) as plain:
        def notify(event):
            events.append(event);tty(event);plain(event)
        campaign=OptimizationCampaign(definition,ROOT/'campaign',evaluator=SyntheticEvaluator(clock),clock=clock)
        campaign.run(300,maximum_candidates=4,progress_callback=notify)
    (ROOT/'terminal_capture.txt').write_text(capture.raw.getvalue())
    (ROOT/'events.json').write_text(json.dumps(events,indent=2)+'\n')


if __name__=='__main__': main()

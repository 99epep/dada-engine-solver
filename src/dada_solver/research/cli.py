"""Configure, evaluate, search and compare traceable research studies."""
import argparse
from datetime import datetime, timezone
from importlib.resources import files
import json
from pathlib import Path
import time

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.evaluator import MachineEvaluator, EvaluationControl
from dada_solver.campaign.history import atomic_json
from dada_solver.campaign.runner import OptimizationCampaign, parse_budget
from .schema import load_study, compile_study, candidate_for_values
from .sixbar import PARAMETERS
from . import report


def initialize(output):
    output = Path(output)
    if output.suffix != '.toml': raise ValueError('Study output must end in .toml.')
    basis = output.with_suffix('.basis.json')
    if output.exists() or basis.exists(): raise ValueError('Study or basis already exists; choose a new output name.')
    resources = files('dada_solver.research').joinpath('data')
    source = resources.joinpath('sixbar-thermo5d.toml').read_text().replace('path = "rank01_basis.json"', f'path = {json.dumps(basis.name, ensure_ascii=False)}')
    output.parent.mkdir(parents=True, exist_ok=True)
    with basis.open('x') as stream: stream.write(resources.joinpath('rank01_basis.json').read_text())
    with output.open('x') as stream: stream.write(source)
    return output


def evaluate(study_path, output, *, assignments=(), reference=False, budget=None):
    output = Path(output)
    if output.exists(): raise ValueError('Evaluation output already exists; choose a new filename.')
    study = load_study(study_path)
    definition = compile_study(study)
    values = {p.name:p.initial for p in study.space.parameters}
    if reference:
        if study.data['schema_version']!=1: raise ValueError('V2 uses configured fixed/initial values; --reference is a V1 regression selector.')
        values = {name:study.basis.data['reference_parameters'][spec[2]] for name,spec in PARAMETERS.items()}
    aliases = {spec[2]:name for name,spec in PARAMETERS.items()}
    owned={p.name:p for p in study.space.parameters}
    assigned = set()
    for assignment in assignments:
        name, sep, value = assignment.partition('=')
        name = aliases.get(name,name)
        if not sep or name not in owned or name in assigned:
            raise ValueError(f'Expected a unique known parameter=value, got {assignment!r}.')
        assigned.add(name)
        from dada_solver.campaign.parameters import IntegerParameter, ChoiceParameter
        values[name] = value if isinstance(owned[name],ChoiceParameter) else int(value) if isinstance(owned[name],IntegerParameter) else float(value)
    candidate = candidate_for_values(definition, values)
    started = time.monotonic()
    seconds = parse_budget(budget) if budget is not None else definition.candidate_budget_seconds
    if seconds <= 0: raise ValueError('Evaluation budget must be positive.')
    result = MachineEvaluator(definition).evaluate_with_control(candidate, EvaluationControl(deadline=started+seconds))
    record = dict(**candidate.payload, **result, candidate_id=candidate.candidate_id,
        duration_seconds=time.monotonic()-started, timestamp=datetime.now(timezone.utc).isoformat(),
        cache_hit=False, evaluation_number=0, sequence_index=None, phase_id=0)
    artifact = dict(artifact_type='research_evaluation_v1', name=study.data['study']['name'],
        definition=dict(definition.identity, definition_id=definition.definition_id), record=record,
        selection='reference_with_explicit_overrides' if reference else 'configured_initial_values_with_explicit_overrides')
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(output, artifact)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    init = commands.add_parser('init', help='Create an editable study and its portable basis')
    init.add_argument('preset', choices=['sixbar-thermo5d','kinematics','structured-c2-3952','external-stream-refrigeration','external-stream-motor'])
    from .families import FAMILIES
    init.add_argument('--small',choices=FAMILIES,default='harmonic')
    init.add_argument('--large',choices=FAMILIES,default='harmonic')
    init.add_argument('--coupling',choices=['independent','shared_crank'],default='independent')
    init.add_argument('--output', type=Path, required=True)
    validate = commands.add_parser('validate', help='Validate configuration and geometry without integration')
    validate.add_argument('study', type=Path)
    validate.add_argument('--json', action='store_true', help='Print full validation and runtime metadata')
    run = commands.add_parser('run', help='Start a bounded Sobol campaign')
    run.add_argument('study', type=Path)
    run.add_argument('--directory', type=Path, required=True)
    resume = commands.add_parser('resume', help='Continue a compatible stored campaign')
    resume.add_argument('directory', type=Path)
    resume.add_argument('--retry-incomplete', action='store_true')
    for p in (run,resume):
        p.add_argument('--budget', help='Wall-clock budget such as 2m or 1h30m')
        p.add_argument('--max-candidates', type=int, help='Maximum attempts in this invocation')
    single = commands.add_parser('evaluate', help='Evaluate one exact configuration without a search')
    single.add_argument('study', type=Path)
    single.add_argument('--output', required=True, type=Path)
    single.add_argument('--reference', action='store_true', help='Start from the stored regression candidate instead of TOML initials')
    single.add_argument('--set', action='append', default=[], metavar='PARAMETER=VALUE')
    single.add_argument('--budget')
    refinement = commands.add_parser('refine', help='Create a portable center-first local Sobol study')
    refinement.add_argument('sources',nargs='+',type=Path)
    refinement.add_argument('--candidate',action='append',default=[],help='ID, unique prefix or best; repeat for several centers')
    refinement.add_argument('--radius',type=float,required=True,help='Half-width in global normalized coordinates, in (0,1]')
    refinement.add_argument('--output',type=Path,required=True)
    resize = commands.add_parser('rescale', help='Create a new capacity-scaled candidate study and portable basis; never integrate')
    resize.add_argument('source', type=Path)
    resize.add_argument('--candidate', required=True, help='Exact ID or unambiguous prefix')
    resize.add_argument('--factor', required=True, type=float)
    resize.add_argument('--mode', choices=['capacity'], default='capacity')
    resize.add_argument('--output', required=True, type=Path)
    for name in ('status','report','compare'):
        p = commands.add_parser(name, help='Read stored results; never run integration')
        p.add_argument('paths', nargs='+' if name=='compare' else 1, type=Path)
        p.add_argument('--candidate', action='append', default=[], help='Exact ID, unambiguous prefix or best; repeat to compare')
        p.add_argument('--json', action='store_true', help='Print the complete inspection dataset')
        p.add_argument('--list-candidates', action='store_true', help='List selected candidates in the terminal')
        if name != 'status':
            p.add_argument('--html', type=Path, help='Write or regenerate standalone HTML (report default: CAMPAIGN/report.html)')
            p.add_argument('--plots', choices=['volumes'], help='Sample selected candidate volumes without thermodynamic integration')
    args = parser.parse_args(argv)
    try:
        if args.command == 'init':
            if args.preset=='sixbar-thermo5d':
                if (args.small,args.large,args.coupling)!=('harmonic','harmonic','independent'):
                    raise ValueError('Use the kinematics preset to choose cylinder families.')
                path = initialize(args.output)
            elif args.preset.startswith('external-stream-'):
                from .presets import initialize_v3
                path=initialize_v3(args.output,args.small,args.large,mode=args.preset.removeprefix('external-stream-'))
            else:
                from .presets import initialize_v2
                path=initialize_v2(args.output,args.small,args.large,coupling=args.coupling,champion=args.preset=='structured-c2-3952')
            print(f'Created {path} and {path.with_suffix(".basis.json")}. Next: dada-research validate {path}')
        elif args.command == 'validate':
            study = load_study(args.study)
            definition = compile_study(study)
            validation = dict(valid=True, study_id=study.study_id, protocol=study.data['study']['protocol'],
                fixed=getattr(study,'fixed_parameters',study.data.get('fixed',{})), parameters=study.data['parameters'], backend=definition.numerical_settings['wall_backend'],
                integration_started=False)
            if args.json:
                print(json.dumps(validation, indent=2))
            else:
                print(f"Valid study: {study.data['study']['name']}\nStudy: {study.study_id}")
                if study.data['schema_version']==1:
                    print(f"Fixed six-bar pair; air inlets {study.data['fixed']['cold_air_inlet_K']} / {study.data['fixed']['hot_air_inlet_K']} K; no integration started.")
                else:
                    print(f"Families: {study.settings['small']['family']} / {study.settings['large']['family']}; {len(study.space.parameters)} active coordinates; no integration started.")
                for p in study.data['parameters']:
                    description=f"fixed {p['value']}" if 'value' in p else f"choices {p['choices']}; initial {p['initial']}" if p.get('kind')=='choice' else f"{p['lower']} .. {p['upper']}; initial {p['initial']}"
                    print(f"  {p['name']} [{p['unit']}]: {description}")
                if definition.wall_backend.name=='numba' and not validation['backend'].get('numba_available'):
                    print('Numba is unavailable; the existing Python fallback will be recorded in each result.')
        elif args.command in ('run','resume'):
            if args.budget is not None: parse_budget(args.budget)
            if args.max_candidates is not None and args.max_candidates < 0: raise ValueError('Maximum candidates must be nonnegative.')
            if args.command == 'run':
                if args.directory.exists() and any(args.directory.iterdir()): raise ValueError('Study directory is not empty; use resume or choose a new directory.')
                definition = compile_study(load_study(args.study))
                if not definition.space.parameters: raise ValueError('No active parameters. Use evaluate, or replace a fixed value by initial/lower/upper before starting a search.')
                campaign = OptimizationCampaign(definition,args.directory)
            else:
                definition = CampaignDefinition.resume(args.directory)
                if not hasattr(definition,'study'): raise ValueError('Use dada-optimize to resume a legacy campaign.')
                campaign = OptimizationCampaign(definition,args.directory)
            budget = parse_budget(args.budget or definition.study.data['execution']['default_budget'])
            from .progress import CLIProgress
            # Retire the historical 16-attempt Research default without editing snapshots.
            limit = 512 if args.max_candidates is None and definition.maximum_candidates == 16 else args.max_candidates
            with CLIProgress(scientific=definition.study.scientific) as progress:
                campaign.run(budget, maximum_candidates=limit, progress_callback=progress,
                             retry_incomplete=getattr(args,'retry_incomplete',False))
        elif args.command == 'refine':
            from .refine import refine
            path=refine(args.sources,args.candidate,args.radius,args.output)
            print(f'Created {path}. Centers will be evaluated before local Sobol sampling.')
        elif args.command == 'rescale':
            from .rescale import rescale
            path = rescale(args.source, args.candidate, args.factor, args.output, mode=args.mode)
            print(f'Created {path} and {path.with_suffix(".basis.json")}. Constraints remain unchanged; evaluate to verify scaling.')
        elif args.command == 'evaluate':
            record = evaluate(args.study, args.output, assignments=args.set, reference=args.reference, budget=args.budget)
            print(f"{record['candidate_id']} {record['status']}; saved {args.output}")
            print(json.dumps(record['metrics'],indent=2))
        else:
            data = report.compare(args.paths,args.candidate, plots=getattr(args,'plots',None))
            destination = getattr(args,'html',None)
            if args.command == 'report' and destination is None:
                source = args.paths[0]
                destination = source/'report.html' if source.is_dir() else source.with_suffix('.html')
            if destination is not None: report.render_html(data,destination)
            if args.json:
                print(json.dumps(data,indent=2))
            else:
                print(report.text_report(data, list_candidates=args.list_candidates),end='')
                if destination is not None: print(f'HTML: {destination}')
    except KeyboardInterrupt:
        parser.exit(130, 'Interrupted. Campaign pending state is preserved; use resume. Standalone evaluations must be requested again.\n')
    except (ValueError, OSError, KeyError, TypeError, RuntimeError) as error:
        parser.exit(2, f'Research error: {error}\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

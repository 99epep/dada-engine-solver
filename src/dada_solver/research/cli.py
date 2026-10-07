"""Configure, evaluate, search and compare traceable research studies."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
import sys

from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.evaluator import MachineEvaluator, EvaluationControl
from dada_solver.campaign.history import atomic_json
from dada_solver.campaign.runner import OptimizationCampaign, parse_budget
from .schema import load_study, compile_study, candidate_for_values
from . import report


def evaluate(study_path, output, *, assignments=(), budget=None):
    output = Path(output)
    if output.exists(): raise ValueError('Evaluation output already exists; choose a new filename.')
    study = load_study(study_path)
    definition = compile_study(study)
    values = {p.name:p.initial for p in study.space.parameters}
    owned={p.name:p for p in study.space.parameters}
    assigned = set()
    for assignment in assignments:
        name, sep, value = assignment.partition('=')
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
        selection='configured_initial_values_with_explicit_overrides')
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(output, artifact)
    return record


class _ConsistentSelection(argparse.Action):
    """Reject contradictory repeated selectors instead of silently taking the last."""
    def __call__(self, parser, namespace, values, option_string=None):
        previous=getattr(namespace,self.dest,None)
        if previous is not None and previous!=values:
            parser.error(f'{option_string} was repeated with conflicting values.')
        setattr(namespace,self.dest,values)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    init = commands.add_parser('init', help='Create an editable study and its portable basis')
    init.add_argument('preset', choices=['kinematics','external-stream-refrigeration','external-stream-motor'])
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
    run.add_argument('--directory', type=Path, help='Campaign directory (default: campaign beside the study TOML)')
    resume = commands.add_parser('resume', help='Continue a compatible stored campaign')
    resume.add_argument('directory', type=Path, help='Campaign directory, or study TOML to use its sibling campaign directory')
    resume.add_argument('--retry-incomplete', action='store_true')
    for p in (run,resume):
        p.add_argument('--budget', help='Wall-clock budget such as 2m or 1h30m')
        p.add_argument('--max-candidates', type=int, help='Maximum attempts in this invocation')
    single = commands.add_parser('evaluate', help='Evaluate one exact configuration without a search')
    single.add_argument('study', type=Path)
    single.add_argument('--output', required=True, type=Path)
    single.add_argument('--set', action='append', default=[], metavar='PARAMETER=VALUE')
    single.add_argument('--budget')
    refinement = commands.add_parser('refine', help='Create a portable center-first local Sobol study')
    refinement.add_argument('sources',nargs='+',type=Path)
    refinement.add_argument('--candidate',action='append',default=[],help='ID, unique prefix, best or second; repeat for several centers')
    refinement.add_argument('--radius',type=float,required=True,help='Half-width in global normalized coordinates, in (0,1]')
    refinement.add_argument('--output',type=Path,required=True)
    resize = commands.add_parser('rescale', help='Create a new capacity-scaled candidate study and portable basis; never integrate')
    resize.add_argument('source', type=Path)
    resize.add_argument('--candidate', required=True, help='Exact ID or unambiguous prefix')
    resize.add_argument('--factor', required=True, type=float)
    resize.add_argument('--mode', choices=['capacity'], default='capacity')
    resize.add_argument('--output', required=True, type=Path)
    target = commands.add_parser('motion-target', help='Extract periodic study-angle motion without integration')
    target.add_argument('source', type=Path)
    target.add_argument('--candidate')
    target.add_argument('--samples', type=int, default=721)
    target.add_argument('--output', type=Path, required=True)
    refit = commands.add_parser('motion-refit', help='Fit structured_c2_15p and create a motion-only Research study')
    refit.add_argument('source', type=Path)
    refit.add_argument('--candidate')
    refit.add_argument('--validate-only', action='store_true')
    refit.add_argument('--output', type=Path, help='New schema-3 study TOML; requires a complete Research source')
    refit.add_argument('--report', type=Path, help='Stable geometric refit report; sufficient for MotionTarget-only sources')
    refit.add_argument('--radius', type=float, default=.1, help='Local structured coordinate search radius in (0, 0.5]; default 0.1')
    refit.add_argument('--plot', action='store_true', help='Open static target/structured comparison with events and extrema')
    refit.add_argument('--no-show', action='store_true')
    mechanism = commands.add_parser('mechanism', help='Validate synthesis protocols and inspect physical artifacts')
    mechanical = mechanism.add_subparsers(dest='mechanism_command', required=True)
    pairing = mechanical.add_parser('pair', help='Select independent SMALL and LARGE mechanisms; no integration')
    pairing.add_argument('library', type=Path, nargs='?', help='Common library; cannot combine with side-library options')
    pairing.add_argument('--small-library', type=Path, action=_ConsistentSelection)
    pairing.add_argument('--large-library', type=Path, action=_ConsistentSelection)
    pairing.add_argument('--small', required=True, action=_ConsistentSelection, metavar='FAMILY_ID')
    pairing.add_argument('--large', required=True, action=_ConsistentSelection, metavar='FAMILY_ID')
    pairing.add_argument('--output', required=True, type=Path, action=_ConsistentSelection)
    adapt = mechanical.add_parser('adapt', help='Generate a local paired thermodynamic study; no integration')
    adapt.add_argument('source', type=Path, help='Research study, evaluation, or campaign')
    adapt.add_argument('--candidate', help='Required for a campaign source')
    adapt.add_argument('--library', type=Path, required=True)
    adapt.add_argument('--family-id', required=True, help='Selected complete SMALL/LARGE member')
    adapt.add_argument('--output', type=Path, required=True, help='New portable study TOML')
    adapt.add_argument('--radius', type=float, default=.1, help='Local normalized Sobol radius in (0, 0.5]; default 0.1')
    retune = mechanical.add_parser('retune', help='Generate a local hardware study with the adapted mechanical pair fixed; no integration')
    retune.add_argument('source', type=Path, help='Paired thermodynamic Research study, evaluation, campaign or retuning descendant')
    retune.add_argument('--candidate', help='Exact ID, prefix or best; required for campaign sources')
    retune.add_argument('--scope', choices=['source-active'], default=None, help='Reopen originally active non-kinematic parameters (default when no targeted selection)')
    from .hardware_retuning import GROUPS
    retune.add_argument('--group', choices=tuple(GROUPS), action='append', default=[])
    retune.add_argument('--parameter', action='append', default=[], metavar='NAME')
    retune.add_argument('--radius', type=float, default=.1, help='Radius in original normalized parameter domains, in (0, 1]; default 0.1')
    retune.add_argument('--source-study', type=Path, help='Original Research source if its domains were not recorded; scientific identity is checked')
    retune.add_argument('--output', type=Path, required=True)
    synthesize = mechanical.add_parser('synthesize', help='Discover or polish physical mechanism families without thermodynamic integration')
    synthesize.add_argument('source', type=Path)
    synthesize.add_argument('--candidate')
    from .families import PHYSICAL_FAMILIES
    synthesize.add_argument('--family', choices=PHYSICAL_FAMILIES, required=True)
    synthesize.add_argument('--stage', action='append', required=True)
    synthesize.add_argument('--side', choices=['small','large','both'], default='large')
    synthesize.add_argument('--validate-only', action='store_true')
    synthesize.add_argument('--output',type=Path,help='New MechanismLibrary JSON')
    synthesize.add_argument('--html',type=Path,help='Standalone human catalogue (default: beside library JSON)')
    synthesize.add_argument('--library',type=Path,help='Input library for downstream, polish or opposite-piston stages')
    synthesize.add_argument('--family-id',action='append',default=[],help='Retained parent member; repeat to preserve several families')
    synthesize.add_argument('--config',type=Path,help='TOML search policy, bounds, categories and mechanical constraints')
    synthesize.add_argument('--budget',help='Cooperative geometry-search budget, e.g. 2m; not a thermodynamic budget')
    for option in ('islands','population','generations','seed','max-evaluations'):
        synthesize.add_argument('--'+option,type=int,default=None)
    visual = mechanical.add_parser('visualize', help='Open interactive production geometry and motion, without integration')
    visual.add_argument('artifact', type=Path)
    visual.add_argument('--family-id', help='Required when inspecting a library')
    visual.add_argument('--side', choices=['small','large'], default='small')
    visual.add_argument('--target', type=Path, help='MotionTarget or current study TOML')
    visual.add_argument('--static', action='store_true')
    visual.add_argument('--no-show', action='store_true', help='Construct and validate the view without opening a window')
    catalogue = mechanical.add_parser('catalogue', help='List all retained families without ranking')
    catalogue.add_argument('library', type=Path)
    catalogue.add_argument('--plot', action='store_true')
    catalogue.add_argument('--no-show', action='store_true')
    for name in ('status','report','compare'):
        p = commands.add_parser(name, help='Inspect stored results; report curves can replay one saved-state cycle')
        p.add_argument('paths', nargs='+' if name=='compare' else 1, type=Path)
        p.add_argument('--candidate', action='append', default=[], help='Exact ID, unambiguous prefix, best or second; repeat to compare')
        p.add_argument('--json', action='store_true', help='Print the complete inspection dataset')
        p.add_argument('--list-candidates', action='store_true', help='List selected candidates in the terminal')
        if name != 'status':
            p.add_argument('--html', type=Path, help='Write or regenerate standalone HTML (report default: CAMPAIGN/report.html)')
            p.add_argument('--plots', action='append', metavar='NAMES', help='Comma-separated/repeated: positions, volumes, pressures, heat, temperatures, flows, velocity, acceleration, mechanisms; default, all or none. Report defaults to positions,heat,pressures,temperatures for the best two.')
    args = parser.parse_args(argv)
    try:
        if args.command == 'init':
            if args.preset.startswith('external-stream-'):
                from .presets import initialize_external_stream
                path=initialize_external_stream(args.output,args.small,args.large,mode=args.preset.removeprefix('external-stream-'))
            else:
                from .presets import initialize_kinematics
                path=initialize_kinematics(args.output,args.small,args.large,coupling=args.coupling)
            print(f'Created {path} and {path.with_suffix(".basis.json")}. Next: dada-research validate {path}')
        elif args.command == 'validate':
            study = load_study(args.study)
            definition = compile_study(study)
            validation = dict(valid=True, study_id=study.study_id, protocol=study.data['study']['protocol'],
                fixed=study.fixed_parameters, parameters=study.data['parameters'], backend=definition.numerical_settings['wall_backend'],
                integration_started=False)
            if args.json:
                print(json.dumps(validation, indent=2))
            else:
                print(f"Valid study: {study.data['study']['name']}\nStudy: {study.study_id}")
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
                args.directory = args.directory if args.directory is not None else args.study.parent/'campaign'
                if args.directory.exists() and any(args.directory.iterdir()): raise ValueError('Study directory is not empty; use resume or choose a new directory.')
                definition = compile_study(load_study(args.study))
                if not definition.space.parameters: raise ValueError('No active parameters. Use evaluate, or replace a fixed value by initial/lower/upper before starting a search.')
                campaign = OptimizationCampaign(definition,args.directory)
            else:
                if args.directory.suffix.lower() == '.toml' and not args.directory.is_dir():
                    args.directory = args.directory.parent/'campaign'
                if not args.directory.is_dir():
                    raise ValueError(f'Campaign directory not found: {args.directory}. Start it with run STUDY.toml first, or provide an existing campaign directory.')
                definition = CampaignDefinition.resume(args.directory)
                if not hasattr(definition,'study'): raise ValueError('Use dada-optimize to resume an optimization campaign.')
                campaign = OptimizationCampaign(definition,args.directory)
            budget = parse_budget(args.budget or definition.study.data['execution']['default_budget'])
            from .progress import CLIProgress
            limit = args.max_candidates
            with CLIProgress(scientific=definition.study.scientific) as progress:
                campaign.run(budget, maximum_candidates=limit, progress_callback=progress,
                             retry_incomplete=getattr(args,'retry_incomplete',False))
        elif args.command == 'refine':
            from .refine import refine
            path=refine(args.sources,args.candidate,args.radius,args.output)
            print(f'Created {path}\nCenters will be evaluated before local Sobol sampling.')
        elif args.command == 'rescale':
            from .rescale import rescale
            path = rescale(args.source, args.candidate, args.factor, args.output, mode=args.mode)
            print(f'Created {path} and {path.with_suffix(".basis.json")}. Constraints remain unchanged; evaluate to verify scaling.')
        elif args.command in ('motion-target', 'motion-refit'):
            from .motion_target import load_motion_target
            target = load_motion_target(args.source, candidate=args.candidate,
                                        samples=getattr(args, 'samples', 1441))
            if args.command == 'motion-target':
                target.save(args.output)
                print(f'Motion target {target.content_hash}; saved {args.output}; no integration started.')
            else:
                from .motion_refit import MotionRefitRequest
                request = MotionRefitRequest(target)
                if args.validate_only:
                    print(json.dumps(dict(target_hash=target.content_hash,family=request.destination_family,
                                          parameter_count=15,implemented=True,integration_started=False)))
                else:
                    if args.output is None and args.report is None:
                        raise ValueError('Motion refit requires --output STUDY.toml or --report REFIT.json.')
                    if args.report is not None and args.report.exists():
                        raise ValueError('Refit report already exists; choose a new path.')
                    if args.output is not None:
                        from .refit_study import refit_study
                        path,result = refit_study(args.source,args.output,candidate=args.candidate,report=args.report,
                            radius=args.radius)
                        print(f'Created {path} and {path.with_suffix(".basis.json")}; 15 active structured coordinates; no integration started.')
                    else:
                        result = request.execute()
                        result.save(args.report)
                    for side,fitted in result.data['scientific']['sides'].items():
                        d,n = fitted['diagnostics'],fitted['initial_diagnostics']
                        print(f"{side.upper()}: position RMS {n['position_rms']:.6g} -> {d['position_rms']:.6g}; max {d['maximum_absolute_position_error']:.6g}; {d['extrema_count']} extrema")
                    if args.report is not None: print(f'Refit report: {args.report}')
                    if args.plot:
                        from .motion_refit import plot_refit
                        import matplotlib.pyplot as plt
                        figure = plot_refit(result)
                        if args.no_show:
                            figure.canvas.draw(); plt.close(figure)
                        else: plt.show()
        elif args.command == 'mechanism':
            from .artifacts import MechanismArtifact, MechanismLibrary
            if args.mechanism_command == 'pair':
                if args.library is not None:
                    if args.small_library is not None or args.large_library is not None:
                        raise ValueError('Use a common positional library OR both --small-library and --large-library, never both forms.')
                    small_library=large_library=args.library
                else:
                    if args.small_library is None or args.large_library is None:
                        raise ValueError('Provide a common LIBRARY or both --small-library and --large-library.')
                    small_library,large_library=args.small_library,args.large_library
                from .mechanism_pairing import pair_mechanisms
                paired=pair_mechanisms(small_library=small_library,small=args.small,
                    large_library=large_library,large=args.large,output=args.output)
                member=paired.members[0]
                print(f"Created {args.output}; family_id {member['family_id']}; "
                      f"SMALL {member['metadata']['families']['small']}, LARGE {member['metadata']['families']['large']}; no integration started.")
            elif args.mechanism_command == 'adapt':
                from .mechanism_adaptation import paired_thermodynamic
                output = paired_thermodynamic(args.source,args.library,args.family_id,args.output,
                                               candidate=args.candidate,radius=args.radius)
                generated = load_study(output)
                print(f'Created {output} with {len(generated.space.parameters)} active mechanical coordinates; no integration started.')
            elif args.mechanism_command == 'retune':
                from .hardware_retuning import hardware_retuning
                output=hardware_retuning(args.source,args.output,candidate=args.candidate,scope=args.scope,
                    groups=args.group,parameters=args.parameter,radius=args.radius,source_study=args.source_study)
                generated=load_study(output)
                print(f'Created {output} with {len(generated.space.parameters)} active hardware parameters and 0 active mechanism coordinates; no integration started.')
            elif args.mechanism_command == 'synthesize':
                from .motion_target import load_motion_target
                from .synthesis import SynthesisRequest, SynthesisPlan, release_coordinates
                target=load_motion_target(args.source,candidate=args.candidate)
                from .synthesis_search import SearchPolicy
                import tomllib
                settings={} if args.config is None else tomllib.loads(args.config.read_text())
                constraints=tuple(settings.pop('mechanical_constraints',()))
                declared_limits=settings.pop('volume_limits',{})
                from dada_solver.geometry import CylinderVolumeLimits
                if set(declared_limits)-{'small','large'} or any(set(row)!={'minimum','maximum'} for row in declared_limits.values()):
                    raise ValueError('Volume limits require small/large sides and exactly minimum/maximum.')
                limits={side:CylinderVolumeLimits(row['minimum'],row['maximum']) for side,row in declared_limits.items()}
                for option in ('islands','population','generations','seed','max_evaluations'):
                    value=getattr(args,option)
                    if value is not None: settings[option]=value
                if args.budget is not None: settings['budget_seconds']=parse_budget(args.budget)
                policy=SearchPolicy(**settings)
                request=SynthesisRequest(target.content_hash,args.family,'design_exploitation',tuple(args.stage),
                                         retained_family_ids=tuple(args.family_id),mechanical_constraints=constraints)
                plan=SynthesisPlan(target,request)
                sides=('small','large') if args.side=='both' else (args.side,)
                released={stage:release_coordinates(stage,sides,family=args.family) for stage in args.stage}
                supported=(args.family in ('slider_crank','four_bar') and tuple(args.stage) in (('global_discovery',),('full_local_polish',))) or (args.family=='six_bar' and len(args.stage)==1 and args.stage[0] in ('primary_discovery','downstream_fit','full_local_polish','mirror_initialization','opposite_local_adaptation'))
                if args.validate_only:
                    policy.coordinates(args.family); policy.choices(args.family)
                    print(json.dumps(dict(target_hash=target.content_hash,family=args.family,
                        stages=args.stage,released_coordinates=released,implemented=supported,integration_started=False)))
                else:
                    if not supported: plan.execute(policy=policy,sides=sides)
                    if args.output is None: raise ValueError('Synthesis requires --output for its MechanismLibrary.')
                    catalogue=args.html if args.html is not None else args.output.with_suffix('.html')
                    if args.output.exists() or catalogue.exists() or args.output.resolve()==catalogue.resolve():
                        raise ValueError('Synthesis library/catalogue paths must be distinct and must not exist.')
                    library=None if args.library is None else MechanismLibrary.load(args.library)
                    import matplotlib  # Verify the catalogue dependency before starting search.
                    result=plan.execute(policy=policy,sides=sides,library=library,volume_limits=limits)
                    from .synthesis_catalogue import render_synthesis_catalogue
                    render_synthesis_catalogue(result,target,catalogue,library_path=args.output)
                    result.save(args.output)
                    print(f'Created {args.output} and {catalogue}; {len(result.members)} retained families; no thermodynamic integration started.')
            elif args.mechanism_command == 'catalogue':
                from .synthesis import mechanism_catalogue
                library = MechanismLibrary.load(args.library)
                print(json.dumps(mechanism_catalogue(library), indent=2))
                if args.plot:
                    from .mechanism_view import plot_library_catalogue
                    import matplotlib.pyplot as plt
                    figure = plot_library_catalogue(library)
                    plt.close(figure) if args.no_show else plt.show()
            else:
                raw = json.loads(args.artifact.read_text())
                thermodynamic = None
                if raw.get('artifact_type') == 'mechanism_library':
                    if args.family_id is None:
                        raise ValueError('Choose an explicit --family-id from the library catalogue.')
                    member = MechanismLibrary.from_data(raw).member(args.family_id)
                    if args.side not in member['mechanisms']:
                        raise ValueError('Selected family has no artifact for that piston.')
                    artifact = MechanismArtifact.from_data(member['mechanisms'][args.side])
                    thermodynamic = member['metadata'].get('thermodynamic')
                else:
                    if args.family_id is not None:
                        raise ValueError('--family-id applies only to a mechanism library.')
                    artifact = MechanismArtifact.from_data(raw)
                from .motion_target import load_motion_target
                target = load_motion_target(args.target) if args.target else None
                from .mechanism_view import plot_mechanism, animate_mechanism
                import matplotlib.pyplot as plt
                if args.static:
                    figure = plot_mechanism(artifact, target=target, side=args.side, thermodynamic=thermodynamic)
                else:
                    figure, animation = animate_mechanism(artifact, target=target, side=args.side, thermodynamic=thermodynamic)
                print(f"Viewing {artifact.content_hash}; no thermodynamic integration started.")
                if args.no_show:
                    figure.canvas.draw()
                    plt.close(figure)
                else:
                    plt.show()
        elif args.command == 'evaluate':
            record = evaluate(args.study, args.output, assignments=args.set, budget=args.budget)
            print(f"{record['candidate_id']} {record['status']}; saved {args.output}")
            print(json.dumps(record['metrics'],indent=2))
        else:
            destination = getattr(args,'html',None)
            if args.command == 'report' and destination is None:
                source = args.paths[0]
                destination = source/'report.html' if source.is_dir() else source.with_suffix('.html')
            plots=getattr(args,'plots',None)
            if plots is None and args.command=='report': plots='default'
            cache_directory=destination.parent/'.research-plot-cache' if destination is not None else None
            data = report.compare(args.paths,args.candidate, plots=plots,cache_directory=cache_directory,
                notify=lambda message: print(message,file=sys.stderr,flush=True))
            if destination is not None: report.render_html(data,destination)
            if args.json:
                print(json.dumps(data,indent=2))
            else:
                print(report.text_report(data, list_candidates=args.list_candidates),end='')
                if destination is not None: print(f'HTML: {destination}')
    except KeyboardInterrupt:
        parser.exit(130, 'Interrupted. Campaign pending state is preserved; use resume. Standalone evaluations must be requested again.\n')
    except ImportError as error:
        parser.exit(2, f'Research dependency error: {error}. Interactive views require dada-engine-solver[plot].\n')
    except (ValueError, OSError, KeyError, TypeError, RuntimeError) as error:
        parser.exit(2, f'Research error: {error}\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

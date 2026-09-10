"""Run or resume a persistent Sobol campaign for an approximate time budget."""
import argparse
from pathlib import Path
from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.runner import OptimizationCampaign, parse_budget
from dada_solver.campaign.report import readable_report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path, help='Campaign TOML for creation, or existing directory for resume')
    parser.add_argument('--directory', type=Path, help='Destination for a new campaign')
    parser.add_argument('--budget', required=True, help='For example 30m, 1h30m, 45s')
    parser.add_argument('--max-candidates', type=int, help='Optional per-run cap for smoke work')
    parser.add_argument('--retry-incomplete', action='store_true',
        help='Retry the most recent deadline-interrupted candidate before advancing')
    args = parser.parse_args(argv)
    try:
        if args.campaign.is_dir():
            if args.directory: parser.error('Do not supply --directory when resuming a directory.')
            campaign = OptimizationCampaign.resume(args.campaign)
        else:
            definition = CampaignDefinition(args.campaign)
            destination = args.directory or Path('outputs')/(args.campaign.stem+'_campaign')
            campaign = OptimizationCampaign(definition, destination)
        report = campaign.run(parse_budget(args.budget), maximum_candidates=args.max_candidates,
                              retry_incomplete=args.retry_incomplete)
    except (ValueError, RuntimeError) as error:
        parser.exit(2, f'Campaign error: {error}\n')
    print(readable_report(report))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

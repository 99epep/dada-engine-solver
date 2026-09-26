"""Capture one original thermo-5D evaluation and its actual local read inputs.

This audit-only script runs no search, edits no original result, and does not
implement a Research adapter. Run from the checkout:
MPLCONFIGDIR=/tmp/dada-matplotlib PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:examples \
    python3 docs/research_audit/capture_reference.py
"""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
DEST = Path(__file__).resolve().parent


def main():
    if Path.cwd() != ROOT:
        raise RuntimeError("Run this reference capture from the checkout root.")
    reads = set()

    def audit(event, args):
        if event != "open" or not isinstance(args[0], (str, bytes)):
            return
        path = Path(args[0]).resolve()
        if path.is_relative_to(ROOT) and path.suffix in (".json", ".jsonl", ".toml"):
            if args[1] is not None and "r" in args[1]:
                reads.add(path)

    sys.addaudithook(audit)
    import optimize_sixbar_pairs_thermo5d_3952_hlat25 as legacy
    from dada_solver.campaign.definition import runtime_identity
    from dada_solver.wall_backend import backend_identity

    start = time.monotonic()
    seed, definition, mass, geometry, parameters, identity = legacy.make_3952_basis()
    pair = legacy.load_pair(1)
    report_path = ROOT / "outputs/sixbar_thermo5d_3952/rank_01/report.json"
    report = json.loads(report_path.read_text())
    selected = report["best_feasible"]
    history = legacy.load_history(report_path.with_name("history.jsonl"))
    source = next(r for r in history if r["candidate_id"] == selected["warm_start_source"])
    design = legacy.build_design(seed, geometry, mass, selected["parameters"], pair["small"], pair["large"])
    hardware = legacy.hardware(design)
    initial = legacy.warm_state(source["last_complete_state"], source["hardware"], hardware)
    setup_seconds = time.monotonic() - start
    result = legacy.evaluate(label="research_audit_rank01_reference", design=design,
                             definition=definition, initial=initial,
                             deadline=time.monotonic() + 180, candidate_seconds=180)
    observed = dict(result["compact"])
    fields = ("indicated_thermal_efficiency", "indicated_power_w", "heat_input_w",
              "maximum_pressure_pa", "maximum_temperature_k", "maximum_absolute_mass_flow_kg_s")
    differences = {k: observed[k] - selected["result"][k] for k in fields
                   if observed.get(k) is not None and selected["result"].get(k) is not None}
    inputs = []
    tracked = set(subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0"))
    for path in sorted(reads):
        rel = str(path.relative_to(ROOT))
        inputs.append(dict(path=rel, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                           bytes=path.stat().st_size, tracked=rel in tracked))
    inventory = json.loads((DEST / "example_inventory.json").read_text())
    script = "examples/optimize_sixbar_pairs_thermo5d_3952_hlat25.py"
    row = next(r for r in inventory["examples"] if r["file"] == script)
    source_paths = [script, *row["transitive_example_imports"]]
    snapshot = dict(schema_version=1, purpose="Original evaluator replay; not Research parity or a search",
                    runtime=runtime_identity(), backend=backend_identity(definition.wall_backend),
                    source_commit=inventory["source_commit"], inputs=inputs,
                    example_sources=[dict(path=p, sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()) for p in source_paths],
                    selected_candidate=selected, selected_pair=json.loads(pair["pair_path"].read_text()),
                    source_initial_candidate_id=source["candidate_id"], initial_state=initial.tolist(),
                    geometry=geometry, seed_parameters=parameters, total_gas_mass_kg=mass,
                    resolved_configuration=asdict(design.configuration),
                    resolved_heat_in=asdict(design.heat_in), resolved_heat_out=asdict(design.heat_out),
                    wall_settings=asdict(definition.wall_numerical_settings),
                    fresh_replay=dict(compact=observed, feasible=result["feasible"], reasons=result["reasons"],
                                      safe_retry=result["safe_retry"], backend=result["backend"],
                                      final_state=result["state"].tolist() if result["state"] is not None else None,
                                      setup_seconds=setup_seconds, evaluation_seconds=result["elapsed_seconds"],
                                      differences_from_stored_result=differences))
    (DEST / "rank01_reference_snapshot.json").write_text(json.dumps(snapshot, indent=2, allow_nan=False) + "\n")
    print(json.dumps(snapshot["fresh_replay"], indent=2))


if __name__ == "__main__":
    main()

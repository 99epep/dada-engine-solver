"""Regenerate the local, static migration inventory without importing examples.

Run from any directory with Python 3.11+. This is an audit utility, not a
research implementation. AST evidence is not a claim of runtime parity.
"""
from __future__ import annotations

import ast
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tomllib

ROOT = Path(__file__).resolve().parents[2]
DEST = Path(__file__).resolve().parent
ARTIFACT_SUFFIXES = (".toml", ".json", ".jsonl", ".csv", ".npz", ".npy",
                     ".png", ".svg", ".gif", ".mp4", ".pdf", ".txt")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def classify(path, data, calls):
    """Provisional migration destination; evidence is stored separately."""
    name = path.stem
    if path.suffix == ".toml":
        if "campaign" in data:
            return "campaign preset", "campaign.CampaignDefinition (retain existing preset)"
        if "air_sizing" in data or "tube_bank" in data or ("geometry" in data and "properties" in data):
            return "hardware configuration", "exchangers.hardware.load_hardware_definition"
        return "model or sizing configuration", "existing configuration loaders (retain schema)"
    if name.startswith(("benchmark_", "profile_", "check_numba", "check_solver", "check_stage5",
                        "check_wall", "compare_numba", "compare_solver", "compare_stage5",
                        "diagnose_periodic", "audit_wall", "summarize_solver", "study_stage5",
                        "numba_wall", "experimental_periodic")) or name == "stage5_mechanisms":
        return "developer benchmark", "developer tooling; no research strategy"
    if name.startswith("test_sixbar_robustness"):
        return "robustness", "research tolerance protocol (deferred)"
    if name.startswith(("animate_", "draw_", "visualize_")):
        return "presentation / layout history", "shared mechanism frames + optional layout preset (deferred)"
    if name.startswith(("plot_", "export_", "report_")):
        return "diagnostic export / presentation", "research diagnostics + read-only renderer (deferred)"
    if name.startswith("fit_"):
        return "motion fitting", "motion fitting evaluator (deferred)"
    if name.startswith(("search_", "polish_", "mirror_polish_")):
        return "mechanism synthesis", "mechanical screen / staged synthesis protocol (deferred)"
    if "sixbar_pairs_thermo5d" in name:
        return "six-bar thermo-5D", "first research slice: fixed pair + five thermal coordinates"
    if "sixbar_pairs_thermo" in name or "rank01_compact" in name:
        return "coupled mechanism refinement", "coupled evaluator / separate refinement protocol (deferred)"
    if "temperature" in name:
        return "temperature study", "temperature sweep + child campaigns (deferred)"
    if "valve" in name:
        return "valve topology study", "enumerated topology protocol (deferred)"
    if name.startswith(("optimize_", "refine_")):
        return "motion / hardware search", "family adapter + explicit search protocol (deferred)"
    if any(x in name for x in ("hardware", "air_wall", "doty", "microtube", "parallel", "trials")):
        return "hardware / exchanger trial", "existing hardware and wall evaluators; named study policies"
    if name.startswith(("evaluate_", "scan_", "compare_", "analyze_", "diagnose_", "check_")) or name.endswith("_comparison"):
        return "evaluation / diagnostic study", "evaluation or sweep protocol; preserve study limits"
    if name == "pluggable_kinematics_smoke":
        return "developer smoke", "retain integration regression"
    return "mechanism geometry helper", "six_bar / four_bar geometry after parity review"


def main():
    tracked = set(subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0"))
    paths = sorted(p for p in (ROOT / "examples").rglob("*") if p.suffix in (".py", ".toml"))
    modules = {p.stem: str(p.relative_to(ROOT)) for p in paths if p.suffix == ".py"}
    docs = {str(p.relative_to(ROOT)): p.read_text() for p in (ROOT / "docs").rglob("*.md")
            if "DADA_ENGINE_RESEARCH" not in p.name}
    tests = {str(p.relative_to(ROOT)): p.read_text() for p in (ROOT / "tests").rglob("*.py")}
    rows = []
    for path in paths:
        rel = str(path.relative_to(ROOT))
        source = path.read_text()
        evidence = dict(imports=[], local_imports=[], cli=[], inputs=[], outputs=[],
                        algorithm=[], feasibility=[], rendering=[], path_bindings=[],
                        functions=[], artifact_literals=[], top_level_effects=[])
        data = {}
        if path.suffix == ".py":
            tree = ast.parse(source, filename=rel)
            evidence["description"] = ast.get_docstring(tree) or ""
            def snippet(node):
                return f"L{node.lineno}: {ast.unparse(node)}"
            for node in ast.walk(tree):
                if isinstance(node, (ast.Import, ast.ImportFrom)):
                    evidence["imports"].append(snippet(node))
                    names = [a.name for a in node.names] if isinstance(node, ast.Import) else [node.module or ""]
                    for name in names:
                        key = name.removeprefix("examples.").split(".")[0]
                        if key in modules:
                            evidence["local_imports"].append(modules[key])
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    evidence["functions"].append(f"{node.name}:L{node.lineno}")
                    if re.search(r"feasib|valid|admiss|constraint|eligible|screen|closure", node.name):
                        evidence["feasibility"].append(snippet(node))
                if isinstance(node, ast.Call):
                    function = ast.unparse(node.func)
                    leaf = function.split(".")[-1]
                    if leaf == "add_argument": evidence["cli"].append(snippet(node))
                    if re.search(r"read_|load|genfromtxt|glob", leaf): evidence["inputs"].append(snippet(node))
                    if re.search(r"write|save|dump|append_record|append_jsonl", leaf): evidence["outputs"].append(snippet(node))
                    if leaf == "open":
                        evidence["inputs"].append(snippet(node))
                        evidence["outputs"].append(snippet(node))
                    if re.search(r"Sobol|LatinHypercube|differential_evolution|minimize|least_squares|curve_fit|random_base2|next_point", leaf):
                        evidence["algorithm"].append(snippet(node))
                    if re.search(r"plot|imshow|scatter|FuncAnimation|PillowWriter|FFMpegWriter|savefig|svg", leaf):
                        evidence["rendering"].append(snippet(node))
                if isinstance(node, (ast.Assign, ast.AnnAssign)):
                    expression = ast.unparse(node)
                    if any(x in expression for x in ("Path(", "Path.cwd", "ROOT", "OUTPUT", "SOURCE", "REPORT")):
                        # Avoid duplicating large dictionaries and plotting code.
                        evidence["path_bindings"].append(snippet(node))
                if isinstance(node, ast.Constant) and isinstance(node.value, str):
                    value = node.value
                    if "\n" not in value and (value.endswith(ARTIFACT_SUFFIXES) or value.startswith(("outputs/", "examples/"))):
                        evidence["artifact_literals"].append(f"L{node.lineno}: {value}")
            for node in tree.body:
                if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
                    evidence["top_level_effects"].append(snippet(node))
        else:
            data = tomllib.loads(source)
            evidence["description"] = "TOML tables: " + ", ".join(data)
            evidence["toml"] = data
            for number, line in enumerate(source.splitlines(), 1):
                if re.search(r"configuration\s*=|path\s*=|file\s*=", line):
                    evidence["inputs"].append(f"L{number}: {line.strip()}")
                if re.search(r"objective|constraint|limit|valid", line):
                    evidence["feasibility"].append(f"L{number}: {line.strip()}")
        category, successor = classify(path, data, evidence["algorithm"])
        for key in evidence:
            if isinstance(evidence[key], list): evidence[key] = sorted(set(evidence[key]))
        production = [x for x in evidence["imports"] if "dada_solver" in x]
        rows.append(dict(file=rel, sha256=sha(path), tracked=rel in tracked,
                         category=category, canonical_successor=successor,
                         shared_logic=production + evidence["local_imports"],
                         study_specific_logic=evidence["description"],
                         inputs=evidence["inputs"], outputs=evidence["outputs"],
                         tests=[p for p, s in tests.items() if path.name in s or re.search(r"\b" + re.escape(path.stem) + r"\b", s)],
                         documentation=[p for p, s in docs.items() if path.name in s],
                         migration_status="retained; static inventory verified; semantic parity pending",
                         evidence=evidence))
    graph = {r["file"]: r["evidence"]["local_imports"] for r in rows}
    # Reviewed semantic status is explicit, never inferred from an AST import.
    reviewed=json.loads((DEST/'kinematics_v2_migrations.json').read_text())
    for row in rows:
        if row['file'] in reviewed['examples']:
            row.update(reviewed['examples'][row['file']])
            row['migration_evidence']=reviewed['evidence']
    for row in rows:
        seen = set()
        todo = list(graph[row["file"]])
        while todo:
            p = todo.pop()
            if p not in seen:
                seen.add(p)
                todo.extend(graph.get(p, []))
        row["transitive_example_imports"] = sorted(seen - {row["file"]})
    payload = dict(schema_version=1,
                   method="Full-source AST/TOML inspection; no example execution. Categories are proposed destinations; dynamic paths remain expressions.",
                   source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip(),
                   counts=dict(total=len(rows), python=sum(p.suffix == ".py" for p in paths),
                               toml=sum(p.suffix == ".toml" for p in paths),
                               untracked=sum(not r["tracked"] for r in rows)),
                   categories=dict(sorted(Counter(r["category"] for r in rows).items())), examples=rows)
    (DEST / "example_inventory.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    columns = ["file", "category", "canonical_successor", "shared_logic", "study_specific_logic",
               "inputs", "outputs", "tests", "migration_status", "sha256", "tracked", "documentation"]
    with (DEST.parent / "DADA_ENGINE_RESEARCH_MIGRATION_MATRIX.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(row[k], ensure_ascii=False) if isinstance(row[k], list) else row[k] for k in columns})
    print(json.dumps({"counts": payload["counts"], "categories": payload["categories"]}, indent=2))


if __name__ == "__main__":
    main()

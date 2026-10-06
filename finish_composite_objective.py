#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path.cwd()


def replace_once(text, old, new, label):
    if new in text:
        return text
    if old not in text:
        raise SystemExit(f"Cannot locate {label}; no file was written.")
    return text.replace(old, new, 1)


def patch_report():
    path = ROOT / "src/dada_solver/research/report.py"
    text = path.read_text()

    pattern = re.compile(
        r"def _add_microtube_productivity\(record\):\n"
        r".*?(?=\n\ndef select_records\()",
        re.S,
    )
    m = pattern.search(text)
    if not m:
        raise SystemExit("Cannot locate report._add_microtube_productivity; no file was written.")
    if "cooling_cop_times_power_per_total_microtube_w" not in m.group(0):
        helper = '''def _add_microtube_productivity(record):
    """Derive missing microtube metrics offline, never rewriting storage."""
    from dada_solver.campaign.objectives import (
        cooling_power_per_total_microtube, cooling_cop_times_power_per_total_microtube)
    metrics = record.setdefault('metrics', {})
    physical = record.get('resolved_parameters', {})
    derived = record.setdefault('derived', {})
    if 'cooling_power_per_total_microtube_w' not in metrics:
        counts = [derived.get(f'{side}_microtube_count',
                              physical.get(f'microtube.{side}.tube_count'))
                  for side in ('heat_in', 'heat_out')]
        if all(type(n) is int and n > 0 for n in counts):
            derived.update(heat_in_microtube_count=counts[0], heat_out_microtube_count=counts[1],
                           total_microtube_count=sum(counts))
            refrigeration = (metrics.get('operating_mode') == 'refrigeration'
                or (metrics.get('operating_mode') is None and metrics.get('cooling_cop') is not None))
            metrics['cooling_power_per_total_microtube_w'] = cooling_power_per_total_microtube(
                metrics.get('cooling_power_w') if refrigeration else None, *counts)
            derived['microtube_productivity_source'] = 'offline_stored_values'
    metrics.setdefault('cooling_cop_times_power_per_total_microtube_w',
        cooling_cop_times_power_per_total_microtube(
            metrics.get('cooling_cop'), metrics.get('cooling_power_per_total_microtube_w')))
'''
        text = text[:m.start()] + helper + text[m.end():]

    marker = """        if m.get('cooling_power_per_total_microtube_w') is not None:
            performance += f\"; Qcold/microtube={m['cooling_power_per_total_microtube_w']} W/microtube\"
"""
    addition = marker + """        if m.get('cooling_cop_times_power_per_total_microtube_w') is not None:
            performance += (f\"; COP×Qcold/microtube=\"
                f\"{m['cooling_cop_times_power_per_total_microtube_w']} W/microtube\")
"""
    text = replace_once(text, marker, addition, "report text composite metric")

    text = replace_once(
        text,
        "for k in ('cooling_cop','indicated_thermal_efficiency','cooling_power_per_total_microtube_w')",
        "for k in ('cooling_cop','indicated_thermal_efficiency',\n                                'cooling_power_per_total_microtube_w',\n                                'cooling_cop_times_power_per_total_microtube_w')",
        "report compact metric list",
    )

    old_js = """const productivity=d.scientific.objective.type==='maximize_cooling_power_per_total_microtube';
const ykey=productivity?'cooling_power_per_total_microtube_w':cooling?'cooling_cop':'indicated_thermal_efficiency',yscale=cooling?1:100,ylabel=productivity?'Cooling power per microtube [W/microtube]':cooling?'Cooling COP [1]':'Indicated efficiency [%]';"""
    new_js = """const objective=d.scientific.objective.type;
const productivity=objective==='maximize_cooling_power_per_total_microtube';
const composite=objective==='maximize_cooling_cop_times_power_per_total_microtube';
const ykey=composite?'cooling_cop_times_power_per_total_microtube_w':productivity?'cooling_power_per_total_microtube_w':cooling?'cooling_cop':'indicated_thermal_efficiency',yscale=cooling?1:100,ylabel=composite?'COP×Qcold/microtube [W/microtube]':productivity?'Cooling power per microtube [W/microtube]':cooling?'Cooling COP [1]':'Indicated efficiency [%]';"""
    text = replace_once(text, old_js, new_js, "report objective display selector")

    text = replace_once(
        text,
        "if(cooling){byId('performanceTitle').textContent=productivity?'Cooling power per microtube':'Indicated power and cooling COP';",
        "if(cooling){byId('performanceTitle').textContent=composite?'COP×Qcold/microtube':productivity?'Cooling power per microtube':'Indicated power and cooling COP';",
        "report performance title",
    )

    path.write_text(text)
    print(f"patched {path}")


def patch_tests():
    path = ROOT / "tests/test_microtube_objective.py"
    text = path.read_text()

    text = replace_once(
        text,
        """from dada_solver.campaign.objectives import (
    MaximizeCoolingPowerPerTotalMicrotube, cooling_power_per_total_microtube,
)""",
        """from dada_solver.campaign.objectives import (
    MaximizeCoolingPowerPerTotalMicrotube,
    MaximizeCoolingCopTimesPowerPerTotalMicrotube,
    cooling_power_per_total_microtube, cooling_cop_times_power_per_total_microtube,
)""",
        "microtube-objective imports",
    )

    if "COMPOSITE_NAME = 'maximize_cooling_cop_times_power_per_total_microtube'" not in text:
        m = re.search(r"^NAME\s*=\s*['\"]maximize_cooling_power_per_total_microtube['\"].*$", text, re.M)
        if m:
            text = text[:m.end()] + "\nCOMPOSITE_NAME = 'maximize_cooling_cop_times_power_per_total_microtube'" + text[m.end():]
        else:
            anchor = "from dada_solver.research import report\n"
            if anchor not in text:
                raise SystemExit("Cannot locate safe insertion point for COMPOSITE_NAME.")
            text = text.replace(anchor, anchor + "\nCOMPOSITE_NAME = 'maximize_cooling_cop_times_power_per_total_microtube'\n", 1)

    if "def test_composite_score_objective_unavailability_and_ranking():" not in text:
        anchor = "@pytest.mark.parametrize('objective,unit', OBJECTIVE_UNITS.items())"
        idx = text.find(anchor)
        if idx < 0:
            raise SystemExit("Cannot locate schema objective parametrization.")
        block = '''def test_composite_score_objective_unavailability_and_ranking():
    assert cooling_cop_times_power_per_total_microtube(2.0, 1.5) == 3.0
    assert cooling_cop_times_power_per_total_microtube(None, 1.5) is None
    assert cooling_cop_times_power_per_total_microtube(2.0, None) is None
    objective = MaximizeCoolingCopTimesPowerPerTotalMicrotube()
    evaluation = SimpleNamespace(performance=SimpleNamespace(
        operating_mode=OperatingMode.REFRIGERATION, cooling_power=450., cooling_cop=2.))
    result = objective.evaluate(
        evaluation, heat_in_microtube_count=100, heat_out_microtube_count=200)
    assert result.available and result.value == -3.0
    evaluation.performance.cooling_cop = None
    assert not objective.evaluate(
        evaluation, heat_in_microtube_count=100, heat_out_microtube_count=200).available
    evaluation.performance.cooling_cop = 2.
    assert not objective.evaluate(evaluation).available

    from dataclasses import asdict
    from dada_solver.campaign.report import elite_records
    def score(candidate_id, cop, productivity):
        candidate = SimpleNamespace(performance=SimpleNamespace(
            operating_mode=OperatingMode.REFRIGERATION,
            cooling_power=productivity*100, cooling_cop=cop))
        value = objective.evaluate(
            candidate, heat_in_microtube_count=50, heat_out_microtube_count=50)
        return dict(candidate_id=candidate_id, status='feasible', objective=asdict(value))

    ranked = elite_records([
        score('A', 1.5, 0.40),
        score('B', 2.2, 0.30),
        score('C', 3.0, 0.15),
    ])
    assert ranked[0]['candidate_id'] == 'B'


'''
        text = text[:idx] + block + text[idx:]

    text = re.sub(
        r"def test_motor_rejects_productivity\(tmp_path\):\n"
        r"(\s+)path = initialize_v3\(tmp_path/'study\.toml', mode='motor'\)\n"
        r"\1raw = tomllib\.loads\(path\.read_text\(\)\)\n"
        r"\1raw\['objective'\] = dict\(type=NAME, unit='W/microtube'\)",
        r"@pytest.mark.parametrize('objective', [NAME, COMPOSITE_NAME])\n"
        r"def test_motor_rejects_productivity(tmp_path, objective):\n"
        r"\1path = initialize_v3(tmp_path/'study.toml', mode='motor')\n"
        r"\1raw = tomllib.loads(path.read_text())\n"
        r"\1raw['objective'] = dict(type=objective, unit='W/microtube')",
        text,
        count=1,
    )

    text = text.replace(
        "@pytest.mark.parametrize('objective', [NAME, 'maximize_cooling_cop'])",
        "@pytest.mark.parametrize('objective', [NAME, COMPOSITE_NAME, 'maximize_cooling_cop'])",
        1,
    )

    text = replace_once(
        text,
        """    assert result['metrics']['cooling_power_per_total_microtube_w'] == 1.5
    assert result['objective']['value'] == (-1.5 if objective == NAME else -2.)""",
        """    assert result['metrics']['cooling_power_per_total_microtube_w'] == 1.5
    assert result['metrics']['cooling_cop_times_power_per_total_microtube_w'] == (
        result['metrics']['cooling_cop']
        * result['metrics']['cooling_power_per_total_microtube_w'])
    expected = {NAME: -1.5, COMPOSITE_NAME: -3.0, 'maximize_cooling_cop': -2.0}
    assert result['objective']['value'] == expected[objective]""",
        "built-count composite assertions",
    )

    marker = "    assert 'Cooling power per microtube [W/microtube]' in page\n"
    if "tmp_path/'composite.html'" not in text:
        addition = marker + """    composite_page = report.render_html(dict(name='Composite', records=[], best=[],
        scientific=dict(objective=dict(type=COMPOSITE_NAME, unit='W/microtube'))),
        tmp_path/'composite.html').read_text()
    assert '\"cooling_objective\": true' in composite_page
    assert "composite=objective==='maximize_cooling_cop_times_power_per_total_microtube'" in composite_page
    assert 'COP×Qcold/microtube [W/microtube]' in composite_page
"""
        text = replace_once(text, marker, addition, "composite report assertions")

    text = replace_once(
        text,
        """    parts = metric_parts(dict(cooling_power_per_total_microtube_w=1.5,
                             cooling_cop=2., cooling_power_w=450.))""",
        """    parts = metric_parts(dict(
                             cooling_cop_times_power_per_total_microtube_w=3.0,
                             cooling_power_per_total_microtube_w=1.5,
                             cooling_cop=2., cooling_power_w=450.))""",
        "composite progress metrics",
    )

    if "assert 'COP×Qcold/microtube 3 W/microtube' in parts" not in text:
        text = text.replace(
            "    assert 'Qcold/microtube 1.5 W/microtube' in parts\n",
            "    assert 'COP×Qcold/microtube 3 W/microtube' in parts\n"
            "    assert 'Qcold/microtube 1.5 W/microtube' in parts\n",
            1,
        )

    old_data = "                    cooling_cop=2., cooling_power_w=450., cooling_power_per_total_microtube_w=1.5))],"
    new_data = """                    cooling_cop=2., cooling_power_w=450.,
                    cooling_power_per_total_microtube_w=1.5,
                    cooling_cop_times_power_per_total_microtube_w=3.0))],"""
    text = replace_once(text, old_data, new_data, "text-report composite fixture")

    if "assert 'COP×Qcold/microtube=3.0 W/microtube' in report.text_report(data)" not in text:
        text = text.replace(
            "    assert 'Qcold/microtube=1.5 W/microtube' in report.text_report(data)\n",
            "    assert 'Qcold/microtube=1.5 W/microtube' in report.text_report(data)\n"
            "    assert 'COP×Qcold/microtube=3.0 W/microtube' in report.text_report(data)\n",
            1,
        )

    if "assert inspected['metrics']['cooling_cop_times_power_per_total_microtube_w']" not in text:
        text = text.replace(
            "    assert inspected['metrics']['cooling_power_per_total_microtube_w'] == 450/total\n",
            "    assert inspected['metrics']['cooling_power_per_total_microtube_w'] == 450/total\n"
            "    assert inspected['metrics']['cooling_cop_times_power_per_total_microtube_w'] == 2*450/total\n",
            1,
        )

    path.write_text(text)
    print(f"patched {path}")


patch_report()
patch_tests()

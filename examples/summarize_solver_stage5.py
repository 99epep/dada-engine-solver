"""Separate exact production timings from interpolation experiments."""
import json
from pathlib import Path
import statistics

ROOT=Path(__file__).resolve().parents[1]/'outputs/solver_acceleration_stage5'


def summarize():
    result={}
    for name in ('p1_baseline','p1_cache','p2_shared','p1_sixbar_baseline','p1_sixbar_cache','p2_sixbar','profile'):
        directory=ROOT/name
        if not directory.exists():continue
        cases={}
        for p in directory.glob('*_*.json'):
            row=json.loads(p.read_text())
            if 'case' not in row:continue
            cases.setdefault(row['case'],[]).append(row)
        result[name]={}
        for case,rows in cases.items():
            sections={k:statistics.median(sum(r['sections_seconds'].get(k,[])) for r in rows)
                      for k in {k for r in rows for k in r['sections_seconds']}}
            first=rows[0]
            result[name][case]=dict(repeats=len(rows),median_seconds=statistics.median(r['elapsed_seconds'] for r in rows),
                maximum_seconds=max(r['elapsed_seconds'] for r in rows),sections_median_seconds=sections,
                kinematics_seconds=statistics.median(r['backend_statistics'].get('kinematics_seconds') or 0 for r in rows),
                exact_cache=first['backend_statistics'].get('exact_kinematics_cache'),
                report_materialization_and_other_seconds=statistics.median(r['elapsed_seconds']-r['candidate_seconds']-
                    sum(sum(v) for v in r['sections_seconds'].values()) for r in rows),
                json_serialization_seconds=statistics.median(r['report_serialization_seconds'] for r in rows))
    return result

if __name__=='__main__':print(json.dumps(summarize(),indent=2))

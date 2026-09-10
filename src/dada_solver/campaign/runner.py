"""Budgeted sequential orchestration with durable resume and feasible archives."""
from datetime import datetime, timezone
from pathlib import Path
import math
import re
import time
import numpy as np
from dada_solver.campaign.candidate import Candidate
from dada_solver.campaign.definition import CampaignDefinition
from dada_solver.campaign.evaluator import MachineEvaluator
from dada_solver.campaign.history import CampaignHistory, atomic_json
from dada_solver.campaign.report import elite_records, make_report, readable_report
from dada_solver.campaign.strategy import SobolStrategy


def parse_budget(text):
    text = str(text).lower().replace(' ', '')
    if re.fullmatch(r'\d+(?:\.\d+)?', text): return float(text)
    pieces = re.findall(r'(\d+(?:\.\d+)?)([hms])', text)
    if not pieces or ''.join(number+unit for number,unit in pieces) != text:
        raise ValueError('Budget must be seconds or a duration such as 30m, 1h30m or 45s.')
    return sum(float(number)*{'h':3600,'m':60,'s':1}[unit] for number,unit in pieces)


def estimated_next_seconds(records, initial_seconds):
    durations = [r['duration_seconds'] for r in records if r['integrated'] and not r.get('cache_hit')][-10:]
    return max(1e-6, float(np.percentile(durations, 75))*1.1) if durations else initial_seconds


class OptimizationCampaign:
    def __init__(self, definition, directory, *, evaluator=None, clock=time.monotonic):
        self.definition, self.directory = definition, Path(directory)
        self.history = CampaignHistory(self.directory)
        self.evaluator = evaluator if evaluator is not None else MachineEvaluator(definition)
        self.clock = clock
        with self.history.locked():
            manifest = self.directory/'definition.json'
            if manifest.exists():
                import json
                if json.loads(manifest.read_text())['definition_id'] != definition.definition_id:
                    raise ValueError('Campaign identity changed; use a new directory.')
            else:
                if self.history.path.exists(): raise ValueError('History exists without a campaign definition.')
                (self.directory/'campaign.toml').write_text(definition.source)
                (self.directory/'base.toml').write_text(definition.base_source)
                atomic_json(manifest, dict(definition_id=definition.definition_id, **definition.identity))

    @classmethod
    def resume(cls, directory, *, evaluator=None, clock=time.monotonic):
        return cls(CampaignDefinition.resume(directory), directory, evaluator=evaluator, clock=clock)

    def run(self, budget_seconds, *, maximum_candidates=None):
        if not math.isfinite(budget_seconds) or budget_seconds < 0:
            raise ValueError('Budget must be finite and nonnegative.')
        limit = self.definition.maximum_candidates if maximum_candidates is None else maximum_candidates
        if limit is not None and (not isinstance(limit, int) or limit < 0):
            raise ValueError('Candidate limit must be a nonnegative integer.')
        with self.history.locked():
            return self._run_locked(float(budget_seconds), limit)

    def _run_locked(self, budget, limit):
        started = self.clock()
        records = self.history.load()
        before = list(records)
        state = self.history.state()
        phase_id = max(state.get('phase_counter',0), max((r['phase_id'] for r in records), default=0))+1
        completed_indices = {r['sequence_index'] for r in records}
        pending = state.get('pending')
        if pending and pending['sequence_index'] in completed_indices:
            pending = None
        consumed = max((r['sequence_index']+1 for r in records), default=0)
        index = max(consumed, state.get('search',{}).get('index',0))
        strategy = SobolStrategy(len(self.definition.space.parameters), seed=self.definition.seed,
                                 scramble=self.definition.scramble, index=index)
        cache = {r['candidate_id']:r for r in records if not r.get('cache_hit')}
        phase = []
        def save_state():
            self.history.save_state(dict(schema_version=1, phase_counter=phase_id,
                search=strategy.state(), pending=pending,
                best_feasible_ids=[r['candidate_id'] for r in elite_records(records,self.definition.elite_size)]))
        save_state()
        while limit is None or len(phase)<limit:
            remaining = budget-(self.clock()-started)
            estimate = estimated_next_seconds(records, self.definition.initial_evaluation_seconds)
            if remaining <= 0 or remaining < estimate:
                break
            if pending:
                candidate = Candidate(pending['payload_json'], pending['candidate_id'])
                sequence_index = pending['sequence_index']
            else:
                sequence_index = strategy.index
                candidate = Candidate.create(self.definition.space, strategy.next_point(),
                    families=self.definition.families, numerical_settings=self.definition.numerical_settings,
                    definition_id=self.definition.definition_id)
                pending = dict(candidate_id=candidate.candidate_id, payload_json=candidate.payload_json,
                               sequence_index=sequence_index)
                # Record the in-flight coordinate before integration. Resume retries
                # only this point if no durable completed result was written.
                save_state()
            evaluation_started = self.clock()
            cached = cache.get(candidate.candidate_id)
            if cached is not None:
                reserved = {'candidate_id','timestamp','evaluation_number','sequence_index','phase_id',
                    'duration_seconds','cache_hit','cache_source_evaluation',*candidate.payload.keys()}
                result = {k:v for k,v in cached.items() if k not in reserved}
            else:
                result = self.evaluator.evaluate(candidate)
            record = dict(**candidate.payload, **result,
                candidate_id=candidate.candidate_id, timestamp=datetime.now(timezone.utc).isoformat(),
                evaluation_number=len(records), sequence_index=sequence_index, phase_id=phase_id,
                duration_seconds=self.clock()-evaluation_started, cache_hit=cached is not None,
                cache_source_evaluation=None if cached is None else cached['evaluation_number'])
            self.history.save(record)
            records.append(record); phase.append(record)
            if cached is None: cache[candidate.candidate_id] = record
            pending = None
            save_state()
        report = make_report(self.definition.space, before, phase, requested_seconds=budget,
            elapsed_seconds=self.clock()-started, elite_size=self.definition.elite_size)
        report.update(phase_id=phase_id, next_sequence_index=strategy.index,
            stopping_reason='candidate_limit' if limit is not None and len(phase)>=limit else 'wall_clock_budget',
            estimated_next_evaluation_seconds=estimated_next_seconds(records,self.definition.initial_evaluation_seconds))
        reports = self.directory/'reports'; reports.mkdir(exist_ok=True)
        atomic_json(reports/f'phase_{phase_id:04d}.json', report)
        (reports/f'phase_{phase_id:04d}.txt').write_text(readable_report(report))
        atomic_json(self.directory/'report.json', report)
        (self.directory/'report.txt').write_text(readable_report(report))
        return report

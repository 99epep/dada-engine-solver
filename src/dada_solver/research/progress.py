"""Small stdout progress sink, independent of reports and scientific records."""
import sys


def duration(seconds):
    return f'{seconds/60:g}m' if seconds and seconds % 60 == 0 else f'{seconds:.0f}s'


def best_text(best):
    if not best: return 'best: unavailable'
    objective = best.get('objective') or {}
    parts = [f"best {best['candidate_id'][:12]} objective={objective.get('value')}"]
    m = best.get('metrics') or {}
    for key, label in [('cooling_cop','COP'), ('cooling_power_w','Qcold [W]'),
                       ('indicated_mechanical_input_power_w','indicated Pin [W]')]:
        if m.get(key) is not None: parts.append(f'{label}={m[key]:.6g}')
    return ' · '.join(parts)


class CLIProgress:
    def __init__(self, stream=None):
        self.stream = stream if stream is not None else sys.stdout
        self.tty = self.stream.isatty()

    def __call__(self, event):
        limit = event['maximum_candidates'] if event['maximum_candidates'] is not None else 'unlimited'
        if event['event'] == 'start':
            text = (f"Phase {event['phase_id']} · {'Search step' if event.get('scheduled_search') else 'Sobol index'} {event['sobol_index']} · "
                    f"budget {duration(event['budget_seconds'])} · max {limit} new candidates")
            print(text, file=self.stream, flush=True)
            return
        failures = ', '.join(f'{k}={v}' for k,v in sorted(event['failure_counts'].items(), key=lambda x:(-x[1],x[0]))[:3]) or 'none'
        text = (f"{'Done · ' if event['event']=='finish' else ''}{event['attempted']}/{limit} attempts · "
                f"{duration(event['elapsed_seconds'])}/{duration(event['budget_seconds'])} · "
                f"converged {event['converged']} · feasible {event['feasible']} · "
                f"failures: {failures} · {best_text(event['best'])}")
        print(('\r\033[K' if self.tty else '')+text, file=self.stream,
              end='\n' if not self.tty or event['event']=='finish' else '', flush=True)

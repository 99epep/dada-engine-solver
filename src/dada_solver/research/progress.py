"""Terminal-only presentation: durable evaluation lines and one transient status."""
import math
import os
import re
import shutil
import sys
import unicodedata


def duration(seconds):
    seconds = max(0, round(seconds))
    minutes, seconds = divmod(seconds, 60)
    return f'{minutes}m{seconds:02d}s' if minutes and seconds else f'{minutes}m' if minutes else f'{seconds}s'


def clean(text):
    # Scientific messages are evidence, not terminal control sequences.
    text = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', str(text))
    return ' '.join(''.join(c if unicodedata.category(c)[0] != 'C' else ' ' for c in text).split())


def cells(text):
    return sum(0 if unicodedata.combining(c) else 2 if unicodedata.east_asian_width(c) in ('W','F') else 1 for c in text)


def clipped(text, width):
    if cells(text) <= width: return text
    out = ''
    for char in text:
        if cells(out+char) > max(0,width-1): break
        out += char
    return out+'…' if width else ''


def metric_parts(metrics, *, champion=False):
    m = metrics or {}
    precision = '.6g' if champion else '.5g'
    fields = ([('cooling_cop','COP',precision,''),('cooling_power_w','Qcold','.1f',' W'),
               ('indicated_mechanical_input_power_w','Pin','.1f',' W')]
              if m.get('cooling_cop') is not None else
              [('indicated_thermal_efficiency','eta',precision,''),('indicated_power_w','Pgas','.1f',' W')])
    if m.get('cooling_power_per_total_microtube_w') is not None:
        fields.insert(0, ('cooling_power_per_total_microtube_w', 'Qcold/microtube', precision, ' W/microtube'))
    return [f'{label} {m[key]:{fmt}}{unit}' for key,label,fmt,unit in fields
            if isinstance(m.get(key),(int,float)) and math.isfinite(m[key])]


def failure_text(record, scientific):
    from .cockpit import reason_category
    if record['status']=='converged_infeasible':
        from .margins import enrich_constraints
        violated = [c for c in enrich_constraints(record.get('constraints') or [], scientific.get('constraints',[]))
                    if not c.get('satisfied') or not c.get('available')]
        if violated:
            names = ' + '.join(clean(c['name']) for c in violated[:3])
            if len(violated)>3: names += f' + {len(violated)-3} more'
            for c in violated:
                if c['name']=='maximum_absolute_mass_flow' and c.get('value') is not None and c.get('limit') is not None:
                    names += f"  mdot {c['value']:.3f} > {c['limit']:.3f}"
                    break
            return names
    reason = clean(record.get('reason') or '')
    if record['status']=='integration_failure':
        category = reason_category(record,scientific)
        if category != record['status']: return category
    # Preserve the actual microtube criteria instead of replacing them by a category.
    reason = re.sub(r'^\w+(?:Error|Exception):\s*', '', reason)
    return clipped(re.sub(r'\s*;\s*', ' + ', reason),120) or record['status']


def evaluation_text(record, scientific, new_best=False):
    origin = record.get('search_origin') or {}
    region = clean(origin.get('region_id') or 'global')
    kind = origin.get('kind')
    if kind in ('center','sobol','initial'): region += ':'+kind
    parts = [f"[{record['evaluation_number']+1:04d}] {region}  {clean(record['status'])}"]
    if new_best: parts.append('BEST')
    if record.get('cache_hit'): parts.append('(cached)')
    metrics = record.get('metrics') or {}
    if record['status']=='feasible':
        parts.extend(metric_parts(metrics))
    else:
        parts.append(failure_text(record,scientific))
    flow = metrics.get('maximum_absolute_mass_flow_kg_s')
    if isinstance(flow,(int,float)) and math.isfinite(flow) and not any('mdot ' in p for p in parts):
        parts.append(f'mdot {flow:.3f}')
    return '  '.join(parts)


class CLIProgress:
    def __init__(self, stream=None, *, scientific=None, width=None):
        self.stream = stream if stream is not None else sys.stdout
        self.tty = self.stream.isatty()
        self.scientific = scientific or {}
        self.width = width
        self.active = False

    def __enter__(self): return self

    def __exit__(self, *exception):
        self.close()
        return False

    def close(self):
        if self.active:
            self.stream.write('\r\033[K\n')
            self.stream.flush()
            self.active = False

    def _clear(self):
        if self.active:
            self.stream.write('\r\033[K')
            self.active = False

    def _width(self):
        if self.width is not None: return max(1,self.width-1)
        try: width = os.get_terminal_size(self.stream.fileno()).columns
        except (OSError,ValueError,AttributeError): width = shutil.get_terminal_size().columns
        return max(1,width-1)  # Keep away from the terminal's auto-wrap column.

    def _status(self, event, *, final=False):
        limit = event['maximum_candidates'] if event['maximum_candidates'] is not None else '∞'
        parts = [f"{event['attempted']}/{limit} · {duration(event['elapsed_seconds'])}/{duration(event['budget_seconds'])}",
                 f"feasible {event['feasible']}"]
        best = event.get('best')
        extras = metric_parts(best.get('metrics'),champion=True) if best else []
        if best:
            if extras: parts.append(extras.pop(0))
            else: parts.append(f"objective {best.get('objective',{}).get('value')}")
            parts.append(f"best {best['candidate_id'][:12]}")
            parts.extend(extras)
        else: parts.append('best —')
        parts.append(f"converged {event['converged']}")
        if final: return 'Finished · '+' · '.join(parts)
        width = self._width()
        text = clipped(parts[0],width)
        for part in parts[1:]:
            if cells(text+' · '+part)<=width: text += ' · '+part
            else: break
        return text

    def _draw(self,event):
        self.stream.write('\r\033[K'+self._status(event))
        self.stream.flush()
        self.active = True

    def __call__(self, event):
        kind = event['event']
        if kind=='start':
            limit = event['maximum_candidates'] if event['maximum_candidates'] is not None else 'unlimited'
            print(f"Phase {event['phase_id']} · {'Search step' if event.get('scheduled_search') else 'Sobol index'} {event['sobol_index']} · "
                  f"budget {duration(event['budget_seconds'])} · max {limit} new candidates",file=self.stream,flush=True)
        elif kind=='evaluation':
            self._clear()
            print(evaluation_text(event['evaluation'],self.scientific,event.get('new_best',False)),file=self.stream,flush=True)
        elif kind=='finish':
            self._clear()
            print(self._status(event,final=True),file=self.stream,flush=True)
            return
        if self.tty: self._draw(event)

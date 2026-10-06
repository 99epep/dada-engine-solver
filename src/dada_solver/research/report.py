"""Stored-result comparisons with optional, separately labelled report replays."""
from collections import Counter
import html
import math
import json
from pathlib import Path

from dada_solver.campaign.candidate import content_hash
from dada_solver.campaign.definition import runtime_identity
from dada_solver.campaign.report import elite_records

PAYLOAD_KEYS = ('schema_version','definition_id','normalized','physical','families','numerical_settings')


def verify_record(record):
    if content_hash({k:record[k] for k in PAYLOAD_KEYS}) != record['candidate_id']:
        raise ValueError('Candidate payload does not match its recorded identity.')


def inspect(path):
    """Read journal and completed orphans without creating a lock or repairing files."""
    path = Path(path)
    warnings = []
    if path.is_file():
        artifact = json.loads(path.read_text())
        if artifact.get('artifact_type') != 'research_evaluation_v1':
            raise ValueError(f'Expected a Research evaluation artifact: {path}')
        records = [artifact['record']]
        definition = artifact['definition']
        title = artifact['name']
    else:
        definition = json.loads((path/'definition.json').read_text())
        if definition.get('definition_kind') not in ('research_v1','research_v2','research_v3'):
            raise ValueError('This report requires a Dada-Engine Research study directory.')
        study = json.loads((path/'study.json').read_text())
        if study['study_id'] != definition['study_id'] or study['definition_id'] != definition['definition_id']:
            raise ValueError('Study manifest does not match the campaign definition.')
        title = study['name']
        from dada_solver.campaign.history import journal_path,read_journal
        records,_,torn_offset=read_journal(journal_path(path))
        if torn_offset is not None:
            warnings.append('Incomplete final journal append ignored for inspection; resume performs recovery.')
        from dada_solver.campaign.history import recovery_records
        known = {r['evaluation_number']:r for r in records}
        for record in recovery_records(path):
            verify_record(record)
            if record['evaluation_number'] in known and known[record['evaluation_number']] != record:
                raise ValueError('Recovery conflicts with a journal evaluation.')
            if record['evaluation_number'] not in known:
                records.append(record)
                known[record['evaluation_number']] = record
                warnings.append('Completed record awaiting journal recovery included read-only.')
        records.sort(key=lambda r:r['evaluation_number'])
        state_path = path/'state.json'
        if state_path.exists() and json.loads(state_path.read_text()).get('pending'):
            warnings.append('A pending candidate is recorded; resume will reconcile or retry it.')
    if content_hash({k:v for k,v in definition.items() if k != 'definition_id'}) != definition['definition_id']:
        raise ValueError('Stored definition does not match its identity.')
    if content_hash(definition['scientific']) != definition['study_id']:
        raise ValueError('Stored scientific inputs do not match their study identity.')
    for record in records:
        verify_record(record)
        if record['definition_id'] != definition['definition_id']:
            raise ValueError('A result belongs to a different study definition.')
    compatible = definition['runtime'] == runtime_identity()
    if not compatible: warnings.append('Stored runtime/source differs from the current checkout. Inspection is available; execution resume is not compatible.')
    scientific = definition['scientific']
    from .margins import enrich_constraints, limiting_evidence
    from .cockpit import unique_prefixes
    prefixes = unique_prefixes(r['candidate_id'] for r in records)
    for record in records:
        record['constraints']=enrich_constraints(record.get('constraints',[]),scientific['constraints'])
        from .cockpit import reason_category
        record['failure_category']=reason_category(record,scientific) if record['status']!='feasible' else None
        record['limiting_evidence']=limiting_evidence(record,scientific)
        record['resolved_parameters']=dict(scientific.get('fixed_parameters',{}),**record['physical'])
        _add_microtube_productivity(record)
        derived = record.setdefault('derived', {})
        volumes = [derived.get(f'{side}_gas_volume_m3') for side in ('heat_in', 'heat_out')]
        derived['total_exchanger_gas_volume_m3'] = (
            sum(volumes) if all(isinstance(v, (int, float)) for v in volumes) else None)
        record['parameter_units']={p['name']:p['unit'] for p in scientific['parameters']}
        record['kinematic_families']=record.get('families',{})
        record['report_source'] = str(path)
        record['display_candidate_id'] = prefixes[record['candidate_id']]
        diagnostic = record.get('diagnostics') or {}
        topology = diagnostic.get('topology')
        record['topology_display'] = topology
        if topology and not diagnostic.get('valve_events') and topology.get('classification') != 'unavailable':
            # Presentation correction only: retain the original stored diagnostic.
            record['topology_display'] = dict(classification='unavailable',
                reasons=['legacy_record_has_no_usable_valve_event_sequence'],
                stored_classification=topology.get('classification'))
    from .cockpit import campaign_evidence
    from .local_search import region_summary
    return dict(local_search=region_summary(records,scientific),cockpit=campaign_evidence(records, scientific, str(path), path.is_dir()),
                schema_version=1, name=title, study_id=definition['study_id'], definition_id=definition['definition_id'],
                scientific=scientific, runtime_compatible=compatible, warnings=warnings,
                status_counts=dict(Counter(r['status'] for r in records)), records=records,
                best=elite_records(records,5), source=str(path), replay='not_requested')


def _add_microtube_productivity(record):
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


def select_records(data, selectors=()):
    distinct = {}
    for row in data['records']: distinct[row['candidate_id']] = row
    values = list(distinct.values())
    if not selectors: return values
    selected = []
    for selector in selectors:
        if selector in ('best','second'):
            rank = 1 if selector=='best' else 2
            ranked = elite_records(values,rank)
            if len(ranked)<rank:
                raise ValueError(f'Candidate selector {selector!r} requires at least {rank} ranked feasible candidate(s); found {len(ranked)}.')
            matches = [ranked[rank-1]]
        else:
            matches = [r for r in values if r['candidate_id'].startswith(selector)]
        if len(matches) != 1:
            raise ValueError(f'Candidate selector {selector!r} is absent or ambiguous.')
        selected.append(matches[0])
    return selected


def compare(paths, selectors=(), *, plots=None, cache_directory=None, notify=lambda message: None):
    studies = [inspect(path) for path in paths]
    records = []
    from .plot_data import plot_names, candidate_plots
    requested = plot_names(plots)
    chosen = None
    if selectors:
        chosen = {r['candidate_id'] for r in select_records(
            {'records':[r for s in studies for r in s['records']]}, selectors)}
    for study in studies:
        selected = [r for r in select_records(study) if chosen is None or r['candidate_id'] in chosen]
        if requested:
            from .snapshots import stored_study
            # Default plot scope is the existing objective's two best candidates.
            # Explicit selectors request those exact candidates, without a rerank.
            targets = selected if selectors else elite_records(selected,2)
            if not targets and len(selected)==1 and selected[0]["status"]=="feasible":
                targets = selected  # A standalone legacy result needs no ranking.
            if targets:
                with stored_study(study) as snapshot:
                    for row in targets:
                        row['plots'] = candidate_plots(snapshot,row,requested,
                            cache_directory=cache_directory,notify=notify)
                        row['plots_runtime_compatible'] = study['runtime_compatible']
        study['replay'] = 'derived_report_curves' if requested else 'not_requested'
        records.extend(selected)
    if len(studies) == 1:
        result = studies[0]
        result['selected'] = records
        result['explicit_selection'] = bool(selectors)
        return result
    compatible = len({s['study_id'] for s in studies}) == 1
    result = dict(studies[0], name='Candidate comparison', records=records, selected=records,
                  status_counts=dict(Counter(r['status'] for r in records)), comparison_compatible=compatible, explicit_selection=bool(selectors),
                  cockpits=[s['cockpit'] for s in studies],
                  compared_sources=[dict(source=s['source'], study_id=s['study_id'], scientific=s['scientific']) for s in studies],
                  warnings=[w for s in studies for w in s['warnings']], best=elite_records(records,5) if compatible else [])
    from .local_search import region_summary
    result['local_search']=region_summary(records,result['scientific']) if compatible else None
    if not compatible:
        result['warnings'].append('Different scientific study identities: values are shown side by side; no combined feasible ranking is assigned. Review each source definition before interpreting differences.')
    return result


def text_report(data, *, list_candidates=False):
    lines = [data['name'], f"Study: {data['study_id']}",
             f"Attempts: {len(data['records'])}; statuses: {data['status_counts']}",
             'Power is indicated gas power. Useful output is unavailable; mechanical losses are unknown.']
    for row in data.get('selected', data['records']) if list_candidates else data['best'][:1]:
        m = row.get('metrics',{})
        if m.get('cooling_cop') is not None:
            performance=f"cooling={m.get('cooling_power_w')} W; COP={m.get('cooling_cop')}; indicated input={m.get('indicated_mechanical_input_power_w')} W"
        else:
            performance=f"power={m.get('indicated_power_w')} W; efficiency={m.get('indicated_thermal_efficiency')}"
        if m.get('cooling_power_per_total_microtube_w') is not None:
            performance += f"; Qcold/microtube={m['cooling_power_per_total_microtube_w']} W/microtube"
        if m.get('cooling_cop_times_power_per_total_microtube_w') is not None:
            performance += (f"; COP×Qcold/microtube="
                f"{m['cooling_cop_times_power_per_total_microtube_w']} W/microtube")
        volume = row.get('derived', {}).get('total_exchanger_gas_volume_m3')
        if volume is not None:
            performance += f'; exchanger gas volume={volume:g} m³'
        lines.append(f"{row['candidate_id']} {row['status']}: {performance}")
    local=data.get('local_search')
    if local:
        lines.append(f"Local regions: {len(local['regions'])}; radius={local['radius_fraction']}; round-robin")
        for region in local['regions']:
            best=region['best']; metrics=best.get('metrics',{}) if best else {}
            lines.append(f"{region['id']}: {region['attempts']} attempts, {region['converged']} converged, {region['feasible']} feasible; best objective={best['objective']['value'] if best else None}; best COP={metrics.get('cooling_cop')}; best ID={best['candidate_id'] if best else None}; rejections={region['rejection_counts']}")
    lines.extend(data['warnings'])
    return '\n'.join(lines)+'\n'


def render_html(data, destination):
    destination = Path(destination)
    if destination.suffix.lower() != '.html': raise ValueError('HTML report destination must end in .html.')
    attempts=len(data['records'])
    limit=math.ceil(attempts*.10)
    chosen=data.get('selected',data['records'])
    if data.get('comparison_compatible',True):
        chosen=elite_records(chosen,limit)
    else:
        # Incompatible studies cannot share an objective ranking.
        chosen=list({r['candidate_id']:r for r in chosen if r['status']=='feasible'}.values())[:limit]
    # Explicit comparison/plot requests remain visible even outside the detail cap.
    included={r['candidate_id'] for r in chosen}
    for record in data.get('selected',[]):
        if (record.get('plots') or data.get('explicit_selection')) and record['candidate_id'] not in included:
            chosen.append(record); included.add(record['candidate_id'])
    ranked=chosen[:2] if data.get('comparison_compatible',True) else []
    def brief(record):
        return {k:record.get(k) for k in ('candidate_id','objective','metrics')} if record else None
    from .schema_v2 import COOLING_OBJECTIVES
    data=dict(data,cooling_objective=data['scientific']['objective']['type'] in COOLING_OBJECTIVES,
        selected=chosen,comparison_default_ids=[r['candidate_id'] for r in ranked],
        html_selection=dict(policy='best_distinct_feasible_tenth_v1',attempts=attempts,
                            limit=limit,retained=len(chosen)),
        best=[brief(r) for r in data.get('best',[])],
        records=[dict(candidate_id=r['candidate_id'],status=r['status'],
                      duration_seconds=r.get('duration_seconds'),metrics={k:r.get('metrics',{}).get(k)
                      for k in ('cooling_cop','indicated_thermal_efficiency',
                                'cooling_power_per_total_microtube_w',
                                'cooling_cop_times_power_per_total_microtube_w')}) for r in data['records']])
    if data.get('local_search'):
        data['local_search']=dict(data['local_search'],regions=[dict(r,best=brief(r.get('best')))
            for r in data['local_search']['regions']])
    escaped_data = json.dumps(data, allow_nan=False).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    page = HTML.replace('REPORT_CURVES_SCRIPT', Path(__file__).with_name('report_curves.js').read_text()).replace('TITLE_TEXT', html.escape(data['name'])).replace('EMBEDDED_DATA', escaped_data)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(page)
    return destination


HTML = '''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>TITLE_TEXT</title>
<style>
#candidateViewport{max-height:640px;overflow:auto;padding:0}#candidates td{height:24px;white-space:nowrap}#candidates th{position:sticky;top:0;z-index:1}#candidates tbody tr{cursor:pointer}#candidates tbody tr:hover{background:#edf6fa}#bounds{max-height:420px}#bounds th{position:sticky;top:0}#suggestions{max-height:440px}body{font:15px/1.5 system-ui,sans-serif;color:#182938;background:#f4f6f8;margin:0}main{max-width:1200px;margin:auto;padding:28px}h1{margin-bottom:4px}h2{font-size:21px;margin-top:30px}.card{background:white;border:1px solid #d9e0e6;border-radius:8px;padding:18px;margin:16px 0;overflow:auto}.muted{color:#506475}.warning{border-left:4px solid #a65b00;background:#fff4df;padding:12px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;padding:8px;border-bottom:1px solid #d9e0e6;vertical-align:top}th{background:#edf2f5}code{overflow-wrap:anywhere}th button{font:inherit;font-weight:600;border:0;background:transparent;cursor:pointer;text-align:left;color:inherit}select,input{font:inherit;padding:7px;max-width:100%;margin:5px}label{display:inline-block;margin-right:15px}pre{white-space:pre-wrap;word-break:break-word;font-size:12px}.plots{display:grid;grid-template-columns:1fr 1fr;gap:16px}svg{width:100%;height:auto}.near{background:#fff1cb;color:#775100}.good{color:#087660}.bad{color:#a33b35}a{color:#17658b}@media(max-width:750px){.plots{grid-template-columns:1fr}main{padding:12px}}
</style><main><section class="card" id="localRegions" hidden><h2>Local regions</h2><div id="regionSummary"></div></section>
<h1>TITLE_TEXT</h1><p class="muted">Dada-Engine Research · offline report · stored metrics unchanged · reconstructed curves labelled separately</p>
<div id="warnings"></div><div class="card" id="summary"></div>
<h2>Campaign funnel</h2><div class="card" id="funnel"></div>
<h2>Pressure on parameter bounds</h2><div class="card" id="bounds"></div>
<h2>Suggested next steps</h2><p class="muted">Deterministic observations from stored data. Commands are suggestions for human review; this page never executes them. Run commands from the checkout used to generate this report. Resume still checks runtime compatibility.</p>
<div class="card" id="suggestions"></div><div class="card" id="selectedCommands"></div>
<p id="boundary">Efficiency uses indicated gas work divided by external-stream heat input. Useful shaft power is unavailable; mechanical losses are unknown. External-loop hydraulic losses and pump/fan consumption are excluded from the balance, not physically zero.</p>
<div class="plots"><div class="card"><h2 id="performanceTitle">Power and efficiency</h2><div id="scatter"></div></div><div class="card"><h2>Evaluation progress</h2><div id="progress"></div></div></div>
<h2>Candidates</h2><label>Status <select id="status"><option value="all">All statuses</option></select></label><label><input type="checkbox" id="validOnly">All constraints available and satisfied</label>
<label>Sort by <select id="sortBy"></select></label><label>Direction <select id="sortDirection"><option value="asc">Ascending</option><option value="desc">Descending</option></select></label>
<p class="muted">Click a column heading to sort. Missing values stay last. Additional metrics, active parameters and constraint margins are available in Sort by.</p>
<label>Failure category <select id="failureCategory"><option value="all">All categories</option></select></label>
<label>Violated constraint <select id="violatedConstraint"><option value="all">All constraints</option></select></label>
<button id="clearInspection" type="button">Clear inspection filters</button><p id="inspectionContext" class="muted"></p>
<p class="muted">All records remain available. Scroll for additional candidates; click a row to inspect it as Candidate A.</p>
<div class="card" id="candidateViewport" tabindex="0" aria-label="Scrollable candidate table"><table id="candidates"></table></div>
<div id="additionalPlots"></div><div id="mechanismPlots"></div>
<section id="volumeSection" hidden><h2>Cylinder volumes</h2><p>One cycle in solver angle, with the operation convention applied once. Solid: small cylinder; dashed: large cylinder. No thermodynamic integration.</p><div class="card" id="volumes"></div></section>
<h2>Compare selected candidates</h2><p>Values are absolute; the difference is B − A. Missing data stays unavailable.</p>
<label>Candidate A <select id="left"></select></label><label>Candidate B <select id="right"></select></label>
<div class="card"><table id="comparison"></table></div>
<h2>Constraint values and margins</h2><p>Stored constraints and model-domain evidence are ordered by violation, unavailable evidence, then proximity to a limit. Positive margins satisfy the displayed boundary. Near boundary means within 5% of a nonzero limit: a display cue, not an optimization active-set or safety factor. Categorical verdicts have no relative margin. Trial-state rejections and sampled cycle extrema are labelled separately; missing boundaries remain unavailable.</p>
<div class="card"><h3 id="constraintContextA">Candidate A</h3><table id="constraintsA"></table><h3 id="constraintContextB">Candidate B</h3><table id="constraintsB"></table></div>
<section id="trialHistory" hidden><h2>Rejected trial history</h2><p>Historical trial states, separate from final-cycle constraints. A first trial snapshot is not a measurement of the final periodic cycle, nor necessarily the last failed retry.</p><div class="card" id="trialDiagnostics"></div></section>
<details class="card" id="candidateDetails"><summary>Selected candidate details, constraints and convergence</summary><pre id="detail"></pre></details>
<details class="card"><summary>Scientific definition and provenance</summary><pre id="provenance"></pre></details>
<p class="muted">Closure screens are sampled numerical checks, not continuous feasibility proofs. This validation study does not establish a global optimum. Volume curves are available on request. Thermodynamic trajectories and animations are not generated by this report.</p>
</main><script id="data" type="application/json">EMBEDDED_DATA</script><script>
'use strict';
const d=JSON.parse(document.getElementById('data').textContent),byId=id=>document.getElementById(id);
const fmt=v=>v===null||v===undefined?'unavailable':typeof v==='number'?Number(v.toPrecision(8)).toString():String(v);
const esc=v=>String(v).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const rows=[...new Map((d.selected||d.records).map(r=>[r.candidate_id,r])).values()];
byId('warnings').innerHTML=d.warnings.map(w=>'<p class="warning">'+esc(w)+'</p>').join('');
const fixed=d.scientific.fixed||d.scientific.fixed_parameters||{},cooling=d.cooling_objective;
const objective=d.scientific.objective.type;
const productivity=objective==='maximize_cooling_power_per_total_microtube';
const composite=objective==='maximize_cooling_cop_times_power_per_total_microtube';
const ykey=composite?'cooling_cop_times_power_per_total_microtube_w':productivity?'cooling_power_per_total_microtube_w':cooling?'cooling_cop':'indicated_thermal_efficiency',yscale=cooling?1:100,ylabel=composite?'COP×Qcold/microtube [W/microtube]':productivity?'Cooling power per microtube [W/microtube]':cooling?'Cooling COP [1]':'Indicated efficiency [%]';
byId('summary').innerHTML='<b>Study</b> <code>'+esc(d.study_id)+'</code><p>'+d.records.length+' attempts · '+esc(JSON.stringify(d.status_counts))+'</p><p>Families: '+esc(d.scientific.kinematics?d.scientific.kinematics.small.family+' / '+d.scientific.kinematics.large.family:JSON.stringify(d.scientific.families))+'</p>';
byId('summary').innerHTML+='<p>Detailed candidates: '+d.html_selection.retained+' (base limit '+d.html_selection.limit+': best distinct feasible 10%, rounded up; plot targets and explicit selections also retained). Global statistics and progress include all attempts. Complete records remain in the journal and JSON report.</p>';
if(cooling){byId('performanceTitle').textContent=composite?'COP×Qcold/microtube':productivity?'Cooling power per microtube':'Indicated power and cooling COP';byId('boundary').textContent='Cooling power and COP use heat absorbed at the cold external-stream boundary and indicated mechanical input. External-loop hydraulics and pump/fan consumption are unmodelled. Shaft losses and useful human input remain uncalibrated; no mechanical efficiency is assumed.';}
byId('provenance').textContent=JSON.stringify(d.compared_sources||d.scientific,null,2);
for(const status of Object.keys(d.status_counts)){const o=new Option(status,status);byId('status').add(o);}
for(const [i,r] of rows.entries())for(const id of ['left','right'])byId(id).add(new Option(r.candidate_id.slice(0,12)+' · '+r.status,String(i)));
const defaults=(d.comparison_default_ids||[]).map(id=>rows.findIndex(r=>r.candidate_id===id)).filter(i=>i>=0);
for(let i=0;defaults.length<2&&i<rows.length;i++)if(!defaults.includes(i))defaults.push(i);
byId('left').selectedIndex=defaults[0]??-1;byId('right').selectedIndex=defaults[1]??defaults[0]??-1;
function table(id,head,body){byId(id).innerHTML='<thead><tr>'+head.map(h=>'<th>'+esc(h)+'</th>').join('')+'</tr></thead><tbody>'+body.map(r=>'<tr>'+r.map(v=>'<td>'+esc(fmt(v))+'</td>').join('')+'</tr>').join('')+'</tbody>';}
function commandBox(signal,command,parent){
 const box=document.createElement('div'),p=document.createElement('p'),code=document.createElement('code'),button=document.createElement('button');
 p.textContent=signal;code.textContent=command;button.textContent='Copy';button.type='button';button.style.margin='8px';
 button.onclick=async()=>{try{if(navigator.clipboard?.writeText)await navigator.clipboard.writeText(command);else throw Error('clipboard unavailable');button.textContent='Copied';}
 catch(_){const input=document.createElement('textarea');input.value=command;document.body.append(input);input.select();const ok=document.execCommand('copy');input.remove();button.textContent=ok?'Copied':'Select command to copy';}};
 box.append(p,code,button);parent.append(box);
}
function cockpit(){
 for(const evidence of d.cockpits||[d.cockpit]){
  if(!evidence)continue;
  const f=evidence.funnel,section=document.createElement('section');
  const title=document.createElement('b');title.textContent=evidence.source;section.append(title);
  const counts=document.createElement('p');counts.textContent='Attempted '+f.attempted+' · Integrated '+f.integrated+' · Converged '+f.converged+' · Feasible '+f.feasible+' · Cache hits '+f.cache_hits+' · Distinct '+f.distinct_candidates;section.append(counts);
  const rejects=document.createElement('p');rejects.textContent='Rejections: '+(Object.entries(evidence.rejection_categories).map(([k,v])=>k+'='+v).join(' · ')||'none');section.append(rejects);byId('funnel').append(section);
  const note=document.createElement('p');note.textContent=evidence.source+' — '+evidence.bound_convention;byId('bounds').append(note);
  const tableElement=document.createElement('table');tableElement.innerHTML='<thead><tr><th>Parameter</th><th>Bound</th><th>Limit</th><th>Global near / total</th><th>Global %</th><th>Elites near / total</th></tr></thead><tbody>'+evidence.bounds.map(b=>'<tr class="'+(b.pressed?'near':'')+'">'+[b.parameter,b.side,fmt(b.limit)+' '+b.unit,b.count+' / '+b.total,b.percent.toFixed(1)+'%',b.elite_count+' / '+b.elite_total].map(x=>'<td>'+esc(x)+'</td>').join('')+'</tr>').join('')+'</tbody>';byId('bounds').append(tableElement);
  for(const suggestion of evidence.suggestions){
   if(suggestion.command){commandBox(suggestion.signal,suggestion.command,byId('suggestions'));continue;}
   const box=document.createElement('div'),p=document.createElement('p');p.textContent=suggestion.signal;box.append(p);
   if(suggestion.action){const button=document.createElement('button');button.type='button';button.textContent=suggestion.action.type==='bounds'?'Inspect bound pressure':'Inspect matching candidates';button.onclick=()=>inspectAction(suggestion.action);box.append(button);}
   byId('suggestions').append(box);
  }
 }
}
const shellQuote=value=>/^[a-zA-Z0-9_./:=@%-]+$/.test(value)?value:String.fromCharCode(39)+value.split(String.fromCharCode(39)).join(String.fromCharCode(39,34,39,34,39))+String.fromCharCode(39);
function selectedCommands(a,b){
 const parent=byId('selectedCommands');parent.replaceChildren();
 const selected=[...new Map([a,b].map(r=>[r.candidate_id,r])).values()],groups=new Map();
 for(const r of selected){const source=r.report_source||d.source;if(!groups.has(source))groups.set(source,[]);groups.get(source).push(r);}
 for(const [source,items] of groups){
  const output=source.endsWith('.json')?source.slice(0,-5)+'.selected-volumes.html':source+'/selected-volumes.html';
  const args=['dada-research','report',source];for(const r of items)args.push('--candidate',r.display_candidate_id||r.candidate_id);
  args.push('--plots','volumes','--html',output);
  commandBox('Compare selected candidate volumes without thermodynamic integration.',args.map(shellQuote).join(' '),parent);
 }
}
function scatter(id,points,xlabel,ylabel){
 const w=500,h=280,l=68,b=48,t=15,r=20;
 if(!points.length){byId(id).textContent='No available converged metrics.';return;}
 let xs=points.map(p=>p[0]),ys=points.map(p=>p[1]),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
 let dx=(xmax-xmin)||Math.max(1,Math.abs(xmax)*.05),dy=(ymax-ymin)||1;xmin-=dx*.07;xmax+=dx*.07;ymin-=dy*.07;ymax+=dy*.07;
 const x=v=>l+(v-xmin)/(xmax-xmin)*(w-l-r),y=v=>h-b-(v-ymin)/(ymax-ymin)*(h-b-t);
 const tick=(v,span)=>v.toFixed(Math.max(0,Math.min(9,1-Math.floor(Math.log10(span/4)))));
 let s='<svg role="img" aria-label="'+esc(xlabel+' versus '+ylabel)+'" viewBox="0 0 '+w+' '+h+'"><path d="M'+l+' '+t+' V'+(h-b)+' H'+(w-r)+'" fill="none" stroke="#526575"/>';
 for(let i=0;i<=4;i++){let xv=xmin+(xmax-xmin)*i/4,yv=ymin+(ymax-ymin)*i/4;s+='<text x="'+x(xv)+'" y="'+(h-b+18)+'" font-size="10" text-anchor="middle">'+tick(xv,xmax-xmin)+'</text><text x="'+(l-6)+'" y="'+y(yv)+'" font-size="10" text-anchor="end">'+tick(yv,ymax-ymin)+'</text>';}
 for(const p of points)s+='<circle cx="'+x(p[0])+'" cy="'+y(p[1])+'" r="4" fill="'+(p[3]==='feasible'?'#087660':'#b34d39')+'"><title>'+esc(p[2])+'</title></circle>';
 s+='<text x="280" y="273" text-anchor="middle" font-size="12">'+esc(xlabel)+'</text><text transform="translate(14 130) rotate(-90)" text-anchor="middle" font-size="12">'+esc(ylabel)+'</text></svg>';byId(id).innerHTML=s;
}
const refluxMin=r=>{const v=Object.values(r.derived?.local_reflux?.minimum_signed_flows_kg_s||{}).filter(Number.isFinite);return v.length?Math.min(...v):null;};
if(d.local_search){
 byId('localRegions').hidden=false;
 const root=byId('regionSummary');
 const intro=document.createElement('p');intro.textContent=d.local_search.regions.length+' regions · normalized radius '+d.local_search.radius_fraction+' · round-robin';root.append(intro);
 for(const r of d.local_search.regions){
  const p=document.createElement('p');p.textContent=r.id+': '+r.attempts+' attempts, '+r.converged+' converged, '+r.feasible+' feasible; best objective '+(r.best?.objective?.value??'unavailable')+'; best COP '+(r.best?.metrics?.cooling_cop??'unavailable')+'; best '+(r.best?.candidate_id??'none');root.append(p);
  const detail=document.createElement('details');const title=document.createElement('summary');title.textContent='Center and rejection evidence';detail.append(title);
  const pre=document.createElement('pre');pre.textContent=JSON.stringify({source_candidate_id:r.source_candidate_id,center:r.center,rejection_rates:r.rejection_rates},null,2);detail.append(pre);root.append(detail);
 }
}
const columns=[
 ['basin','Basin',r=>r.search_origin?.region_id||'global'],
 ['origin','Search origin',r=>r.search_origin?.kind||'global Sobol'],
 ['id','Candidate ID',r=>r.candidate_id],['status','Status',r=>r.status],
 ['input','Indicated input [W]',r=>r.metrics.indicated_mechanical_input_power_w],
 ['power','Indicated gas power [W]',r=>r.metrics.indicated_power_w],
 ['mass','Gas inventory [kg]',r=>r.metrics.total_mass_kg],
 ['cooling','Cooling power [W]',r=>r.metrics.cooling_power_w],
 ['microtube','Cooling power per microtube [W/microtube]',r=>r.metrics.cooling_power_per_total_microtube_w],
 ['composite','COP × cooling power per microtube [W/microtube]',r=>r.metrics.cooling_cop_times_power_per_total_microtube_w],
 ['hxvolume','Total exchanger gas volume [m³]',r=>r.derived?.total_exchanger_gas_volume_m3],
 ['cop','Cooling COP [1]',r=>r.metrics.cooling_cop],
 ['efficiency','Efficiency [1]',r=>r.metrics.indicated_thermal_efficiency],
 ['topology','Topology',r=>r.topology_display?.classification],
 ['reflux','Reflux detected',r=>r.derived?.local_reflux?.detected],
 ['minimumFlow','Worst signed flow [kg/s]',refluxMin],
 ['cycles','Cycles',r=>r.periodic_cycle_count],['duration','Duration [s]',r=>r.duration_seconds]];
const sortColumns=[...columns,
 ['pressure','Pressure max [Pa]',r=>r.metrics.maximum_pressure_pa],
 ['temperature','Temperature max [K]',r=>r.metrics.maximum_temperature_k],
 ['flow','Mass flow max [kg/s]',r=>r.metrics.maximum_absolute_mass_flow_kg_s]];
for(const key of new Set(rows.flatMap(r=>Object.keys(r.physical||{}))))sortColumns.push(['parameter:'+key,'Parameter: '+key,r=>r.physical?.[key]]);
for(const key of new Set(rows.flatMap(r=>(r.constraints||[]).map(c=>c.name)))){
 for(const field of ['margin','relative_margin'])sortColumns.push(['constraint:'+key+':'+field,key+' '+field,r=>{const c=r.constraints?.find(c=>c.name===key);return c?.available?c[field]:null;}]);
}
for(const [key,label] of sortColumns)byId('sortBy').add(new Option(label,key));
byId('sortBy').value=composite?'composite':productivity?'microtube':cooling?'cop':'efficiency';byId('sortDirection').value='desc';
function sortedRows(visible){
 const get=sortColumns.find(c=>c[0]===byId('sortBy').value)[2],sign=byId('sortDirection').value==='asc'?1:-1;
 const missing=v=>v===null||v===undefined||(typeof v==='number'&&!Number.isFinite(v));
 return [...visible].sort((a,b)=>{const x=get(a),y=get(b);if(missing(x)||missing(y))return missing(x)===missing(y)?a.candidate_id.localeCompare(b.candidate_id):missing(x)?1:-1;
 const order=typeof x==='number'||typeof x==='boolean'?Number(x)-Number(y):String(x).localeCompare(String(y));return sign*order||a.candidate_id.localeCompare(b.candidate_id);});
}
let inspectionSource=null;
for(const category of [...new Set(rows.map(r=>r.failure_category).filter(Boolean))].sort())byId('failureCategory').add(new Option(category,category));
for(const name of [...new Set(rows.flatMap(r=>(r.constraints||[]).filter(c=>c.available&&!c.satisfied).map(c=>c.name)))].sort())byId('violatedConstraint').add(new Option(name,name));
function matchesInspection(r){return (!inspectionSource||r.report_source===inspectionSource)&&
 (byId('failureCategory').value==='all'||r.failure_category===byId('failureCategory').value)&&
 (byId('violatedConstraint').value==='all'||r.constraints.some(c=>c.name===byId('violatedConstraint').value&&c.available&&!c.satisfied));}
function clearInspection(){inspectionSource=null;byId('failureCategory').value='all';byId('violatedConstraint').value='all';byId('status').value='all';byId('validOnly').checked=false;}
function inspectAction(action){
 if(action.type==='bounds'){byId('bounds').scrollIntoView({block:'center'});return;}
 clearInspection();inspectionSource=action.source;
 const select=byId(action.field==='failure_category'?'failureCategory':'violatedConstraint');
 if(![...select.options].some(option=>option.value===action.value))select.add(new Option(action.value,action.value));
 select.value=action.value;
 show();byId('candidateViewport').scrollIntoView({block:'center'});
 const first=sortedRows(rows.filter(matchesInspection))[0];
 if(first){byId('left').value=String(rows.indexOf(first));compare();byId('candidateDetails').open=true;}
}
function show(){
 const visible=sortedRows(rows.filter(r=>matchesInspection(r)&&(byId('status').value==='all'||r.status===byId('status').value)&&(!byId('validOnly').checked||r.constraints.length>0&&r.constraints.every(c=>c.available&&c.satisfied))));
 const shown=columns.filter(c=>cooling?c[0]!=='efficiency':!['input','cooling','cop'].includes(c[0]));
 table('candidates',shown.map(c=>c[1]),visible.map(r=>shown.map(c=>c[2](r))));
 byId('inspectionContext').textContent=visible.length+' / '+rows.length+' distinct candidates shown'+(inspectionSource?' · Source: '+inspectionSource:'')+(!visible.length?' · No matching records in this embedded selection. Clear filters or inspect the complete JSON report for omitted candidates.':'');
 const idColumn=shown.findIndex(c=>c[0]==='id');
 [...byId('candidates').querySelectorAll('tbody tr')].forEach((tr,i)=>{tr.cells[idColumn].title=visible[i].candidate_id;tr.cells[idColumn].textContent=visible[i].candidate_id.slice(0,12);
 tr.onclick=()=>{byId('left').value=String(rows.indexOf(visible[i]));compare();byId('candidateDetails').open=true;};});
 [...byId('candidates').querySelectorAll('th')].forEach((th,i)=>{const key=shown[i][0],selected=byId('sortBy').value===key;
 th.setAttribute('aria-sort',selected?(byId('sortDirection').value==='asc'?'ascending':'descending'):'none');
 const button=document.createElement('button');button.textContent=shown[i][1]+(selected?(byId('sortDirection').value==='asc'?' ↑':' ↓'):'');
 button.onclick=()=>{byId('sortDirection').value=selected&&byId('sortDirection').value==='asc'?'desc':'asc';byId('sortBy').value=key;show();};th.replaceChildren(button);
 });
 scatter('scatter',visible.filter(r=>Number.isFinite(r.metrics.indicated_power_w)&&Number.isFinite(r.metrics[ykey])).map(r=>[r.metrics.indicated_power_w,yscale*r.metrics[ykey],r.candidate_id,r.status]),'Indicated gas power [W]',ylabel);
 let elapsed=0,points=[];for(const r of d.records){elapsed+=r.duration_seconds||0;if(Number.isFinite(r.metrics[ykey]))points.push([elapsed,yscale*r.metrics[ykey],r.candidate_id,r.status]);}scatter('progress',points,'Cumulative evaluation time [s]',ylabel);
}
function volumePlots(){
 const plotted=rows.filter(r=>r.plots?.volumes);if(!plotted.length)return;
 byId('volumeSection').hidden=false;
 const w=1000,h=400,l=80,b=55,t=20,right=20,colors=['#17658b','#b74424','#087660','#773baa','#986700'];
 const axis=niceAxis(0,Math.max(...plotted.flatMap(r=>[...r.plots.volumes.small,...r.plots.volumes.large]))*1000,true),ymax=axis.high/1000;
 const x=v=>l+v/360*(w-l-right),y=v=>h-b-v/ymax*(h-b-t);
 let svg='<svg role="img" aria-label="Small and large cylinder volumes over one cycle" viewBox="0 0 '+w+' '+h+'"><path d="M'+l+' '+t+' V'+(h-b)+' H'+(w-right)+'" fill="none" stroke="#526575"/>';
 for(let i=0;i<=6;i++)svg+='<text x="'+x(i*60)+'" y="'+(h-b+22)+'" text-anchor="middle" font-size="13">'+i*60+'</text>';
 for(const v of axis.ticks)svg+='<text x="'+(l-8)+'" y="'+(y(v/1000)+4)+'" text-anchor="end" font-size="13">'+tickLabel(v,axis.step)+'</text>';
 for(const [i,row] of plotted.entries())for(const side of ['small','large']){
  const p=row.plots.volumes,color=colors[i%colors.length],points=p.angle.map((a,j)=>x(a)+','+y(p[side][j])).join(' ');
  svg+='<polyline data-candidate="'+esc(row.candidate_id)+'" data-side="'+side+'" points="'+points+'" fill="none" stroke="'+color+'" stroke-width="2"'+(side==='large'?' stroke-dasharray="7 4"':'')+'><title>'+esc(row.candidate_id+' · '+side)+'</title></polyline>';
 }
 svg+='<text x="530" y="395" text-anchor="middle">Solver cycle angle [deg]</text><text transform="translate(20 190) rotate(-90)" text-anchor="middle">Cylinder volume [L]</text></svg>';
 byId('volumes').innerHTML=svg+'<ul>'+plotted.map((r,i)=>'<li style="color:'+colors[i%colors.length]+'"><code>'+esc(r.candidate_id)+'</code> · solid small / dashed large</li>').join('')+'</ul>';
}
function compare(){
 if(!rows.length){byId('comparison').textContent='No feasible candidates retained in this HTML report.';return;}
 const a=rows[Number(byId('left').value)],b=rows[Number(byId('right').value)],out=[];
 selectedCommands(a,b);
 const add=(label,x,y)=>out.push([label,x,y,typeof x==='number'&&typeof y==='number'?y-x:null]);
 const ap=a.resolved_parameters||a.physical,bp=b.resolved_parameters||b.physical;
 const active=new Set([...Object.keys(a.physical||{}),...Object.keys(b.physical||{})]);
 const parameter=name=>add(name+' ['+(a.parameter_units?.[name]||b.parameter_units?.[name]||'see definition')+']',ap[name],bp[name]);
 for(const name of active)parameter(name);
 for(const side of ['small','large'])add(side+' kinematic family',a.families?.[side]||a.families?.kinematics,b.families?.[side]||b.families?.kinematics);
 for(const name of new Set([...Object.keys(ap),...Object.keys(bp)]))if(!active.has(name))parameter(name);
 add('Total exchanger gas volume [m³]',a.derived?.total_exchanger_gas_volume_m3,b.derived?.total_exchanger_gas_volume_m3);
 for(const [key,unit] of Object.entries({total_mass_kg:'kg',indicated_power_w:'W',indicated_thermal_efficiency:'1',heat_input_w:'W',cooling_power_w:'W',cooling_cop:'1',cooling_power_per_total_microtube_w:'W/microtube',cooling_cop_times_power_per_total_microtube_w:'W/microtube',indicated_mechanical_input_power_w:'W',maximum_pressure_pa:'Pa',maximum_temperature_k:'K',maximum_absolute_mass_flow_kg_s:'kg/s',useful_mechanical_power_w:'W'}))add(key+' ['+unit+']',a.metrics[key],b.metrics[key]);
 for(const side of ['heat_in','heat_out'])for(const key of ['geometry_model','bundle_diameter_m','bundle_face_area_m2','pitch_m','tube_flow_area_m2','conduit_area_m2','conduit_diameter_m','conduit_area_ratio','collector_half_angle_deg','collector_height_m','header_gas_volume_m3','tube_gas_volume_m3','additional_internal_volume_m3','working_gas_volume_m3','valve_model','valve_cda_m2','header_loss_coefficient','header_loss_model']){
  const x=a.derived?.hardware?.[side]?.[key],y=b.derived?.hardware?.[side]?.[key];if(x!=null||y!=null)add(side+' '+key,x,y);
 }
 for(const side of ['heat_in','heat_out'])add(side+' external / peak internal capacity rate [1]',a.derived?.external_air_capacity_diagnostics?.[side]?.external_to_peak_internal_capacity_rate_ratio,b.derived?.external_air_capacity_diagnostics?.[side]?.external_to_peak_internal_capacity_rate_ratio);
 for(const side of ['heat_in','heat_out'])for(const key of ['fluid','inlet_temperature_k','outlet_minimum_k','outlet_maximum_k','mass_flow_kg_s','cp_j_kg_k','capacity_rate_w_k','wall_conductance_w_k','heat_into_machine_per_cycle_j','mean_heat_into_machine_w','external_loop_losses'])add(side+' external '+key,a.derived?.external_streams?.[side]?.[key],b.derived?.external_streams?.[side]?.[key]);
 for(const port of ['small_to_cold','cold_to_large','large_to_hot','hot_to_small'])add(port+' minimum signed flow [kg/s]',a.derived?.local_reflux?.minimum_signed_flows_kg_s?.[port],b.derived?.local_reflux?.minimum_signed_flows_kg_s?.[port]);
 for(const c of a.constraints){const other=b.constraints.find(x=>x.name===c.name);const unit=c.unit||d.scientific.constraints.find(x=>x.type===c.name)?.unit||'1';add(c.name+' margin ['+unit+']',c.available?c.margin:null,other?.available?other.margin:null);}
 for(const [id,candidate] of [['constraintsA',a],['constraintsB',b]]){
  const evidence=candidate.limiting_evidence||candidate.constraints;
  byId(id==='constraintsA'?'constraintContextA':'constraintContextB').textContent=(id==='constraintsA'?'Candidate A':'Candidate B')+' · '+candidate.candidate_id.slice(0,12)+' · '+candidate.status+' · '+(candidate.search_origin?.region_id||'global');
  table(id,['Constraint / boundary','Current value','Relation','Limit','Unit','Signed margin','Relative margin','Near boundary','State','Scope / context','Evidence'],evidence.map(c=>[c.name,c.value,c.relation,c.limit,c.unit,c.margin,c.relative_margin,c.near_active,c.state,[c.scope,c.context].filter(Boolean).join(' / '),c.method]));
  [...byId(id).querySelectorAll('tbody tr')].forEach((tr,i)=>{const c=evidence[i];tr.className=!c.available?'muted':!c.satisfied?'bad':c.near_active?'near':'good';if(c.near_active)tr.setAttribute('aria-label',c.name+' near active limit');});
 }
 const history=byId('trialDiagnostics');history.replaceChildren();
 for(const [label,candidate] of [['A',a],['B',b]]){
  const failure=candidate.diagnostics?.first_microtube_failure;if(!failure)continue;
  const box=document.createElement('div'),title=document.createElement('h3'),note=document.createElement('p'),detail=document.createElement('pre');
  const recovered=candidate.converged||['feasible','converged_infeasible'].includes(candidate.status);
  title.textContent='Candidate '+label+' · '+(recovered?'Recovered trial diagnostic':'Rejected trial diagnostic');
  note.textContent='First rejected trial · safe uniform retry: '+(candidate.safe_retry_used===true?'used':candidate.safe_retry_used===false?'not used':'not recorded')+' · final result: '+candidate.status;
  detail.textContent=JSON.stringify(failure,null,2);box.append(title,note,detail);history.append(box);
 }
 byId('trialHistory').hidden=!history.childElementCount;
 table('comparison',['Quantity','A','B','B − A'],out);byId('detail').textContent=JSON.stringify({A:a,B:b},null,2);
}
byId('clearInspection').onclick=()=>{clearInspection();show();};byId('failureCategory').onchange=show;byId('violatedConstraint').onchange=show;
REPORT_CURVES_SCRIPT
cockpit();byId('sortBy').onchange=show;byId('sortDirection').onchange=show;volumePlots();byId('status').onchange=show;byId('validOnly').onchange=show;byId('left').onchange=compare;byId('right').onchange=compare;show();compare();
</script></html>'''

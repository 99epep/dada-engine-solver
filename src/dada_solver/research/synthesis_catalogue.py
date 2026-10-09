"""Standalone sortable human catalogue, drawn from production mechanisms."""
import html
import io
import json
import shlex
from pathlib import Path

import numpy as np
from .artifacts import MechanismArtifact
from .mechanism_view import MechanismModel
from .motion_target import PERIOD, PeriodicTargetSide


def _preview(artifact,target,side):
    from matplotlib.figure import Figure
    figure=Figure(figsize=(10,3))
    axes=figure.subplots(1,3)
    model=MechanismModel(artifact)
    angles=np.linspace(0.,PERIOD,361)
    reference=PeriodicTargetSide(target,side)
    states=[model.state(float(t)) for t in np.linspace(0.,PERIOD,49)]
    first=states[0]
    for name in first['joints']:
        trajectory=np.array([s['joints'][name] for s in states])
        axes[0].plot(trajectory[:,0],trajectory[:,1],':',alpha=.25)
        point=first['joints'][name]
        axes[0].text(*point,name,fontsize=7)
    for left,right in first['links']:
        points=np.array([first['joints'][left],first['joints'][right]])
        axes[0].plot(points[:,0],points[:,1],'o-',ms=3,lw=1.5)
    axes[0].set_aspect('equal'); axes[0].set_title('Production geometry / crank radius',fontsize=8)
    if model.primary:
        from .synthesis_six_bar import primary_evidence, PrimaryProjection
        axis=primary_evidence(artifact,target,side)['primary_projection']['axis']
        model.law=PrimaryProjection(model.geometry,axis=axis)
    position=model.law.value(angles)-1.
    velocity=model.law.value(angles,1)
    axes[1].plot(angles,reference.position(angles),'k--',label='Target')
    axes[1].plot(angles,position,label='E projection (not piston)' if model.primary else 'Mechanism')
    axes[2].plot(angles,velocity,label='Mechanism')
    if reference.has_velocity: axes[2].plot(angles,reference.velocity(angles),'k--',label='Target')
    for axis,label in zip(axes[1:],('Position / stroke','Velocity / rad')):
        axis.set_title(label,fontsize=8); axis.set_xlim(0.,PERIOD); axis.set_xlabel('Study angle / rad',fontsize=8)
        axis.legend(fontsize=7); axis.grid(alpha=.2)
    figure.tight_layout()
    stream=io.StringIO(); figure.savefig(stream,format='svg')
    return stream.getvalue()[stream.getvalue().index('<svg'):]


def render_synthesis_catalogue(library,target,path,*,library_path=None):
    """Write inspectable evidence and family IDs; no absolute winner is designated."""
    path=Path(path)
    if path.exists(): raise ValueError('Synthesis catalogue already exists; choose a new path.')
    library_path=Path(library_path) if library_path is not None else path.with_suffix('.json')
    common_command='dada-research mechanism pair '+shlex.quote(str(library_path))
    entries=[]
    for member in library.members:
        if member['metadata'].get('target_hash')!=target.content_hash:
            raise ValueError('Catalogue target and library identities disagree.')
        for side,raw in member['mechanisms'].items():
            artifact=MechanismArtifact.from_data(raw)
            evidence=member['metadata']['evidence'][side]
            fit,mechanical=evidence['fit'] or evidence.get('primary_projection'),evidence['mechanical']
            categories=artifact.data['provenance'].get('categories',member['metadata']['provenance'].get('categories',{}))
            labels=[member['family_id'],side,artifact.scientific['settings']['family']+(' primary (E proxy)' if artifact.scientific['settings'].get('component')=='primary' else ''),
                    fit['position_rms'],fit['position_maximum_error'],fit['velocity_rms_per_rad']]
            metrics={k:v for k,v in mechanical['metrics'].items() if k in ('stroke_over_crank',
                'minimum_primary_transmission_sine','minimum_secondary_transmission_sine','minimum_rod_axis_cosine','stroke_over_envelope','closure_margin','E_projection_span_over_crank','E_span_over_crank','E_axis_variance_fraction','E_lateral_rms_over_projection_span','EH_over_crank','H_axis_lateral_rms_over_stroke','H_axis_lateral_span_over_stroke','crank_axis_to_EFH_clearance_over_crank')}
            screen=mechanical.get('primary_mechanical_screen',mechanical.get('mechanical_screen',artifact.data['provenance'].get('mechanical_screen',{})))
            records=[*mechanical['constraints'],*mechanical.get('primary_constraints',[])]
            admissibility=dict(satisfied=all(r['satisfied'] for r in records),samples=mechanical.get('samples'),
                               mechanical_screen=screen,constraints=records,envelope_frame=mechanical.get('envelope_frame'))
            summary=html.escape('; '.join(f'{k}: {v:.5g}' for k,v in metrics.items() if v is not None))
            cells=''.join(f'<td data-sort="{html.escape(str(v))}">{html.escape(str(v) if not isinstance(v,float) else f"{v:.6g}")}</td>' for v in labels)
            detail=html.escape(json.dumps(dict(artifact_hash=artifact.content_hash,geometry=artifact.scientific['geometry'],categories=categories,
                fit=fit,primary_cadence=evidence.get('primary_cadence'),mechanical=mechanical,search=member['metadata'].get('search'),
                provenance=member['metadata']['provenance']),indent=2))
            primary=artifact.scientific['settings'].get('component') is not None
            complete=set(member['mechanisms'])=={'small','large'} and all(
                raw['scientific']['settings'].get('component') is None for raw in member['mechanisms'].values())
            selector=f'--{side}-library {shlex.quote(str(library_path))} --{side} {shlex.quote(member["family_id"])}'
            selection=f'<p>Source library: <code>{html.escape(str(library_path))}</code></p>'
            if primary:
                selection+='<p>Intermediate component: not eligible for pairing or adaptation.</p>'
            else:
                selection+=f'<button data-short="{html.escape(f"--{side} {shlex.quote(member['family_id'])}",quote=True)}" data-side="{side}" data-id="{html.escape(member["family_id"],quote=True)}" onclick="selectMechanism(this)">Select as {side.upper()}</button> '
                selection+=f'<button data-selector="{html.escape(selector,quote=True)}" onclick="copySelector(this)">Copy {side.upper()} selector</button><pre>{html.escape(selector)}</pre>'
            if complete:
                adapt=f'dada-research mechanism adapt path/to/source/campaign --candidate SOURCE_ID --library {shlex.quote(str(library_path))} --family-id {shlex.quote(member["family_id"])} --output path/to/thermo/study.toml'
                selection+='<p>Complete pair: usable directly with mechanism adapt.</p><pre>'+html.escape(adapt)+'</pre>'
            cadence=evidence.get('primary_cadence')
            cadence_summary='Not a primary stage'
            if cadence is not None:
                c=cadence['target']['cadence']
                lines=[f"Policy: {cadence['policy_version']}", f"Internal cadence score: {cadence['score']:.6g}",
                       f"Turnarounds (study rad): {cadence['matched_turnarounds_rad']}",
                       f"Extra reversals: {cadence['extra_crossings']}",
                       f"Long mirror asymmetry RMS: {cadence['long_mirror_asymmetry_rms']}",
                       f"Weighted long symmetry contribution: {cadence['terms']['long_mirror_asymmetry']}",
                       f"Fast/slow speed ratio: {cadence['fast_slow_speed_ratio']}",
                       f"Fast displacement fraction: {cadence['fast_displacement_fraction']}",
                       f"Target cadence: {c}"]
                lines.extend(f"{name}: {value:.6g}" for name,value in cadence['terms'].items())
                cadence_summary='<pre>'+html.escape('\n'.join(lines))+'</pre>'
            elif primary:
                cadence_summary='Archived evidence: no topology/cadence score recorded.'
            entries.append('<tr>'+cells+f'<td>{summary}</td><td>{cadence_summary}</td><td>{html.escape(json.dumps(categories))}</td>'+
                '<td><details><summary>Mechanical admissibility / design screen</summary><pre>'+html.escape(json.dumps(admissibility,indent=2))+'</pre></details>'+selection+'<details><summary>Geometry and target comparison</summary>'+_preview(artifact,target,side)+
                '<details><summary>Evidence / search provenance</summary><pre>'+detail+'</pre></details></details></td></tr>')
    headings=('Family ID','Piston','Mechanism','Position RMS','Max position error','Velocity RMS / rad','Mechanical metrics','Primary topology / cadence','Categories','Selection / preview / evidence')
    headers=''.join(f'<th><button onclick="sortCatalogue({i})">{title}</button></th>' for i,title in enumerate(headings))
    document='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Mechanism synthesis catalogue</title>
<style>body{font:14px system-ui;margin:24px;color:#202830}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd3d8;padding:8px;vertical-align:top}th{background:#eef3f6}button{font:inherit;border:0;background:none;cursor:pointer}svg{width:900px;max-width:100%;height:auto}pre{white-space:pre-wrap;max-width:900px;font-size:12px}summary{cursor:pointer}details{min-width:260px}</style>
<h1>Mechanism families for human selection</h1><p>Mechanical admissibility, target tracking quality and geometric family diversity remain separate from the internal search score.
E-projection position and velocity errors are diagnostics only, excluded from the primary score.
No thermodynamic performance has been inferred. Click headings to sort; sorting does not identify an absolute winner.</p>'''+f'<p>Target hash: <code>{html.escape(target.content_hash)}</code>. Use the Family ID with <code>mechanism visualize --family-id</code> or <code>full_local_polish --family-id</code>.</p>'+f'<p>Source library: <code>{html.escape(str(library_path))}</code></p><p>Select both sides here for the short command, or copy side selectors to combine different libraries.</p><pre id="pair-command" data-base="{html.escape(common_command,quote=True)}">Select a SMALL and a LARGE mechanism.</pre>'+f'<table><thead><tr>{headers}</tr></thead><tbody>{"".join(entries)}</tbody></table>'+'''
<script>
const selectedMechanisms={};
function selectMechanism(button){selectedMechanisms[button.dataset.side]={id:button.dataset.id,short:button.dataset.short};
const box=document.getElementById('pair-command');
box.textContent=selectedMechanisms.small&&selectedMechanisms.large?
box.dataset.base+' '+selectedMechanisms.small.short+' '+selectedMechanisms.large.short+' --output pair.json':
'Selected '+button.dataset.side.toUpperCase()+': '+button.dataset.id+'; select the other side.';}
async function copySelector(button){const text=button.dataset.selector;
try{if(navigator.clipboard){await navigator.clipboard.writeText(text);return;}}catch(error){}
const area=document.createElement('textarea');area.value=text;document.body.appendChild(area);area.select();document.execCommand('copy');area.remove();}
let previous=-1,ascending=true;function sortCatalogue(column){ascending=column===previous?!ascending:true;previous=column;let body=document.querySelector('tbody');let rows=Array.from(body.rows);rows.sort((a,b)=>{let x=a.cells[column].dataset.sort??a.cells[column].textContent,y=b.cells[column].dataset.sort??b.cells[column].textContent;let n=Number(x),m=Number(y);let comparison=Number.isFinite(n)&&Number.isFinite(m)?n-m:x.localeCompare(y);return ascending?comparison:-comparison;});rows.forEach(row=>body.appendChild(row));}</script></html>'''
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as stream: stream.write(document)
    return path

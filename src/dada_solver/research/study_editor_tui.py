"""Small optional terminal menus; all Research transformations belong to the engine."""
import json
import sys
from .study_editor import StudyEditor

LABELS = {
    'kinematics-small':'Mechanism / motion SMALL',
    'kinematics-large':'Mechanism / motion LARGE',
    'exchangers':'Exchangers', 'volumes':'Volumes', 'frequency':'Frequency',
    'charge':'Charge', 'valves':'Valves', 'external-stream':'External streams',
}


def menu_focus_index(choices, focus_key=None, previous_index=0):
    """Restore semantic focus; use the nearest surviving row if a key vanished."""
    if not choices: raise ValueError('A terminal menu requires at least one choice.')
    return next((i for i,(key,_) in enumerate(choices) if key==focus_key),
                min(max(0,previous_index),len(choices)-1))


class TerminalMenus:
    """Arrow/space/enter/escape menus, built on the optional prompt_toolkit extra."""
    def __init__(self):
        if not sys.stdin.isatty() or not sys.stdout.isatty():
            raise ValueError('Interactive study editing requires a terminal. Use study release/freeze for batch editing.')
        try:
            from prompt_toolkit import Application, PromptSession
            from prompt_toolkit.key_binding import KeyBindings
            from prompt_toolkit.layout import Layout, Window
            from prompt_toolkit.layout.controls import FormattedTextControl
            from prompt_toolkit.data_structures import Point
        except ImportError as error:
            raise ValueError('Study editor requires the optional TUI extra: pip install "dada-engine-solver[tui]".') from error
        self.Application,self.KeyBindings,self.Layout,self.Window,self.Control=Application,KeyBindings,Layout,Window,FormattedTextControl
        self.Point=Point
        self.prompt=PromptSession()
        self.menu_positions={}

    def menu(self, title, choices, *, toggle=False, focus_key=None):
        """Return a controller action; no scientific values are changed here."""
        identity=title.split('\n',1)[0]
        old_key,old_index,scroll=self.menu_positions.get(identity,(None,0,0))
        selected=menu_focus_index(choices,focus_key if focus_key is not None else old_key,old_index)
        bindings=self.KeyBindings()
        @bindings.add('up')
        def up(event):
            nonlocal selected
            selected=(selected-1)%len(choices)
        @bindings.add('down')
        def down(event):
            nonlocal selected
            selected=(selected+1)%len(choices)
        @bindings.add('enter')
        def enter(event): event.app.exit(result=('open',choices[selected][0]))
        @bindings.add(' ')
        def space(event):
            if toggle: event.app.exit(result=('toggle',choices[selected][0]))
        @bindings.add('escape')
        @bindings.add('c-c')
        @bindings.add('c-d')
        def back(event): event.app.exit(result=None)
        def text():
            lines=[('bold',title+'\n'),('', '↑/↓ navigate · Enter open/confirm · '+('Space toggle · ' if toggle else '')+'Esc back/cancel\n\n')]
            for i,(_,label) in enumerate(choices):
                lines.append(('reverse' if i==selected else '',('› ' if i==selected else '  ')+label+'\n'))
            return lines
        window=self.Window(self.Control(text,focusable=True,get_cursor_position=lambda:self.Point(x=0,y=title.count('\n')+3+selected)),always_hide_cursor=True)
        window.vertical_scroll=scroll
        app=self.Application(layout=self.Layout(window),key_bindings=bindings,full_screen=False,erase_when_done=True)
        result=app.run()
        self.menu_positions[identity]=(choices[selected][0],selected,window.vertical_scroll)
        return result

    def input(self, label, default=''):
        try: return self.prompt.prompt(label+': ',default=str(default))
        except (KeyboardInterrupt,EOFError): return None

    def message(self, text):
        print(text)

    def confirm(self, text):
        result=self.menu(text,[('yes','Confirm'),('no','Back')])
        return result==('open','yes')


def _scalar(ui, label, value):
    text=ui.input(label,json.dumps(value,allow_nan=False))
    if text is None: return None
    try: return json.loads(text)
    except ValueError as error: raise ValueError('Enter an exact JSON scalar or array; string choices require quotes.') from error


def _release(editor, ui, names):
    domains={}
    for name in names:
        try: editor.resolve_domain(name)
        except ValueError as error:
            ui.message(str(error))
            entry=next(p for p in editor.inventory() if p.name==name)
            if entry.kind=='choice':
                choices=_scalar(ui,f'{name} declared choices',[])
                if choices is None: return
                domains[name]=dict(kind='choice',choices=choices)
            else:
                lower=_scalar(ui,f'{name} lower',entry.value)
                if lower is None: return
                upper=_scalar(ui,f'{name} upper',entry.value)
                if upper is None: return
                domains[name]=dict(kind=entry.kind,lower=lower,upper=upper,
                    **({'encoding':'nearest_even_v1'} if entry.kind=='integer' else {'transform':'linear'}))
    editor.release(parameters=names,domains=domains)


def _parameter(editor, ui, name):
    focused=None
    while True:
        entry=next(p for p in editor.inventory() if p.name==name)
        row=editor._rows().get(name,{})
        choices=[]
        if entry.editable:
            choices.append(('value','Edit '+('initial' if entry.active else 'value')+' = '+repr(entry.value)))
        if entry.active:
            choices+=[('freeze','Freeze at exact current initial')]
            for key in ('lower','upper','transform','choices'):
                if key in row: choices.append((key,f'Edit {key} = {row[key]!r}'))
        elif entry.releasable: choices.append(('release','Activate parameter'))
        choices.append(('back','Back'))
        selected=ui.menu(f'{name} [{entry.unit}]\n{entry.source_of_value}\n{entry.reason}',choices,focus_key=focused)
        if selected is None or selected[1]=='back': return
        action=selected[1]; focused=action
        if action=='freeze': editor.freeze(parameters=[name])
        elif action=='release': _release(editor,ui,[name])
        else:
            field=('initial' if entry.active else 'value') if action=='value' else action
            current=entry.value if action=='value' else row[action]
            if field=='transform':
                selected=ui.menu('Coordinate transform',[('linear','Linear'),('log','Log (positive bounds only)')])
                if selected is None: continue
                value=selected[1]
            elif field in ('value','initial') and entry.kind=='choice':
                choices=row['choices'] if entry.active else editor.specs[name].choices
                selected=ui.menu('Declared scientific choices',[(v,repr(v)) for v in choices])
                if selected is None: continue
                value=selected[1]
            else:
                value=_scalar(ui,field,current)
                if value is None: continue
            editor.edit_parameter(name,**{field:value})


def bounds_preview_text(plan):
    lines=['BOUNDS RECALIBRATION — reference: CURRENT declared domain',
           f"Normalized half-width: {plan['half_width']!r}",
           'Parameter | Old declared range | New declared range']
    for name,change in plan['changes'].items():
        a,b=change['before'],change['after']
        lines.append(f"{name} | [{a['lower']!r}, {a['upper']!r}] | [{b['lower']!r}, {b['upper']!r}]")
    lines.extend(f'Skipped {name}: {reason}' for name,reason in plan['skipped'].items())
    lines.append(f"{len(plan['changes'])} numeric parameters selected; 0 initial values changed.")
    lines.append('Applying again uses these NEW bounds as reference and narrows cumulatively. Undo is available in the group menu.')
    return '\n'.join(lines)


def effective_bounds_text(editor):
    review=editor.review(); rows=editor._rows(); search=review['search']
    lines=['EFFECTIVE SEARCH BOUNDS',f"Global normalized radius: {search.get('radius_fraction','full declared bounds')} (ALL active parameters)"]
    shown=set()
    for group,label in LABELS.items():
        names=[n for n in editor.select(groups=[group]) if n in review['effective_bounds'] and n not in shown]
        if not names: continue
        lines.append('\n'+label+' — declared | initial | effective | clipped')
        for name in names:
            shown.add(name); row=rows[name]
            if row['kind']=='choice':
                lines.append(f"{name}: {row['choices']!r} | {row['initial']!r} | declared choices (no numeric distance)")
                continue
            from dada_solver.campaign.parameters import parameter_from_mapping
            p=parameter_from_mapping({k:v for k,v in row.items() if k!='unit'})
            for i,box in enumerate(review['effective_bounds'][name]):
                clipped=False
                if search['domain']=='local_regions_v1':
                    u=p.encode(search['regions'][i]['center'][name]); r=search['radius_fraction']
                    clipped=u-r<0 or u+r>1
                lines.append(f"{name}: [{row['lower']!r}, {row['upper']!r}] | {row['initial']!r} | [{box['lower']!r}, {box['upper']!r}] | {'yes' if clipped else 'no'}")
    return '\n'.join(lines)


def _recalibrate(editor,ui,group):
    half_width=_scalar(ui,'Normalized half-width (fraction of CURRENT declared domain)',.01)
    if half_width is None: return
    plan=editor.preview_recenter_bounds(groups=[group],half_width=half_width)
    selected=ui.menu(bounds_preview_text(plan),[('apply','Apply'),('back','Back')])
    if selected is None or selected[1]!='apply': return
    editor.recenter_bounds(groups=[group],half_width=half_width)
    if not plan['changes']: return
    mode,radius=editor.search_request or ('local' if editor.raw['search']['domain']=='local_regions_v1' else 'global',editor.raw['search'].get('radius_fraction',1.))
    if mode=='local' and radius!=1. and ui.confirm(
            f'Current Sobol radius: {radius!r}. Set it to 1.0 to explore the full newly declared windows? '
            'This affects ALL active parameters, including unchanged hardware; bounds themselves are not altered.'):
        editor.configure_search('local',radius=1.,recenter=editor.recenter)


def edit_group(editor,ui,group):
    """Stay in this group after actions, preserving semantic parameter focus."""
    focused=None
    while True:
        names=editor.select(groups=[group]); members=[p for p in editor.inventory() if p.name in names]
        items=[(p.name,f'[{"x" if p.active else " " if p.releasable else "-"}] {p.name} = {p.value!r} [{p.unit}]') for p in members]
        items.append(('recalibrate','Recalibrate bounds for this group'))
        if editor.bounds_history: items.append(('undo-bounds','Undo last bounds recalibration (session-wide)'))
        items.append(('back','Back'))
        selected=ui.menu(f"{LABELS[group]}\nActive: {sum(p.active for p in members)} / {len(members)}",items,toggle=True,focus_key=focused)
        if selected is None or selected[1]=='back': return
        action,name=selected; focused=name
        try:
            if name=='recalibrate':
                if action=='open': _recalibrate(editor,ui,group)
            elif name=='undo-bounds':
                if action=='open' and ui.confirm('Restore the domains preceding the last bounds recalibration? Scheduler radius is unchanged.'):
                    editor.undo_recenter_bounds()
            elif action=='toggle':
                entry=next(p for p in members if p.name==name)
                if entry.active: editor.freeze(parameters=[name])
                elif entry.releasable: _release(editor,ui,[name])
                else: ui.message(entry.reason)
            else: _parameter(editor,ui,name)
        except (ValueError,TypeError,KeyError,OSError,RuntimeError) as error:
            ui.message(f'Study edit error: {error}')


def _search(editor,ui):
    selected=ui.menu('SEARCH CONFIGURATION — radius affects ALL active parameters',[
        ('local','Local Sobol, center-first (normalized radius)'),
        ('global','Full declared parameter bounds'),('bounds','Show effective bounds'),('back','Back')])
    if selected is None or selected[1]=='back': return
    if selected[1]=='bounds':
        ui.message(effective_bounds_text(editor)); return
    mode=selected[1]
    radius=.05
    if mode=='local':
        radius=_scalar(ui,'Normalized local radius',editor.raw['search'].get('radius_fraction',.05))
        if radius is None: return
    recenter=False
    if len(editor.raw['search'].get('regions',[]))>1:
        recenter=ui.confirm('Replace multiple regions with one UNEVALUATED center at current initials?' if mode=='local' else 'Remove multiple local regions and explore full declared bounds?')
        if not recenter: return
    editor.configure_search(mode,radius=radius,recenter=recenter)


def _execution(editor,ui):
    choices=[('budget','Default budget'),('max_candidates','Maximum candidates'),
             ('seed','Sobol seed'),('scramble','Scramble'),('back','Back')]
    selected=ui.menu('EXECUTION\n'+json.dumps(dict(editor.raw['execution'],seed=editor.raw['search']['seed'],scramble=editor.raw['search']['scramble']),indent=2),choices)
    if selected is None or selected[1]=='back': return
    field=selected[1]
    if field=='budget': value=ui.input('Default budget',editor.raw['execution']['default_budget'])
    else:
        current=(editor.raw['execution']['default_max_candidates'] if field=='max_candidates' else editor.raw['search'][field])
        value=_scalar(ui,field,current)
    if value is not None: editor.edit_execution(**{field:value})



def review_text(review):
    """Readable confirmation summary; detailed domains remain on their own page."""
    lines=['STUDY EDIT — REVIEW', f"Source: {review['source']}",
           f"Active parameters: {review['active_before']} → {review['active_after']}"]
    for key,label in (('activated_parameters','Newly released'),('frozen_parameters','Frozen')):
        lines.append(label+':')
        lines.extend('  '+name for name in review[key])
        if not review[key]: lines.append('  None')
    lines.append('Changed values:')
    for name,change in review['changed_initial_values'].items():
        lines.append(f"  {name}: {change['before']!r} → {change['after']!r}")
    if not review['changed_initial_values']: lines.append('  None')
    lines.append(f"Changed declared domains: {len(review['changed_bounds'])} (Search configuration → Show effective bounds)")
    search=review['search']
    if search['domain']=='local_regions_v1':
        lines.append(f"Search: local Sobol · normalized radius {search['radius_fraction']!r} for ALL active parameters · {len(search['regions'])} region(s)")
        lines.append('Centers: exact physical values shown in Search configuration; new centers are unevaluated.')
    else: lines.append('Search: full declared parameter bounds')
    lines.append(f"Sobol seed: {search['seed']} · scramble: {search['scramble']}")
    for name,change in review['execution_configuration_change'].items():
        lines.append(f"Execution {name}: {change['before']!r} → {change['after']!r}")
    for side,change in review['artifact_changes'].items():
        lines.append(f"Artifact {side}: {change['before']} → {change['after']}")
    lines.extend('Warning: '+warning for warning in review['warnings'])
    lines.append('Objectives, constraints, physics policies and unedited scientific values are preserved. No thermodynamic evaluation.')
    return '\n'.join(lines)

def run_editor(source, *, ui=None, editor=None):
    """Testable controller: menus choose actions, the pure engine applies them."""
    ui=TerminalMenus() if ui is None else ui
    editor=StudyEditor(source) if editor is None else editor
    focused=None
    while True:
        inventory=editor.inventory(); active=sum(p.active for p in inventory)
        choices=[]
        for group,label in LABELS.items():
            names=editor.select(groups=[group]); members=[p for p in inventory if p.name in names]
            eligible=[p for p in members if p.releasable]
            if not members: continue
            count=sum(p.active for p in eligible)
            marker='x' if eligible and count==len(eligible) else '~' if count else ' '
            locked=sum(not p.releasable for p in members)
            choices.append((group,f'[{marker}] {label:<30} {count} / {len(eligible)}'+(f' · {locked} fixed by contract' if locked else '')))
        choices += [('individual','Edit individual parameters'),('search','Search configuration'),
                    ('execution','Execution settings'),('review','Review changes'),('save','Save as new study'),('cancel','Cancel')]
        result=ui.menu(f'DADA RESEARCH — STUDY EDITOR\nSource: {source}\nActive: {active} · available: {len(inventory)}',choices,toggle=True,focus_key=focused)
        if result is None or result[1]=='cancel':
            ui.message('Cancelled; no files written.'); return None
        action,selection=result; focused=selection
        if action=='toggle' and selection not in LABELS:
            continue  # Space changes activity only; Enter opens action pages.
        try:
            if selection in LABELS:
                names=editor.select(groups=[selection],releasable_only=True)
                if action=='toggle':
                    if names and all(p.active for p in inventory if p.name in names): editor.freeze(parameters=names)
                    elif names: _release(editor,ui,names)
                    else: ui.message('This group has no releasable coordinates under the current Research contract.')
                else:
                    edit_group(editor,ui,selection)
            elif selection=='individual':
                query=ui.input('Filter parameter name (empty = all)','')
                if query is None: continue
                selected=ui.menu('INDIVIDUAL PARAMETERS',[(p.name,f'{p.name} = {p.value!r}') for p in inventory if query in p.name]+[('back','Back')])
                if selected and selected[1]!='back': _parameter(editor,ui,selected[1])
            elif selection=='search': _search(editor,ui)
            elif selection=='execution': _execution(editor,ui)
            elif selection in ('review','save'):
                try: review=editor.review()
                except ValueError as error:
                    if 'multiple local regions' not in str(error): raise
                    if not ui.confirm(str(error)+' Recenter?'): continue
                    editor.recenter=True; review=editor.review()
                ui.message(review_text(review))
                if selection=='save':
                    output=ui.input('New study path','')
                    if not output: continue
                    if ui.confirm(f'Save validated portable study to {output}? Source files will remain unchanged.'):
                        path=editor.save(output); ui.message(f'Created {path}; no integration started.'); return path
        except (ValueError,TypeError,KeyError,OSError,RuntimeError) as error:
            ui.message(f'Study edit error: {error}')

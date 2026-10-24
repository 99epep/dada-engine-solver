"""Read-only declared and effective Research search windows."""
import copy
import math
import shutil
import textwrap
from dada_solver.campaign.parameters import parameter_from_mapping, ChoiceParameter
from .hardware_retuning import GROUPS, _effective_bounds


GROUP_LABELS = {
    'kinematics-small':'Mechanism / motion SMALL',
    'kinematics-large':'Mechanism / motion LARGE',
    'exchangers':'Exchangers', 'volumes':'Volumes', 'frequency':'Frequency',
    'charge':'Charge', 'valves':'Valves', 'external-stream':'External streams',
}


def effective_search_bounds(rows, search):
    """Reuse normalized-domain calculation; preserve legacy scheduler choices."""
    result={}
    for row in rows:
        if 'initial' not in row: continue
        if search['domain']!='local_regions_v1':
            result[row['name']]=[copy.deepcopy(row)];continue
        boxes=[]
        for region in search['regions']:
            box=_effective_bounds(row,region['center'][row['name']],search['radius_fraction'])
            if row['kind']=='choice' and search.get('choice_scope')!='declared_choices':
                p=parameter_from_mapping({k:v for k,v in row.items() if k!='unit'})
                u=p.encode(region['center'][row['name']]);r=search['radius_fraction']
                low,high=max(0.,u-r),min(1.,u+r)
                first=p.choices.index(p.decode(low))
                last=p.choices.index(p.decode(math.nextafter(high,low)))
                box['choices']=list(p.choices[first:last+1])
                box['policy']='legacy_encoded_window'
            boxes.append(box)
        result[row['name']]=boxes
    return result


def search_window(rows, search):
    effective=effective_search_bounds(rows,search)
    result=[]
    local=search['domain']=='local_regions_v1'
    for row in rows:
        if 'initial' not in row: continue
        name=row['name']
        group=next((g for g,prefixes in GROUPS.items() if name.startswith(prefixes)),
                   'kinematics-'+name.split('.')[1] if name.startswith('kinematics.') else 'other')
        parameter=parameter_from_mapping({k:v for k,v in row.items() if k!='unit'})
        boxes=[]
        for index,box in enumerate(effective[name]):
            region=search['regions'][index] if local else None
            clipped=False
            if local and not isinstance(parameter,ChoiceParameter):
                u=parameter.encode(region['center'][name]);r=search['radius_fraction']
                clipped=u-r<0 or u+r>1
            boxes.append(dict(region=region['id'] if region else None,
                              center=region['center'][name] if region else row['initial'],
                              clipped=clipped,**box))
        result.append(dict(group=group,name=name,unit=row['unit'],kind=row['kind'],initial=row['initial'],
            declared={k:row[k] for k in ('lower','upper','choices','transform','encoding') if k in row},windows=boxes))
    return dict(domain=search['domain'],radius_fraction=search.get('radius_fraction'),
                choice_scope=search.get('choice_scope','legacy_encoded_window') if local else 'declared_choices',parameters=result,integration_started=False)


def _display_number(value):
    if isinstance(value,float): return format(value,'.5g')
    return repr(value)


def _display_pair(values, width, *, digits=5):
    """Keep both numeric bounds on one line; round display only to fit."""
    candidates=[]
    for precision in range(digits,0,-1):
        low,high=(format(value,f'.{precision}g') for value in values)
        text=low+('  ' if high.startswith('-') else '   ')+high
        candidates.append(text)
        if len(text)<=width: return text
    return min(candidates,key=len)


def _aligned_table(headers, rows, width):
    """Keep numeric pairs intact and wrap other cells within fixed columns."""
    widths=[13,5,11,14,14,11,8]
    minimum=[len(header)+(column>0) for column,header in enumerate(headers)]
    for column in (3,4):
        minimum[column]=max([minimum[column],*[len(_display_pair(row[column],0))+(row[column][0]>=0)
                           for row in rows if isinstance(row[column],tuple)]])
    widths=[max(size,low) for size,low in zip(widths,minimum)]
    while sum(widths)+len(widths)-1>width:
        index=max(range(len(widths)),key=lambda i:widths[i]-minimum[i])
        if widths[index]<=minimum[index]: break
        widths[index]-=1
    extra=max(0,width-sum(widths)-len(widths)+1)
    # Prioritize full-precision bound pairs, then actual name/unit lengths. Do
    # not turn all remaining screen space into an empty parameter column.
    for column in (3,4,0,1):
        desired=max([minimum[column],*[len(_display_pair(row[column],1000))+(row[column][0]>=0) if isinstance(row[column],tuple)
                    else len(str(row[column]))+(column>0 and not str(row[column]).startswith('-')) for row in rows]])
        growth=min(extra,max(0,desired-widths[column]))
        widths[column]+=growth;extra-=growth
    lines=[]
    for index,row in enumerate([headers,*rows]):
        cells=[]
        for column,(cell,size) in enumerate(zip(row,widths)):
            if isinstance(cell,tuple):
                padding=column>0 and not format(cell[0],'.5g').startswith('-')
                parts=[_display_pair(cell,size-padding)]
            else:
                if str(cell).startswith('-') and len(str(cell))<=size: parts=[str(cell)]
                else: parts=textwrap.wrap(str(cell),size-(column>0),break_long_words=True,break_on_hyphens=False) or ['']
            cells.append([(' ' if column>0 and not part.startswith('-') else '')+part for part in parts])

        for line in range(max(map(len,cells))):
            lines.append('|'.join((cell[line] if line<len(cell) else '').ljust(size)
                                  for cell,size in zip(cells,widths)).rstrip())
        if index==0: lines.append('+'.join('-'*size for size in widths))
    return lines


def format_search_window(window, *, width=None):
    width=max(55,width if width is not None else shutil.get_terminal_size(fallback=(80,24)).columns)
    radius=_display_number(window['radius_fraction']) if window['radius_fraction'] is not None else 'full declared bounds'
    lines=['SEARCH WINDOW',f"Domain: {window['domain']}",
           f"Normalized radius: {radius} (ALL active parameters)",
           'Display: up to 5 significant digits; JSON retains full precision.']
    lines=[part for line in lines for part in textwrap.wrap(line,width)]
    if not window['parameters']: lines.append('No active parameters; use evaluate rather than run.')
    for group in dict.fromkeys(row['group'] for row in window['parameters']):
        members=[row for row in window['parameters'] if row['group']==group]
        lines.extend(['',GROUP_LABELS.get(group,group)])
        for index in range(len(members[0]['windows'])):
            region=members[0]['windows'][index]['region']
            if region: lines.extend(textwrap.wrap('Region: '+region,width))
            rows=[];notes=[]
            for row in members:
                box=row['windows'][index];declared=row['declared'];name=row['name']
                if name.startswith('kinematics.'): name=name.split('.',2)[-1]
                elif name.startswith(('microtube.','thermal.')): name=name.split('.',1)[-1]
                display=repr if row['kind']=='choice' else _display_number
                if row['kind']=='choice':
                    bounds=repr(declared['choices']);explored=repr(box['choices'])
                    notes.append(name+': no numeric distance; '+box.get('policy','declared_choices').replace('_',' '))
                else:
                    bounds=(declared['lower'],declared['upper'])
                    explored=(box['lower'],box['upper'])
                unit='r' if row['unit']=='crank_radius' else row['unit']
                rows.append([name,unit,display(row['initial']),bounds,explored,
                             display(box['center']),'yes' if box['clipped'] else 'no'])
            lines.extend(_aligned_table(['Parameter','Unit','Initial','Declared','Explored','Center','clipped'],rows,width))
            if any(row['unit']=='crank_radius' for row in members): notes.append('r = crank_radius')
            for note in dict.fromkeys(notes): lines.extend(textwrap.wrap(note,width))
    return '\n'.join(lines)

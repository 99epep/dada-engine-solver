"""Source locations for semantic TOML validation, without changing schemas."""
from contextlib import contextmanager
from pathlib import Path
import re
import tomllib


class StudyValidationError(ValueError):
    """A validation failure annotated with the source declaration."""


def _key_path(fragment):
    data=tomllib.loads(fragment)
    result=[]
    while isinstance(data,dict) and data:
        key,data=next(iter(data.items()));result.append(key)
        if isinstance(data,list): data=data[0]
    return tuple(result)


def _multiline_state(line,state):
    """Track TOML string delimiters so embedded headers cannot fake locations."""
    i=0
    while i<len(line):
        if state is not None:
            end=line.find(state,i)
            if end<0: return state
            if state=='"""' and (len(line[:end])-len(line[:end].rstrip('\\')))%2:
                i=end+3;continue
            i=end+3;state=None;continue
        if line[i]=='#': break
        if line[i:i+3] in ('"""',"'''"):
            state=line[i:i+3];i+=3;continue
        if line[i] in ('"',"'"):
            quote=line[i];i+=1
            while i<len(line):
                if quote=='"' and line[i]=='\\': i+=2;continue
                if line[i]==quote: i+=1;break
                i+=1
        else: i+=1
    return state


class SourceLocations:
    def __init__(self,path):
        self.lines=Path(path).read_text().splitlines()
        self.locations={};table=();counts={};state=None
        for number,line in enumerate(self.lines,1):
            inside=state is not None
            state=_multiline_state(line,state)
            if inside: continue
            stripped=line.strip()
            if stripped.startswith('['):
                header=re.match(r'^(\[\[.*?\]\]|\[.*?\])\s*(?:#.*)?$',stripped)
                if header:
                    try: table=_key_path(header[1])
                    except tomllib.TOMLDecodeError: continue
                    if header[1].startswith('[['):
                        index=counts.get(table,0);counts[table]=index+1;table=(*table,index)
                    self.locations[table]=number
            else:
                assignment=re.match(r'^\s*((?:"(?:[^"\\]|\\.)*"|\'[^\']*\'|[A-Za-z0-9_-]+)(?:\s*\.\s*(?:"(?:[^"\\]|\\.)*"|\'[^\']*\'|[A-Za-z0-9_-]+))*)\s*=',line)
                if assignment:
                    key=_key_path(assignment[1]+' = 0')
                    self.locations[(*table,*key)]=number

    def locate(self,address,message):
        # Prefer a key explicitly named by the error within the supplied scope.
        direct=[key for key in self.locations if key[:-1]==address and isinstance(key[-1],str)]
        mentioned=[key for key in direct if re.search(r'(?<!\w)'+re.escape(key[-1])+r'(?!\w)',message)]
        if mentioned: address=max(mentioned,key=lambda key:len(key[-1]))
        elif (*address,'name') in self.locations: address=(*address,'name')
        while address not in self.locations and address: address=address[:-1]
        return self.locations.get(address),address


@contextmanager
def validation_location(path,address=(),*,label=None):
    """Enrich failures only; successful validation performs no extra file reads."""
    try:
        yield
    except StudyValidationError:
        raise
    except (ValueError,TypeError,KeyError,ArithmeticError) as error:
        locations=SourceLocations(path)
        line,resolved=locations.locate(tuple(address),str(error))
        context=label or '.'.join(map(str,address)) or 'study'
        detail='.'.join(map(str,resolved))
        if line is None and isinstance(error,tomllib.TOMLDecodeError):
            match=re.search(r'at line (\d+)',str(error))
            if match: line=int(match[1])
        if line is None:
            raise StudyValidationError(f'{path}: {context}: {error}') from error
        snippet=locations.lines[line-1].strip()
        related=''
        if label and label.endswith('.geometry') and len(address)==2:
            raw=tomllib.loads(Path(path).read_text())
            declarations=[]
            for index,row in enumerate(raw.get('parameters',())):
                name=row.get('name','')
                if name.startswith(f'kinematics.{address[1]}.'):
                    at=locations.locations.get(('parameters',index,'name'))
                    declarations.append(f'    {name}: {path}:{at}' if at is not None else f'    {name}')
            if declarations: related='\n  Coupled geometry parameters (not a single identified culprit):\n'+'\n'.join(declarations)
        raise StudyValidationError(f'{path}:{line}: {context} [{detail}]: {error}\n  {line} | {snippet}{related}') from error

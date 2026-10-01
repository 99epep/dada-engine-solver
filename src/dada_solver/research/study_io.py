"""Small deterministic TOML writer for generated study templates."""
import json
import math


def dumps(data):
    def value(x):
        if isinstance(x,str): return json.dumps(x,ensure_ascii=False)
        if type(x) is bool: return str(x).lower()
        if isinstance(x,(int,float)) and not isinstance(x,bool):
            if not math.isfinite(x): raise ValueError('TOML numeric values must be finite.')
            return repr(x)
        if isinstance(x,(list,tuple)) and not any(isinstance(v,dict) for v in x): return '['+', '.join(value(v) for v in x)+']'
        raise ValueError(f'Unsupported TOML scalar: {type(x).__name__}')
    lines=[]
    def table(row,path=(),array=False):
        if path: lines.extend(['',('[[ ' if array else '[ ')+'.'.join(json.dumps(k) for k in path)+(' ]]' if array else ' ]')])
        for k,v in row.items():
            if not isinstance(v,dict) and not (isinstance(v,list) and v and isinstance(v[0],dict)):
                lines.append(json.dumps(k)+' = '+value(v))
        for k,v in row.items():
            if isinstance(v,dict): table(v,(*path,k))
            elif isinstance(v,list) and v and isinstance(v[0],dict):
                for item in v: table(item,(*path,k),True)
    table(data)
    return '# Dada-Engine Research: fixed values use value; active coordinates use initial with lower/upper or choices.\n'+'\n'.join(lines)+'\n'

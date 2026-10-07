"""Validation errors identify source declarations without thermodynamic replay."""
import re
import tomllib
import pytest
from dada_solver.research.cli import main
from dada_solver.research.presets import initialize_kinematics
from dada_solver.research.study_io import dumps
from dada_solver.research.validation_errors import SourceLocations


def invalid_study(tmp_path,mutate):
    path=initialize_kinematics(tmp_path/'study.toml','slider_crank','slider_crank')
    raw=tomllib.loads(path.read_text());mutate(raw);path.write_text(dumps(raw))
    return path


def error_message(path,capsys):
    with pytest.raises(SystemExit) as error: main(['validate',str(path)])
    assert error.value.code==2
    return capsys.readouterr().err


@pytest.mark.parametrize('field,value',[('value',False),('value',float('nan')),('unit','invalid')])
def test_parameter_name_value_and_exact_line(tmp_path,capsys,monkeypatch,field,value):
    from dada_solver.campaign.evaluator import MachineEvaluator
    monkeypatch.setattr(MachineEvaluator,'evaluate_with_control',lambda *a,**k:pytest.fail('Unexpected integration'))
    name='kinematics.small.rod_over_crank'
    def mutate(raw):
        next(row for row in raw['parameters'] if row['name']==name)[field]=value
    if isinstance(value,float) and value!=value:
        path=invalid_study(tmp_path,lambda raw:None)
        text=path.read_text();start=text.index('\"name\" = \"'+name+'\"');at=text.index('\"value\" = ',start);end=text.index('\n',at)
        path.write_text(text[:at]+'\"value\" = nan'+text[end:])
    else:
        path=invalid_study(tmp_path,mutate)
    raw=tomllib.loads(path.read_text());index=next(i for i,row in enumerate(raw['parameters']) if row['name']==name)
    line=SourceLocations(path).locations[('parameters',index,field)]
    message=error_message(path,capsys)
    assert f'{path}:{line}:' in message and name in message
    assert path.read_text().splitlines()[line-1].strip() in message


def test_active_bounds_error_names_coordinate(tmp_path,capsys):
    name='kinematics.small.rod_over_crank'
    def mutate(raw):
        row=next(row for row in raw['parameters'] if row['name']==name)
        row.pop('value');row.update(kind='continuous',transform='linear',lower=8.,upper=2.,initial=4.)
    path=invalid_study(tmp_path,mutate)
    message=error_message(path,capsys)
    assert name in message and 'bounds must be ordered' in message
    assert re.search(re.escape(str(path))+r':\d+:',message)


@pytest.mark.parametrize('section,key,value',[('numerical','maximum_cycles',0),('kinematics','count',3)])
def test_table_error_locates_key(tmp_path,capsys,section,key,value):
    def mutate(raw):
        if section=='kinematics': raw[section]['small']={'family':'free_spline','count':value,'representation':'shape_coordinates'}
        else: raw[section][key]=value
    path=invalid_study(tmp_path,mutate)
    address=(section,'small',key) if section=='kinematics' else (section,key)
    line=SourceLocations(path).locations[address]
    assert f'{path}:{line}:' in error_message(path,capsys)


def test_mechanical_initial_failure_reports_metric_value_limit_and_line(tmp_path,capsys):
    def mutate(raw):
        raw['mechanical_constraints']=[dict(side='small',metric='minimum_rod_axis_cosine',relation='minimum',limit=2.,unit='1')]
    path=invalid_study(tmp_path,mutate)
    line=SourceLocations(path).locations[('mechanical_constraints',0,'limit')]
    message=error_message(path,capsys)
    assert f'{path}:{line}:' in message
    assert 'small.minimum_rod_axis_cosine.minimum' in message
    assert 'value=' in message and 'limit=2.0' in message


def test_location_map_quoted_headers_comments_multiline_and_duplicate_rows(tmp_path):
    path=tmp_path/'source.toml'
    path.write_text('''[study]
purpose = """
[[parameters]]
name = "fake"
"""
[[ "parameters" ]] # first
name = "first"
value = 1
[[parameters]]
name = "second"
value = 2
''')
    locations=SourceLocations(path)
    assert locations.locations[('parameters',0,'name')]==7
    assert locations.locations[('parameters',1,'value')]==11
    assert len([key for key in locations.locations if key[0]=='parameters' and key[-1]=='name'])==2


def test_toml_parse_error_includes_filename(tmp_path,capsys):
    path=tmp_path/'bad.toml';path.write_text('schema_version = 3\n[ broken\n')
    message=error_message(path,capsys)
    assert str(path) in message and 'line 2' in message


def test_coupled_geometry_error_lists_relevant_parameters(tmp_path,capsys):
    def mutate(raw):
        next(row for row in raw['parameters'] if row['name']=='kinematics.small.rod_over_crank')['value']=.5
    path=invalid_study(tmp_path,mutate)
    message=error_message(path,capsys)
    assert 'small.geometry' in message and 'Coupled geometry parameters' in message
    assert 'kinematics.small.rod_over_crank' in message
    assert 'kinematics.small.offset_over_crank' in message

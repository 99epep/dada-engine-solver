"""HTML retains ranked feasible details while global evidence remains complete."""
import copy
import json
import math
import pytest
from dada_solver.research.report import inspect,render_html
from dada_solver.campaign.report import elite_records
from tests.test_research_cockpit import campaign


def embedded(data,path):
    text=render_html(data,path).read_text()
    return json.loads(text.split('<script id="data" type="application/json">')[1].split('</script>')[0])


@pytest.mark.parametrize('count,feasible',[(0,0),(4,4),(21,2),(100,70),(4096,1000)])
def test_detail_cap_ranking_and_complete_statistics(tmp_path,count,feasible):
    c,_=campaign(tmp_path,1);data=inspect(c.directory);seed=data['records'][0]
    records=[]
    for i in range(count):
        r=copy.deepcopy(seed);r['candidate_id']=f'{i:064x}'
        r['status']='feasible' if i<feasible else 'invalid_exchanger'
        r['objective']['value']=-i;r['diagnostics']={'sentinel':'detail_'+str(i)}
        records.append(r)
    data['records']=records;data['selected']=records
    data['status_counts']={'feasible':feasible,'invalid_exchanger':count-feasible}
    data['cockpit']={'sentinel':'all-attempt evidence'}
    before=copy.deepcopy(data);result=embedded(data,tmp_path/'report.html')
    cap=math.ceil(count*.1);expected=elite_records(records,cap)
    assert [r['candidate_id'] for r in result['selected']]==[r['candidate_id'] for r in expected]
    assert len(result['records'])==count
    assert result['status_counts']==data['status_counts'] and result['cockpit']==data['cockpit']
    assert all(r['status']=='feasible' for r in result['selected'])
    assert result['comparison_default_ids']==[r['candidate_id'] for r in expected[:2]]
    assert data==before


def test_duplicates_are_not_multiple_retained_candidates(tmp_path):
    c,_=campaign(tmp_path,3);data=inspect(c.directory)
    data['records']=[copy.deepcopy(r) for r in data['records'] for _ in range(20)]
    result=embedded(data,tmp_path/'report.html')
    assert len(result['selected'])==3 and len(result['records'])==60

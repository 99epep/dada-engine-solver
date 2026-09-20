"""Disk-cache experiment tests, separate from backend numerical validation."""
import json
import os
from pathlib import Path
import subprocess
import sys
import pytest


def test_cache_namespace_includes_source_and_runtime():
    from dada_solver.wall_backend import cache_namespace
    identity=dict(source_sha256='a',python='3',llvm=[15],settings=dict(profile=False))
    a=cache_namespace(identity,b'equations one')
    assert a!=cache_namespace(identity,b'equations two')
    assert a!=cache_namespace(dict(identity,llvm=[16]),b'equations one')
    assert a==cache_namespace(dict(identity,settings=dict(profile=True)),b'equations one')


def test_disk_cache_reused_in_fresh_process(tmp_path):
    pytest.importorskip('numba')
    code='''import json,sys
from tests.test_wall_backend import values
from tests.test_solver_acceleration import variable_wrapper
from dada_solver.wall_backend import WallBackendSettings,WallRHS
w=variable_wrapper();r=WallRHS(w,WallBackendSettings('numba',disk_cache=True,cache_directory=sys.argv[1]))
y=r(0.,values(w));s=r.snapshot()
print(json.dumps(dict(values=y.tolist(),hits=s['dispatcher_cache_hits'],misses=s['dispatcher_cache_misses'],
    namespace=s['cache_namespace'],first=s['first_call_seconds'])))
'''
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONPATH='src')
    def run(): return json.loads(subprocess.check_output([sys.executable,'-c',code,str(tmp_path)],env=env,text=True))
    first=run();second=run()
    assert first['misses']==1 and first['hits']==0
    assert second['hits']==1 and second['misses']==0
    assert first['values']==second['values']
    assert first['namespace']==second['namespace']
    assert list(tmp_path.rglob('*.nbc'))


def test_unwritable_cache_falls_back_to_uncached_compilation(tmp_path):
    pytest.importorskip('numba')
    from dada_solver.wall_backend import WallRHS,WallBackendSettings
    from tests.test_solver_acceleration import variable_wrapper
    # A file where a directory is needed is deterministic even as root.
    path=tmp_path/'file';path.write_text('not a directory')
    r=WallRHS(variable_wrapper(),WallBackendSettings('numba',disk_cache=True,cache_directory=str(path)))
    assert r.snapshot()['actual_backend']=='numba'
    assert r.snapshot()['cache_error'] and r.snapshot()['cache_directory'] is None

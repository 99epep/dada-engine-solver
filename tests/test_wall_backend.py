"""Production backend identity, shared physics, optional dependency and safeguards."""
from dataclasses import replace
import math
import subprocess
import sys
import numpy as np
import pytest
from dada_solver.wall_backend import WallBackendSettings,WallRHS
from dada_solver import numerical_primitives as numeric
from dada_solver.exchangers.air_wall import AirWallMotor
from dada_solver.exchangers.gas_transport import DiluteGasTransport
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.state import UniformCharge
from tests.test_solver_acceleration import variable_wrapper


def values(w,angle=0.,temperature=350.,pressure=2e5):
    x=UniformCharge(pressure,temperature).create_state(w.model.gas,w.model.volumes(angle)).as_array()
    return np.r_[x,w.heat_in.wall_capacity_j_k*450,w.heat_out.wall_capacity_j_k*330,np.zeros(5)]


def with_species(w,species):
    tr=DiluteGasTransport(species=species)
    def link(x):return replace(x,gas_model=replace(x.gas_model,transport=tr))
    factor=2.5 if species in ('argon','helium') else 3.5
    gas=CaloricallyPerfectGas(tr.gas_constant,factor*tr.gas_constant,(factor-1)*tr.gas_constant)
    m=replace(w.model,gas=gas,small_cold_link=link(w.model.small_cold_link),large_hot_link=link(w.model.large_hot_link),
        cold_large_valve=replace(w.model.cold_large_valve,flow_model=link(w.model.cold_large_valve.flow_model)),
        hot_small_valve=replace(w.model.hot_small_valve,flow_model=link(w.model.hot_small_valve.flow_model)))
    def wall(x): return replace(x,gas_film=replace(x.gas_film,model=replace(x.gas_film.model,transport=tr)))
    return replace(w,model=m,heat_in=wall(w.heat_in),heat_out=wall(w.heat_out))


@pytest.mark.parametrize('species',numeric.SPECIES)
@pytest.mark.parametrize('temperature',[200.,298.15,499.999,500.,598.15,699.999,700.,1000.])
def test_species_properties_and_rhs(species,temperature):
    pytest.importorskip('numba')
    w=with_species(variable_wrapper(),species);r=WallRHS(w,WallBackendSettings('numba'))
    x=values(w,temperature=temperature)
    x[:2]*=1.00001  # Nonzero transport at fixed upstream temperature.
    try: expected=w.derivative(0.,x)
    except ValueError as error:
        # Reconstructing U/(m*Cv) at an exact boundary can round outside it.
        with pytest.raises(type(error)) as actual: r(0.,x)
        assert str(actual.value)==str(error)
    else:
        np.testing.assert_allclose(r(0.,x),expected,rtol=2e-11,atol=1e-11)
        assert r.snapshot()['fallback_calls']==0


@pytest.mark.parametrize('delta',[-1e-5,-1e-12,0.,1e-12,1e-5])
def test_breakpoints_reversal_zero_and_conservation(delta):
    pytest.importorskip('numba')
    w=variable_wrapper();r=WallRHS(w,WallBackendSettings('numba'));original=AirWallMotor.derivative
    for edge in (0.,*w.model.kinematics.breakpoint_angles(),2*math.pi):
        for a in (edge-1e-10,edge,edge+1e-10):
            x=values(w,a);x[1]*=1+delta;x[7]*=1-delta
            y=r(a,x);np.testing.assert_allclose(y,w.derivative(a,x),rtol=2e-11,atol=1e-11)
            assert abs(sum(y[:8:2]))<1e-15
            assert abs(sum(y[1:8:2])+sum(y[8:10])-sum(y[10:12])+y[14])<1e-9
    assert AirWallMotor.derivative is original
    assert r.snapshot()['fallback_calls']==0


@pytest.mark.parametrize('temperature,pressure',[(2000.,2e5),(350.,100.)])
def test_error_fallback_keeps_authoritative_message(temperature,pressure):
    pytest.importorskip('numba')
    w=variable_wrapper();r=WallRHS(w,WallBackendSettings('numba'));x=values(w,temperature=temperature,pressure=pressure)
    with pytest.raises(ValueError) as a:w.derivative(0.,x)
    with pytest.raises(type(a.value)) as b:r(0.,x)
    assert str(a.value)==str(b.value)
    assert r.snapshot()['fallback_calls']==1


def test_unsupported_family_is_explicit():
    pytest.importorskip('numba')
    w=variable_wrapper();w=replace(w,heat_in=replace(w.heat_in,gas_film=None))
    r=WallRHS(w,WallBackendSettings('numba'));x=values(w)
    np.testing.assert_array_equal(r(0.,x),w.derivative(0.,x))
    assert r.snapshot()['fallback_calls']==1
    assert r.snapshot()['actual_backend']=='python'
    assert r.snapshot()['reason'].startswith('unsupported_model:')


def test_without_numba_in_fresh_process():
    code='''import sys,importlib.abc
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname.split('.')[0] in ('numba','llvmlite'): raise ModuleNotFoundError(fullname)
sys.meta_path.insert(0,Block())
from tests.test_solver_acceleration import variable_wrapper
from dada_solver.wall_backend import WallRHS,WallBackendSettings
w=variable_wrapper()
assert WallRHS(w).snapshot()['actual_backend']=='python'
assert 'numba' not in sys.modules
r=WallRHS(w,WallBackendSettings('numba'))
assert r.snapshot()['reason']=='numba_unavailable'
'''
    subprocess.run([sys.executable,'-c',code],check=True)


def test_research_backend_is_identity_owned(tmp_path):
    from dada_solver.research.presets import initialize_v3
    from dada_solver.research.schema import load_study, compile_study
    from tests.test_research_v3 import rewrite
    path = initialize_v3(tmp_path/'study.toml')
    rewrite(path, lambda raw: raw['numerical'].update(backend='python'))
    reference = compile_study(load_study(path))
    rewrite(path, lambda raw: raw['numerical'].update(backend='numba'))
    compiled = compile_study(load_study(path))
    assert reference.wall_backend.name == 'python'
    assert compiled.definition_id != reference.definition_id
    identity = compiled.identity['wall_backend']
    assert identity['settings']['name'] == 'numba'
    assert identity['source_sha256'] and identity['python']
    if identity['numba_available']:
        assert identity['numba'] and identity['llvm'] and identity['llvmlite']


def test_first_progress_interruption_preserves_no_endpoint():
    from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor
    from dada_solver.integration import IntegrationInterrupted
    w=variable_wrapper();records=[]
    def stop(_):raise IntegrationInterrupted('Deadline')
    r=solve_periodic_wall_motor(w,values(w)[:10],maximum_cycles=2,
        backend=WallBackendSettings('numba'),progress_callback=stop,statistics_callback=records.append)
    assert r.status=='interrupted' and r.last_complete_state is None
    assert any(x["phase"]=="integration_preflight" for x in records)


def test_interruption_retains_completed_endpoint(monkeypatch):
    pytest.importorskip('numba')
    from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor
    from dada_solver.integration import IntegrationInterrupted
    w = variable_wrapper()
    ends = []

    def integrate(self, state, *, rhs, progress_callback, **kwargs):
        # Exercise the selected RHS without a costly convergence campaign.
        rhs(0., np.r_[state, np.zeros(5)])
        progress_callback({})
        end = state.copy()
        end[8] *= 1.01  # Keep the first complete cycle non-periodic.
        trajectory = np.zeros((15, 2))
        trajectory[:10, 0] = state
        trajectory[:10, 1] = end
        return np.array([0., 2*math.pi]), trajectory

    monkeypatch.setattr(AirWallMotor, 'integrate_cycle', integrate)
    def stop(_):
        if ends:
            raise IntegrationInterrupted('After one cycle')
    result = solve_periodic_wall_motor(
        w, values(w)[:10], maximum_cycles=3,
        backend=WallBackendSettings('numba'), progress_callback=stop,
        cycle_callback=lambda cycle, end, error, item: ends.append(end.copy()),
    )
    assert result.status == 'interrupted'
    np.testing.assert_array_equal(result.last_complete_state, ends[0])
    assert result.backend_statistics['actual_backend'] == 'numba'
    assert result.backend_statistics['calls'] > 0


@pytest.mark.parametrize('name',['cuda','NUMBA',''])
def test_invalid_backend(name):
    with pytest.raises(ValueError):WallBackendSettings(name)


def test_invalid_preflight_precedes_progress_in_both_backends():
    from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor
    from dada_solver.integration import IntegrationInterrupted
    w=variable_wrapper();errors=[]
    def stop(_): raise IntegrationInterrupted('Should not precede invalid preflight')
    for name in ('python','numba'):
        with pytest.raises(ValueError) as error:
            solve_periodic_wall_motor(w,values(w,temperature=2000.)[:10],maximum_cycles=2,
                backend=WallBackendSettings(name),progress_callback=stop)
        errors.append((type(error.value),str(error.value)))
    assert errors[0]==errors[1]

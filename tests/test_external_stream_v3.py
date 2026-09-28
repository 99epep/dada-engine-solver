"""External boundaries preserve historical algebra and signed cycle accounting."""
from dataclasses import replace, asdict
import json
from pathlib import Path
import numpy as np
import pytest
from dada_solver.exchangers.air_wall import AirWallExchanger, AirWallMotor
from dada_solver.exchangers.external_stream import (ExternalFluidStream, ExternalStreamWallExchanger,
    ExternalStreamWallMachine)
from dada_solver.exchangers.wall_cycle import solve_periodic_wall_machine, wall_cycle_performance
from dada_solver.wall_backend import WallBackendSettings
from dada_solver.performance import OperatingMode, classify_cycle


def neutral(old):
    return ExternalStreamWallExchanger(old.gas_wall_conductance_w_k,old.wall_capacity_j_k,
        old.external_stream,old.gas_film)


@pytest.mark.parametrize('flow',[0.,.01,.2])
@pytest.mark.parametrize('conductance',[0.,1.,5000.])
def test_declared_air_instantaneous_parity(flow,conductance):
    old=AirWallExchanger(13.,conductance,100.,flow,1005.,400.)
    new=neutral(old)
    for t in (250.,400.,600.):
        a=old.thermal_point(t,35000.);b=new.thermal_point(t,35000.)
        for name in ('gas_heat_w','external_heat_w','wall_energy_rate_w','external_outlet_temperature_k','wall_temperature_k'):
            assert getattr(a,name)==getattr(b,name)
        assert b.external_heat_w-b.gas_heat_w==b.wall_energy_rate_w
        assert 'air_heat_w' not in new.rates(t,35000.)
    assert 'air_heat_w' in old.rates(300.,35000.)
    assert AirWallMotor is ExternalStreamWallMachine


@pytest.mark.parametrize('field,value',[('inlet_temperature_k',0.),('mass_flow_kg_s',-.1),('cp_j_kg_k',0.),('wall_conductance_w_k',-1.),('fluid',''),('mass_flow_kg_s',float('nan'))])
def test_stream_validation(field,value):
    data=dict(inlet_temperature_k=300.,mass_flow_kg_s=.1,cp_j_kg_k=4180.,wall_conductance_w_k=100.,fluid='declared liquid')
    data[field]=value
    with pytest.raises(ValueError): ExternalFluidStream(**data)


def test_pre_v3_source_trajectory_and_neutral_roundtrip():
    from dada_solver.research.schema import load_study,compile_study
    root=Path(__file__).resolve().parents[1]
    ref=json.loads((root/'tests/data/pre_v3_air_wall.json').read_text())
    d=compile_study(load_study(root/'src/dada_solver/research/data/sixbar-thermo5d.toml'))
    old=d.adapter.build(dict(d.fixed_parameters,**ref['parameters'])).build()
    new=replace(old,heat_in=neutral(old.heat_in),heat_out=neutral(old.heat_out))
    x=np.array(ref['initial'])
    np.testing.assert_allclose(new.derivative(0,x),ref['rhs'],rtol=2e-11,atol=1e-11)
    for w in (old,new):
        r=solve_periodic_wall_machine(w,x,maximum_cycles=1,settings=d.wall_numerical_settings,backend=WallBackendSettings('numba'))
        np.testing.assert_allclose(r.trajectory[:,-1],ref['end'],rtol=2e-11,atol=1e-11)
        sampled=np.array([np.interp(ref['angles'],r.angles,v) for v in r.trajectory])
        np.testing.assert_allclose(sampled,ref['trajectory'],rtol=2e-11,atol=1e-11)
        p=asdict(wall_cycle_performance(w,r.trajectory));p['operating_mode']=p['operating_mode'].value
        for key in ('cold_heat_per_cycle','hot_heat_per_cycle','gas_work_per_cycle'):
            assert p[key]==pytest.approx(ref['performance'][key],rel=2e-11,abs=1e-11)
        # The original stored periodic rank01 result is also checked by test_research_sixbar.
        for a in (0.,):
            contexts=w.flow_contexts(a,x)
            for i,(o,n) in enumerate(zip((old.heat_in,old.heat_out),(new.heat_in,new.heat_out))):
                lhs=o.thermal_point(350.,x[8+i],context=contexts[i])
                rhs=n.thermal_point(350.,x[8+i],context=contexts[i])
                assert lhs.external_outlet_temperature_k==rhs.external_outlet_temperature_k


@pytest.mark.parametrize('speed,qi,qo,w,mode',[(1,2,-3,-1,'refrigeration'),(-1,3,-2,1,'motor'),(1,2,3,-1,'non_refrigeration'),(1,-2,-3,-1,'non_refrigeration'),(-1,2,-3,-1,'non_refrigeration')])
def test_signed_classification(speed,qi,qo,w,mode):
    assert classify_cycle(speed,qi,qo,w).value==mode

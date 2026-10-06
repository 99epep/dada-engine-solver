"""External boundaries preserve historical algebra and signed cycle accounting."""
from dataclasses import replace, asdict
import json
from pathlib import Path
import numpy as np
import pytest
from dada_solver.exchangers.air_wall import AirWallExchanger
from dada_solver.exchangers.external_stream import ExternalStreamWallMachine
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
        assert 'external_heat_w' in new.rates(t,35000.)
    assert 'external_heat_w' in old.rates(300.,35000.)


@pytest.mark.parametrize('field,value',[('inlet_temperature_k',0.),('mass_flow_kg_s',-.1),('cp_j_kg_k',0.),('wall_conductance_w_k',-1.),('fluid',''),('mass_flow_kg_s',float('nan'))])
def test_stream_validation(field,value):
    data=dict(inlet_temperature_k=300.,mass_flow_kg_s=.1,cp_j_kg_k=4180.,wall_conductance_w_k=100.,fluid='declared liquid')
    data[field]=value
    with pytest.raises(ValueError): ExternalFluidStream(**data)




@pytest.mark.parametrize('speed,qi,qo,w,mode',[(1,2,-3,-1,'refrigeration'),(-1,3,-2,1,'motor'),(1,2,3,-1,'non_refrigeration'),(1,-2,-3,-1,'non_refrigeration'),(-1,2,-3,-1,'non_refrigeration')])
def test_signed_classification(speed,qi,qo,w,mode):
    assert classify_cycle(speed,qi,qo,w).value==mode

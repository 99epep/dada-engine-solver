"""Equivalence to frozen pre-refactor algebra, including physical edge cases."""
from dataclasses import replace
from pathlib import Path
import math

import numpy as np
import pytest

from dada_solver.configuration import load_simulation_configuration
from dada_solver.dynamics import ValveTopology
from dada_solver.factory import build_model, initial_valve_topology
from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.exchangers.hardware import connect_hardware, load_hardware_definition
from dada_solver.state import ThermodynamicState, UniformCharge
from dada_solver.valves import ValveState
from tests import solver_acceleration_reference as reference
from tests.test_dynamics import create_model, LinearInstantKinematics

ROOT = Path(__file__).resolve().parents[1]


def assert_rates_equal(actual, expected):
    np.testing.assert_array_equal(actual.state_derivative, expected.state_derivative)
    for field in ('flows','cold_heat_rate','hot_heat_rate','gas_work_rate',
                  'mass_residual_rate','energy_residual_rate'):
        assert getattr(actual,field) == getattr(expected,field)


@pytest.mark.parametrize('continuous',[False,True])
@pytest.mark.parametrize('stored',[ValveState.OPEN,ValveState.CLOSED])
@pytest.mark.parametrize('delta',[-.2,-1e-12,0.,1e-12,.2])
def test_point_and_balance_match_frozen_reference(ideal_gas, continuous, stored, delta):
    model=replace(create_model(ideal_gas,LinearInstantKinematics(2e-4,5e-4,-2e-6,3e-6)),
                  continuous_ideal_diodes=continuous)
    state=UniformCharge(2e5,300).create_state(ideal_gas,model.volumes(0))
    values=state.as_array(); values[1]*=1+delta; values[7]*=1-delta
    state=ThermodynamicState.from_array(values); topology=ValveTopology(stored,stored)
    expected=reference.evaluate(model,0,state,topology)
    assert_rates_equal(model.evaluate(0,state,topology),expected)
    point=model.instantaneous_point(0,state,topology)
    np.testing.assert_array_equal(point.temperatures,state.temperatures(model.gas))
    np.testing.assert_array_equal(point.pressures,state.pressures(model.gas,model.volumes(0)))
    assert point.volumes == model.volumes(0)
    assert point.volume_rates == model.cylinder_volume_rates(0)
    assert point.topology == reference.effective_topology(model,0,state,topology)
    assert point.flows == expected.flows
    assert_rates_equal(model.assemble_rates(point,expected.cold_heat_rate,expected.hot_heat_rate),expected)


def variable_wrapper():
    config=load_simulation_configuration(ROOT/'examples/motor_demonstrator_original_325c.toml')
    model=build_model(config)
    limits=config.machine_volumes
    model=replace(model,kinematics=FourStageVolumeKinematics(limits.small_cylinder,limits.large_cylinder,
        .25,.45,.75,.1,.8,.2,.9))
    _,bank,hi,ho=load_hardware_definition(
        (ROOT/'examples/motor_hardware_parallel_325c_variable_gas.toml').read_text(),model.gas.heat_capacity_cp)
    wrapper,_=connect_hardware(model,bank,bank,hi,ho,
        heat_in_valve_cda_m2=config.hydraulics.cold_to_large_valve_cda,
        heat_out_valve_cda_m2=config.hydraulics.hot_to_small_valve_cda)
    return wrapper


@pytest.mark.parametrize('temperature',[298.15,598.15])
@pytest.mark.parametrize('delta',[-1e-5,0,1e-5])
def test_wall_rhs_matches_reference_around_breakpoints(temperature,delta):
    wrapper=variable_wrapper()
    for edge in (0.,*wrapper.model.kinematics.breakpoint_angles(),2*math.pi):
        for angle in (edge-1e-10,edge,edge+1e-10):
            gas=UniformCharge(2e5,temperature).create_state(wrapper.model.gas,wrapper.model.volumes(angle))
            values=np.r_[gas.as_array(),wrapper.heat_in.wall_capacity_j_k*450,
                         wrapper.heat_out.wall_capacity_j_k*330,np.zeros(5)]
            values[1]*=1+delta; values[7]*=1-delta
            assert wrapper.flow_contexts(angle,values)==reference.flow_contexts(wrapper,angle,values)
            np.testing.assert_array_equal(wrapper.derivative(angle,values),reference.derivative(wrapper,angle,values))


def test_wall_rhs_solves_hydraulics_once(monkeypatch):
    wrapper=variable_wrapper();cls=type(wrapper.model);original=cls.instantaneous_point;calls=[]
    def counted(self,*args,**kwargs):
        calls.append(1);return original(self,*args,**kwargs)
    monkeypatch.setattr(cls,'instantaneous_point',counted)
    state=UniformCharge(2e5,350).create_state(wrapper.model.gas,wrapper.model.volumes(.3))
    values=np.r_[state.as_array(),wrapper.heat_in.wall_capacity_j_k*450,
                 wrapper.heat_out.wall_capacity_j_k*330,np.zeros(5)]
    wrapper.derivative(.3,values)
    assert len(calls)==1


def test_invalid_transport_state_keeps_exception_semantics():
    wrapper=variable_wrapper()
    state=UniformCharge(2e5,2000).create_state(wrapper.model.gas,wrapper.model.volumes(.3))
    values=np.r_[state.as_array(),wrapper.heat_in.wall_capacity_j_k*450,
                 wrapper.heat_out.wall_capacity_j_k*330,np.zeros(5)]
    with pytest.raises(ValueError) as expected: reference.derivative(wrapper,.3,values)
    with pytest.raises(type(expected.value),match=str(expected.value)):
        wrapper.derivative(.3,values)


def test_tube_validity_reuses_complete_report(monkeypatch):
    from types import SimpleNamespace
    from dada_solver.campaign.evaluator import MachineEvaluator
    import dada_solver.exchangers.gas_diagnostics as diagnostics
    def unexpected(*args):
        raise AssertionError('The completed diagnostic report was recomputed.')
    monkeypatch.setattr(diagnostics,'cycle_microtube_diagnostics',unexpected)
    wrapper=SimpleNamespace(heat_in=SimpleNamespace(requires_flow_context=True),heat_out=None)
    report=dict(passages={'Hi.inlet':dict(hydraulic_upstream_ranges=dict(
        reynolds=dict(maximum=123.),mach=dict(maximum=.02)))})
    assert MachineEvaluator(None)._tube_validity(wrapper,None,None,gas_domains=report)==(123.,.02)


def test_prepared_geometry_is_candidate_owned_and_reporting_is_independent():
    from dataclasses import asdict
    from dada_solver.exchangers.microtube_geometry import MicrotubeBank
    bank=MicrotubeBank(309,.127,.00033,.0001524,.00125,.001)
    original=bank.dimensions();modified=bank.dimensions();modified['tube_flow_area_m2']=0
    assert bank.dimensions()==original
    larger=replace(bank,tube_count=618)
    assert larger.tube_flow_area_m2==2*bank.tube_flow_area_m2
    assert MicrotubeBank(**asdict(bank)).dimensions()==original


def test_prepared_air_constants_refresh_on_replace():
    from dada_solver.exchangers.air_wall import AirWallExchanger
    from dataclasses import asdict
    wall=AirWallExchanger(10,20,100,.05,1005,450)
    changed=replace(wall,air_mass_flow_kg_s=.1,air_wall_conductance_w_k=40)
    assert changed.rates(300,35000)['air_heat_w']==2*wall.rates(300,35000)['air_heat_w']
    assert AirWallExchanger(**asdict(changed)).rates(300,35000)==changed.rates(300,35000)


def test_unsupported_rarefied_state_keeps_domain_error():
    from dada_solver.exchangers.gas_correlations import MicrotubeDomainError
    wrapper=variable_wrapper()
    gas=UniformCharge(100,350).create_state(wrapper.model.gas,wrapper.model.volumes(.3))
    values=np.r_[gas.as_array(),wrapper.heat_in.wall_capacity_j_k*450,
                 wrapper.heat_out.wall_capacity_j_k*330,np.zeros(5)]
    with pytest.raises(MicrotubeDomainError) as expected:
        reference.derivative(wrapper,.3,values)
    with pytest.raises(MicrotubeDomainError) as actual:
        wrapper.derivative(.3,values)
    assert str(actual.value)==str(expected.value)


def test_segment_counters_and_interruption_are_explicit():
    from dada_solver.integration import IntegrationInterrupted
    wrapper=variable_wrapper()
    state=UniformCharge(2e5,350).create_state(wrapper.model.gas,wrapper.model.volumes(0))
    values=np.r_[state.as_array(),wrapper.heat_in.wall_capacity_j_k*450,
                 wrapper.heat_out.wall_capacity_j_k*330]
    records=[]
    def interrupt(_): raise IntegrationInterrupted('Test deadline')
    with pytest.raises(IntegrationInterrupted,match='Test deadline'):
        wrapper.integrate_cycle(values,progress_callback=interrupt,
            progress_interval_seconds=1e-12,statistics_callback=records.append)
    assert records[0]['phase']=='integration_preflight'
    assert records[-1]['status']=='interrupted'
    assert records[-1]['rhs_calls']==0
    assert records[-1]['nfev'] is None
    assert records[-1]['accepted_samples'] is None


def test_complete_segment_counters_preserve_integration(monkeypatch):
    from types import SimpleNamespace
    from dada_solver.exchangers.air_wall import AirWallMotor
    from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor
    wrapper=variable_wrapper()
    # Zero derivative isolates integrator accounting from expensive thermal work.
    monkeypatch.setattr(AirWallMotor,'derivative',lambda self,a,y:np.zeros(15))
    records=[]
    result=solve_periodic_wall_motor(wrapper,np.ones(10),maximum_cycles=1,
        statistics_callback=records.append,measure_rhs_time=True)
    assert result.converged
    segments=[r for r in records if r['phase']=='solver_segment']
    assert len(segments)==4
    assert sum(r['accepted_steps'] for r in segments)==len(result.angles)-4
    assert all(r['nfev']==r['rhs_calls'] for r in segments)
    assert all(r['njev']==r['nlu']==0 for r in segments)
    assert all(r['rhs_elapsed_seconds']>=0 for r in segments)
    assert all(r['maximum_step_angle']<=math.pi/360+1e-14 for r in segments)
    assert [r['ends_at_kinematic_breakpoint'] for r in segments]==[True,True,True,False]
    assert tuple(records)==result.solver_statistics
    np.testing.assert_array_equal(result.last_complete_state,np.ones(10))


def test_candidate_measurement_keeps_failure_record():
    from dada_solver.campaign.evaluator import EvaluationControl
    records=[];control=EvaluationControl(statistics_callback=records.append)
    def fail(): raise ValueError('Test preflight failure')
    with pytest.raises(ValueError,match='Test preflight failure'):
        control.measure('build',fail)
    assert records[0]['phase']=='build'
    assert records[0]['status']=='failed'
    assert records[0]['elapsed_seconds']>=0


def test_interrupt_after_extrapolation_retains_physical_endpoint():
    from types import SimpleNamespace
    from dada_solver.exchangers.wall_cycle import solve_periodic_wall_motor, WallCycleNumericalSettings
    from dada_solver.integration import IntegrationInterrupted
    class Wrapper:
        heat_in=SimpleNamespace(wall_capacity_j_k=1)
        heat_out=SimpleNamespace(wall_capacity_j_k=1)
        calls=0
        endpoint=None
        def integrate_cycle(self,state,**kwargs):
            self.calls+=1
            if self.calls==11:
                assert np.any(state[8:10]!=self.endpoint[8:10])
                raise IntegrationInterrupted('Deadline after initial-guess acceleration')
            end=np.asarray(state).copy();end[8:10]=(end[8:10]+200)/2
            self.endpoint=end.copy()
            trajectory=np.zeros((15,2));trajectory[:10,0]=state;trajectory[:10,1]=end
            return np.array([0.,2*math.pi]),trajectory
    wrapper=Wrapper()
    result=solve_periodic_wall_motor(wrapper,np.r_[np.ones(8),100,100],maximum_cycles=20,
        settings=WallCycleNumericalSettings(accelerate_walls=True))
    assert result.status=='interrupted'
    assert result.history[-1]['wall_initial_guess_extrapolated']
    np.testing.assert_array_equal(result.last_complete_state,wrapper.endpoint)


def test_fast_segments_still_check_deadlines_at_boundaries(monkeypatch):
    from dada_solver.exchangers.air_wall import AirWallMotor
    from dada_solver.integration import IntegrationInterrupted
    wrapper=variable_wrapper();calls=[];records=[]
    monkeypatch.setattr(AirWallMotor,'derivative',lambda self,a,y:np.zeros(15))
    def check(progress):
        calls.append(progress)
        if len(calls)==2: raise IntegrationInterrupted('Boundary deadline')
    with pytest.raises(IntegrationInterrupted,match='Boundary deadline'):
        wrapper.integrate_cycle(np.ones(10),progress_callback=check,
            progress_interval_seconds=1e9,statistics_callback=records.append)
    assert len(calls)==2
    assert calls[0]['current_angle']==0
    assert calls[1]['current_angle']>0
    assert records[-1]['status']=='interrupted'
    assert records[-1]['rhs_calls']==0

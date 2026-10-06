"""Current conservative execution, compiled parity and interruption contracts."""
from dataclasses import replace
import math

import numpy as np
import pytest

from dada_solver.four_stage_kinematics import FourStageVolumeKinematics
from dada_solver.kinematics import ReversedVolumeKinematics
from dada_solver.exchangers.hardware import connect_hardware, HardwareInputs
from dada_solver.exchangers.microtube_geometry import MicrotubeBank
from dada_solver.exchangers.gas_correlations import MicrotubeGasModel
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.state import UniformCharge
from tests.test_dynamics import create_model, LinearInstantKinematics


def variable_wrapper(heat_in='downstream', heat_out='downstream'):
    """Synthetic supported tube/wall machine, independent of application studies."""
    gas = CaloricallyPerfectGas(287.05, 1005., 717.95)
    model = create_model(gas, LinearInstantKinematics(2e-4, 5e-4, 0., 0.))
    limits = model.machine_volumes
    model = replace(
        model, continuous_ideal_diodes=True, study_crank_direction=-1,
        heat_in_valve_placement=heat_in, heat_out_valve_placement=heat_out,
        kinematics=ReversedVolumeKinematics(FourStageVolumeKinematics(
            limits.small_cylinder, limits.large_cylinder,
            .2, .4, .6, .1, .7, .2, .8)),
    )
    bank = MicrotubeBank(256, .2, .001, .0001, .0015, .002)
    inputs = HardwareInputs(
        metal_conductivity_w_m_k=200., metal_density_kg_m3=2700.,
        metal_cp_j_kg_k=900., extra_wall_capacity_j_k=0.,
        gas_conductivity_w_m_k=.026, gas_viscosity_pa_s=1.8e-5,
        gas_nusselt=3.66, air_conductivity_w_m_k=.026,
        air_viscosity_pa_s=1.8e-5, air_nusselt=3.66,
        air_density_kg_m3=1.2, air_cp_j_kg_k=1005.,
        air_mass_flow_kg_s=.02, air_inlet_temperature_k=400.,
        air_poiseuille_number=64., air_minor_loss_coefficient=0.,
        fan_total_efficiency=.5, core_loss_multiplier=1.,
        header_loss_coefficient=0., gas_model=MicrotubeGasModel(),
    )
    wrapper, _ = connect_hardware(
        model, bank, bank, inputs, replace(inputs, air_inlet_temperature_k=300.),
        heat_in_valve_cda_m2=1e-5, heat_out_valve_cda_m2=1e-5,
    )
    return wrapper


@pytest.mark.parametrize('heat_in,heat_out,expected', [
    ('downstream','downstream',(False,False,True,True)),
    ('upstream','downstream',(False,True,True,False)),
    ('downstream','upstream',(True,False,False,True)),
    ('upstream','upstream',(True,True,False,False)),
])
def test_all_valve_placements_match_python_and_compiled_rhs(heat_in, heat_out, expected):
    from dada_solver.wall_backend import WallRHS, WallBackendSettings
    wrapper=variable_wrapper(heat_in,heat_out)
    compiled=WallRHS(wrapper,WallBackendSettings('numba'))
    assert tuple(compiled.implementation.one_way) == expected
    for angle,scale in ((.2,.99999),(.9,1.),(2.1,1.00001)):
        gas=UniformCharge(2e5,350).create_state(wrapper.model.gas,wrapper.model.volumes(angle))
        values=np.r_[gas.as_array(),wrapper.heat_in.wall_capacity_j_k*450,
                     wrapper.heat_out.wall_capacity_j_k*330,np.zeros(5)]
        values[1]*=scale; values[7]/=scale
        np.testing.assert_allclose(compiled(angle,values),wrapper.derivative(angle,values),
                                   rtol=2e-11,atol=1e-11)
    assert compiled.snapshot()['fallback_calls'] == 0


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

"""Conservative table reconstruction, compiled transport energy and strict domains."""
from dataclasses import replace
import numpy as np
import pytest
from dada_solver.fluids import CaloricallyPerfectGas, FluidState
from dada_solver.tabulated_fluid import TabulatedFluid, FluidDomainError, ideal_validation_table
from dada_solver.wall_backend import WallRHS,WallBackendSettings
from dada_solver.state import ThermodynamicState,UniformCharge
from dada_solver.hydraulics import CompressibleOrifice,FlowResult
from tests.test_solver_acceleration import variable_wrapper
from tests.test_wall_backend import values


@pytest.fixture
def fluid(): return ideal_validation_table(CaloricallyPerfectGas(287.,1004.5,717.5))


def test_state_identity_roundtrip_and_immutable_storage(fluid):
    other=TabulatedFluid.from_data(fluid.to_data())
    assert other.content_hash==fluid.content_hash
    for a in (other.properties,other.rho_axis,other.valid_cells):
        with pytest.raises(ValueError): a.flags.writeable=True
    changed=replace(fluid,provenance=fluid.provenance+' different source')
    assert changed.content_hash!=fluid.content_hash
    changed=replace(fluid,limits=(201.,999.,fluid.limits[2],fluid.limits[3]))
    assert changed.content_hash!=fluid.content_hash


@pytest.mark.parametrize('rho',[.1,.25,1.,3.,20.])
def test_interpolation_conservative_exact_ideal(fluid,rho):
    for t in (200.,280.,510.,999.,1000.):
        u=fluid.ideal_reference.heat_capacity_cv*t
        p=fluid.state_from_rho_u(rho,u);reference=fluid.ideal_reference.state_from_rho_u(rho,u)
        for key in ('temperature','pressure','specific_enthalpy','cp','cv','compressibility_factor'):
            assert getattr(p,key)==pytest.approx(getattr(reference,key),rel=8e-15)


@pytest.mark.parametrize('rho,u',[(.099,3e5),(21.,3e5),(1.,1.),(1.,1e7),(float('nan'),3e5)])
def test_table_never_extrapolates(fluid,rho,u):
    with pytest.raises(FluidDomainError): fluid.state_from_rho_u(rho,u)


def test_invalid_single_phase_cell_and_pressure_domain(fluid):
    mask=np.array(fluid.valid_cells);mask[1,1]=False
    masked=replace(fluid,valid_cells=mask)
    with pytest.raises(FluidDomainError): masked.state_from_rho_u(2.,350.*717.5)
    narrow=replace(fluid,limits=(200.,1000.,1e6,2e6))
    with pytest.raises(FluidDomainError): narrow.state_from_rho_u(1.,350.*717.5)


@pytest.mark.parametrize('field', ['pressure','specific_enthalpy','temperature'])
def test_unverified_table_cannot_claim_ideal_hydraulics(fluid,field):
    from dada_solver.tabulated_fluid import PROPERTIES
    p=np.array(fluid.properties);p[:,:,PROPERTIES.index(field)]*=1.01
    with pytest.raises(ValueError,match='declared analytic ideal'): replace(fluid,properties=p)
    unverified=replace(fluid,properties=p,ideal_reference=None)
    with pytest.raises(ValueError,match='hydraulic closure'):
        CompressibleOrifice(1e-4).directed_flow(2e5,1e5,300.,unverified)
    w=variable_wrapper()
    with pytest.raises(ValueError,match='hydraulic closure'): replace(w.model,gas=unverified)


def test_compiled_table_rhs_no_python_properties(monkeypatch):
    pytest.importorskip('numba')
    w=variable_wrapper();fluid=ideal_validation_table(w.model.gas)
    new=replace(w,model=replace(w.model,gas=fluid))
    rhs=WallRHS(new,WallBackendSettings('numba'))
    references=[]
    for a in (0.,.2,1.,3.):
        x=values(w,a);x[1]*=1.00001;x[7]*=.99999
        references.append((a,x,w.derivative(a,x),new.derivative(a,x)))
    def forbidden(*args): raise AssertionError('Python property lookup in compiled RHS')
    monkeypatch.setattr(TabulatedFluid,'state_from_rho_u',forbidden)
    for a,x,exact,python in references:
        actual=rhs(a,x)
        np.testing.assert_allclose(actual,exact,rtol=2e-10,atol=2e-10)
        np.testing.assert_allclose(actual,python,rtol=2e-11,atol=1e-11)
    assert rhs.snapshot()['actual_backend']=='numba_tabulated'
    assert rhs.snapshot()['fallback_calls']==0
    assert rhs.snapshot()['identity']['fluid']['table_sha256']==fluid.content_hash


def test_compiled_out_of_domain_rejects():
    pytest.importorskip('numba')
    w=variable_wrapper();fluid=ideal_validation_table(w.model.gas)
    new=replace(w,model=replace(w.model,gas=fluid));x=values(w);x[1]*=100
    with pytest.raises(FluidDomainError): WallRHS(new,WallBackendSettings('numba'))(0.,x)


def test_nonideal_interface_uses_state_enthalpy():
    # A deterministic interface double, not a claimed thermodynamic dataset.
    class Fluid:
        def state_from_rho_u(self,rho,u):
            return FluidState(rho,u,u/717.5,rho*287.*u/717.5,u+12345.*rho)
    class StateFlow:
        def flow_from_states(self,first,second,fluid,*,one_way):
            return FlowResult(1e-5,False)
    w=variable_wrapper();m=w.model;flow=StateFlow()
    model=replace(m,gas=Fluid(),large_hot_link=flow,small_cold_link=flow,
        hot_small_valve=replace(m.hot_small_valve,flow_model=flow),
        cold_large_valve=replace(m.cold_large_valve,flow_model=flow))
    x=values(w);x[1]*=1.001
    state=ThermodynamicState.from_array(x[:8])
    from dada_solver.factory import initial_valve_topology
    point=model.instantaneous_point(0.,state,initial_valve_topology())
    expected=np.array([s.specific_enthalpy for s in state.fluid_states(model.gas,point.volumes)])
    np.testing.assert_array_equal(point.specific_enthalpies,expected)
    rates=model.assemble_rates(point,0.,0.)
    assert abs(rates.energy_residual_rate)<1e-10
    shifted=replace(point,specific_enthalpies=expected+np.arange(4)*10000.)
    assert not np.allclose(model.assemble_rates(shifted,0.,0.).state_derivative,rates.state_derivative)

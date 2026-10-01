"""Offline CoolProp oracle, species domains and shared compiled transport."""
import json
from pathlib import Path
import numpy as np
import pytest
from dada_solver import numerical_primitives as numeric
from dada_solver.exchangers.gas_transport import DiluteGasTransport, TransportDomainError

DATA=Path(__file__).parent/'data'
REFERENCE=json.loads((DATA/'coolprop8_dilute_transport.json').read_text())


@pytest.mark.parametrize('species,minimum',[('air',100),('nitrogen',200),('argon',200),('helium',50)])
def test_species_domain_and_restrictions(species,minimum):
    tr=DiluteGasTransport(species)
    assert tr.minimum_temperature==minimum
    for method in (tr.viscosity,tr.conductivity,tr.cp):
        assert method(minimum)>0 and method(1000)>0
        for t,category in ((minimum-.1,'below'),(1000.1,'above')):
            with pytest.raises(TransportDomainError) as caught: method(t)
            assert caught.value.category==f'transport_temperature_{category}_domain'
            assert caught.value.diagnostics['species']==species
        with pytest.raises(TransportDomainError): method(float('nan'))
    restricted=DiluteGasTransport(species,minimum+20,900)
    with pytest.raises(TransportDomainError): restricted.cp(minimum)
    for lo,hi in ((minimum-1,1000),(minimum,1001)):
        with pytest.raises(ValueError): DiluteGasTransport(species,lo,hi)


def test_air_observed_boundary():
    assert DiluteGasTransport().viscosity(199.903)>0


@pytest.mark.parametrize('row',REFERENCE['rows'],ids=lambda r:f"{r['species']}-{r['temperature_k']}")
def test_frozen_coolprop8_oracle(row):
    tr=DiluteGasTransport(row['species']);t=row['temperature_k']
    assert tr.viscosity(t)==pytest.approx(row['viscosity_pa_s'],rel=2e-11)
    assert tr.conductivity(t)==pytest.approx(row['conductivity_w_m_k'],rel=2e-11)
    # Approximate 79/21 Shomate air vs pseudo-pure EOS; He keeps rounded DADA R.
    assert tr.cp(t)==pytest.approx(row['cp_j_kg_k'],rel=.008 if row['species']=='air' else 1e-4)
    assert max(row['density_reduction_relative_change'].values())<2e-11


def test_old_range_changes_are_bounded_and_cp_is_unchanged():
    legacy=json.loads((DATA/'legacy_dilute_transport_v1.json').read_text())['rows']
    for row in legacy:
        tr=DiluteGasTransport(row['species']);t=row['temperature_k']
        assert tr.cp(t)==row['cp_j_kg_k']
        for method,key in ((tr.viscosity,'viscosity_pa_s'),(tr.conductivity,'conductivity_w_m_k')):
            tolerance={'air':.041,'helium':.006,'nitrogen':0,'argon':0}[row['species']]
            assert abs(method(t)/row[key]-1)<=tolerance


def test_python_numba_transport_parity():
    numba=pytest.importorskip('numba')
    from dada_solver.wall_backend import compiled_dispatcher
    compiled_dispatcher()  # Register the same shared primitives as the wall kernel.
    functions=[numba.njit(fastmath=False)(f) for f in (numeric.viscosity,numeric.conductivity,numeric.transport_cp)]
    for row in REFERENCE['rows']:
        t=row['temperature_k'];species=numeric.SPECIES.index(row['species'])
        for original,compiled in zip((numeric.viscosity,numeric.conductivity,numeric.transport_cp),functions):
            assert compiled(t,species)==pytest.approx(original(t,species),rel=5e-14)


@pytest.mark.parametrize('temperature',[99.9,100.,150.,199.903])
def test_compiled_rhs_low_air_domain_and_python_fallback(temperature):
    pytest.importorskip('numba')
    from dada_solver.wall_backend import WallRHS,WallBackendSettings
    from dada_solver.state import UniformCharge
    from tests.test_solver_acceleration import variable_wrapper
    wrapper=variable_wrapper()
    state=UniformCharge(2e5,temperature).create_state(wrapper.model.gas,wrapper.model.volumes(.3))
    values=np.r_[state.as_array(),wrapper.heat_in.wall_capacity_j_k*350,
                 wrapper.heat_out.wall_capacity_j_k*350,np.zeros(5)]
    compiled=WallRHS(wrapper,WallBackendSettings('numba'))
    if temperature<100:
        with pytest.raises(TransportDomainError): wrapper.derivative(.3,values)
        with pytest.raises(TransportDomainError): compiled(.3,values)
        assert compiled.implementation.fallbacks==1
    else:
        np.testing.assert_allclose(compiled(.3,values),wrapper.derivative(.3,values),rtol=2e-11,atol=1e-11)

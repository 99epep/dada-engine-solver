from dataclasses import dataclass
from dataclasses import replace

import numpy as np
import pytest

from dada_solver.dynamics import ThermodynamicModel, ValveTopology
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.geometry import CylinderVolumeLimits, MachineVolumes
from dada_solver.heat_transfer import ReservoirHeatTransfer
from dada_solver.hydraulics import CompressibleOrifice
from dada_solver.hydraulics import FlowResult
from dada_solver.state import ThermodynamicState, UniformCharge
from dada_solver.valves import PassiveCheckValve, ValveState
from dada_solver.verification import (
    adiabatic_donor_log_rate_residual,
    closed_adiabatic_log_rate_residual,
)


@dataclass(frozen=True)
class LinearInstantKinematics:
    small_volume: float
    large_volume: float
    small_derivative: float
    large_derivative: float

    def small_cylinder_volume(self, theta: float) -> float:
        return self.small_volume

    def small_cylinder_volume_derivative(self, theta: float) -> float:
        return self.small_derivative

    def large_cylinder_volume(self, theta: float) -> float:
        return self.large_volume

    def large_cylinder_volume_derivative(self, theta: float) -> float:
        return self.large_derivative


@dataclass(frozen=True)
class ZeroBidirectionalFlow:
    def bidirectional_flow(self, *args, **kwargs) -> FlowResult:
        return FlowResult(0.0, False)


def create_model(
    gas: CaloricallyPerfectGas,
    kinematics: LinearInstantKinematics,
) -> ThermodynamicModel:
    volumes = MachineVolumes(
        small_cylinder=CylinderVolumeLimits(1.0e-4, 3.0e-4),
        large_cylinder=CylinderVolumeLimits(2.0e-4, 6.0e-4),
        cold_heat_exchanger=0.8e-4,
        hot_heat_exchanger=1.2e-4,
    )
    return ThermodynamicModel(
        gas=gas,
        machine_volumes=volumes,
        kinematics=kinematics,
        angular_speed=10.0,
        cold_heat_transfer=ReservoirHeatTransfer(5.0, 285.0),
        hot_heat_transfer=ReservoirHeatTransfer(7.0, 315.0),
        large_hot_link=CompressibleOrifice(1.0e-6),
        small_cold_link=CompressibleOrifice(1.0e-6),
        hot_small_valve=PassiveCheckValve(
            CompressibleOrifice(1.0e-6), 100.0, 0.0
        ),
        cold_large_valve=PassiveCheckValve(
            CompressibleOrifice(1.0e-6), 100.0, 0.0
        ),
    )


def test_instantaneous_global_mass_and_energy_balances_are_exact(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    kinematics = LinearInstantKinematics(
        small_volume=2.0e-4,
        large_volume=5.0e-4,
        small_derivative=-2.0e-6,
        large_derivative=3.0e-6,
    )
    model = create_model(ideal_gas, kinematics)
    charge = UniformCharge(pressure=2.0e5, temperature=300.0)
    state = charge.create_state(ideal_gas, model.volumes(theta=0.0))
    perturbed = state.as_array()
    perturbed[1] *= 1.05
    perturbed[7] *= 1.10
    from dada_solver.state import ThermodynamicState

    state = ThermodynamicState.from_array(perturbed)

    rates = model.evaluate(
        theta=0.0,
        state=state,
        topology=ValveTopology(ValveState.OPEN, ValveState.OPEN),
    )

    assert rates.mass_residual_rate == pytest.approx(0.0, abs=1.0e-15)
    assert rates.energy_residual_rate == pytest.approx(0.0, abs=1.0e-10)
    assert np.sum(rates.state_derivative[0::2]) == pytest.approx(0.0, abs=1.0e-15)


def test_simultaneous_valve_opening_is_not_prevented(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    kinematics = LinearInstantKinematics(2.0e-4, 5.0e-4, 0.0, 0.0)
    model = create_model(ideal_gas, kinematics)
    state = UniformCharge(2.0e5, 300.0).create_state(
        ideal_gas, model.volumes(0.0)
    )
    values = state.as_array()
    values[7] *= 1.02
    values[5] *= 1.02
    from dada_solver.state import ThermodynamicState

    state = ThermodynamicState.from_array(values)

    topology = model.valve_transitions(
        0.0,
        state,
        ValveTopology(ValveState.CLOSED, ValveState.CLOSED),
    )

    assert topology.hot_to_small is ValveState.OPEN
    assert topology.cold_to_large is ValveState.OPEN


def test_continuous_diodes_ignore_stored_valve_memory(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model = create_model(
        ideal_gas, LinearInstantKinematics(2.0e-4, 5.0e-4, 0.0, 0.0)
    )
    object.__setattr__(model, "continuous_ideal_diodes", True)
    state = UniformCharge(2.0e5, 300.0).create_state(
        ideal_gas, model.volumes(0.0)
    )
    values = state.as_array()
    values[7] *= 1.02  # H pressure above S pressure.
    from dada_solver.state import ThermodynamicState

    topology = model.effective_topology(
        0.0,
        ThermodynamicState.from_array(values),
        ValveTopology(ValveState.CLOSED, ValveState.OPEN),
    )

    assert topology.hot_to_small is ValveState.OPEN
    assert topology.cold_to_large is ValveState.CLOSED


def test_instantaneous_closed_adiabatic_cylinder_recovers_p_v_gamma_invariant(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    kinematics = LinearInstantKinematics(
        small_volume=2.0e-4,
        large_volume=5.0e-4,
        small_derivative=-3.0e-6,
        large_derivative=2.0e-6,
    )
    model = create_model(ideal_gas, kinematics)
    object.__setattr__(
        model,
        "cold_heat_transfer",
        ReservoirHeatTransfer(0.0, 300.0),
    )
    object.__setattr__(
        model,
        "hot_heat_transfer",
        ReservoirHeatTransfer(0.0, 300.0),
    )
    object.__setattr__(model, "small_cold_link", ZeroBidirectionalFlow())
    object.__setattr__(model, "large_hot_link", ZeroBidirectionalFlow())
    state = UniformCharge(2.0e5, 300.0).create_state(
        ideal_gas, model.volumes(0.0)
    )

    rates = model.evaluate(
        0.0,
        state,
        ValveTopology(ValveState.CLOSED, ValveState.CLOSED),
    )
    pressures = state.pressures(ideal_gas, model.volumes(0.0))
    small_volume_rate, large_volume_rate = model.cylinder_volume_rates(0.0)
    gamma = ideal_gas.heat_capacity_ratio
    for pressure, volume, volume_rate, energy, energy_rate in (
        (
            pressures[0],
            kinematics.small_volume,
            small_volume_rate,
            state.energy_small,
            rates.state_derivative[1],
        ),
        (
            pressures[1],
            kinematics.large_volume,
            large_volume_rate,
            state.energy_large,
            rates.state_derivative[3],
        ),
    ):
        pressure_log_rate = energy_rate / energy - volume_rate / volume
        invariant_log_rate = closed_adiabatic_log_rate_residual(
            pressure, pressure * pressure_log_rate, volume, volume_rate, gamma
        )
        assert invariant_log_rate == pytest.approx(0.0, abs=1.0e-12)


def test_adiabatic_outflow_only_donor_matches_appendix_a5(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    kinematics = LinearInstantKinematics(2.0e-4, 5.0e-4, 0.0, -2.0e-6)
    model = create_model(ideal_gas, kinematics)
    object.__setattr__(
        model, "cold_heat_transfer", ReservoirHeatTransfer(0.0, 300.0)
    )
    object.__setattr__(
        model, "hot_heat_transfer", ReservoirHeatTransfer(0.0, 300.0)
    )
    state = UniformCharge(2.0e5, 300.0).create_state(
        ideal_gas, model.volumes(0.0)
    )
    values = state.as_array()
    values[3] *= 1.20
    from dada_solver.state import ThermodynamicState

    state = ThermodynamicState.from_array(values)
    rates = model.evaluate(
        0.0,
        state,
        ValveTopology(ValveState.CLOSED, ValveState.CLOSED),
    )
    assert rates.flows.large_to_hot > 0.0
    mass_rate = rates.state_derivative[2]
    energy_rate = rates.state_derivative[3]
    temperature = state.temperatures(ideal_gas)[1]
    temperature_rate = (
        energy_rate - ideal_gas.heat_capacity_cv * temperature * mass_rate
    ) / (ideal_gas.heat_capacity_cv * state.mass_large)
    _, large_volume_rate = model.cylinder_volume_rates(0.0)

    residual = adiabatic_donor_log_rate_residual(
        state.mass_large,
        mass_rate,
        temperature,
        temperature_rate,
        kinematics.large_volume,
        large_volume_rate,
        ideal_gas.heat_capacity_ratio,
    )

    assert residual == pytest.approx(0.0, abs=1.0e-12)


@pytest.mark.parametrize('heat_in,heat_out', [
    ('downstream','downstream'), ('upstream','downstream'),
    ('downstream','upstream'), ('upstream','upstream'),
])
def test_valve_half_links_block_reverse_flow_and_free_halves_accept_it(
    ideal_gas, heat_in, heat_out,
) -> None:
    model=replace(create_model(ideal_gas,
        LinearInstantKinematics(2e-4,5e-4,0.,0.)),continuous_ideal_diodes=True,
        heat_in_valve_placement=heat_in,heat_out_valve_placement=heat_out)
    topology=ValveTopology(ValveState.CLOSED,ValveState.CLOSED)
    pairs=((1,3,heat_out=='upstream'),(0,2,heat_in=='upstream'),
           (3,0,heat_out=='downstream'),(2,1,heat_in=='downstream'))
    names=('large_to_hot','small_to_cold','hot_to_small','cold_to_large')
    v=model.volumes(0);volumes=np.asarray((v.small_cylinder,v.large_cylinder,
        v.cold_heat_exchanger,v.hot_heat_exchanger))
    for name,(source,destination,one_way) in zip(names,pairs):
        pressures=np.full(4,2e5);pressures[source]=1.9e5;pressures[destination]=2.1e5
        temperature=350.;masses=pressures*volumes/(ideal_gas.gas_constant*temperature)
        energies=masses*ideal_gas.heat_capacity_cv*temperature
        values=np.empty(8);values[0::2]=masses;values[1::2]=energies
        rates=model.evaluate(0,ThermodynamicState.from_array(values),topology)
        flow=getattr(rates.flows,name)
        assert (flow==0) if one_way else (flow<0)
        assert rates.mass_residual_rate == pytest.approx(0,abs=1e-15)
        assert rates.energy_residual_rate == pytest.approx(0,abs=1e-10)

import numpy as np
import pytest

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.geometry import InstantaneousVolumes
from dada_solver.state import STATE_NAMES, ThermodynamicState, UniformCharge


def test_uniform_charge_reconstructs_uniform_pressure_and_temperature(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    volumes = InstantaneousVolumes(
        small_cylinder=1.2e-4,
        large_cylinder=5.4e-4,
        cold_heat_exchanger=0.8e-4,
        hot_heat_exchanger=1.1e-4,
    )
    charge = UniformCharge(pressure=4.0e5, temperature=298.0)

    state = charge.create_state(ideal_gas, volumes)

    np.testing.assert_allclose(state.temperatures(ideal_gas), 298.0)
    np.testing.assert_allclose(state.pressures(ideal_gas, volumes), 4.0e5)
    assert state.total_mass == pytest.approx(charge.total_mass(ideal_gas, volumes))


def test_state_array_uses_published_ordering() -> None:
    values = np.arange(1.0, 9.0)
    state = ThermodynamicState.from_array(values)

    assert STATE_NAMES == ("m_S", "U_S", "m_L", "U_L", "m_C", "U_C", "m_H", "U_H")
    np.testing.assert_array_equal(state.as_array(), values)


def test_invalid_state_is_rejected() -> None:
    values = np.ones(8)
    values[4] = 0.0

    with pytest.raises(ValueError, match="strictly positive"):
        ThermodynamicState.from_array(values)

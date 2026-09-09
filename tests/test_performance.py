import pytest

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.performance import OperatingMode, calculate_cycle_performance
from tests.test_periodic import create_static_cycle


def test_static_cycle_has_consistent_zero_energy_balance(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model, state, cycle = create_static_cycle(ideal_gas)

    performance = calculate_cycle_performance(cycle, state, model.angular_speed)

    assert performance.cold_heat_per_cycle == pytest.approx(0.0)
    assert performance.hot_heat_per_cycle == pytest.approx(0.0)
    assert performance.gas_work_per_cycle == pytest.approx(0.0)
    assert performance.cooling_cop is None
    assert performance.operating_mode is OperatingMode.NON_REFRIGERATION
    assert performance.conservation.absolute_mass_residual == pytest.approx(0.0)
    assert performance.conservation.absolute_energy_residual == pytest.approx(0.0)


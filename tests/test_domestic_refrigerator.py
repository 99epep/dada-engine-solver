import pytest

from dada_solver.domestic_refrigerator import (
    AnnualEnergyBudget,
    RefrigeratorTemperaturePoint,
    SteadyCabinetLoad,
    annualize_cycling_appliance,
    assess_appliance_power,
)


def test_refrigerator_temperature_point_reports_carnot_bound() -> None:
    point = RefrigeratorTemperaturePoint("fresh_food_32C", 305.15, 277.15)

    assert point.temperature_lift == pytest.approx(28.0)
    assert point.carnot_cooling_cop == pytest.approx(277.15 / 28.0)


def test_appliance_boundary_keeps_auxiliaries_separate() -> None:
    point = RefrigeratorTemperaturePoint("fresh_food_32C", 305.15, 277.15)
    cabinet = SteadyCabinetLoad(thermal_conductance=1.0, internal_heat_load=2.0)

    result = assess_appliance_power(point, cabinet, 2.5, 0.8, 5.0)

    assert result.cabinet_cooling_power == pytest.approx(30.0)
    assert result.thermodynamic_input_power == pytest.approx(12.0)
    assert result.shaft_input_power == pytest.approx(15.0)
    assert result.total_appliance_input_power == pytest.approx(20.0)


def test_annual_energy_budget_converts_to_average_and_on_cycle_power() -> None:
    budget = AnnualEnergyBudget(annual_energy_kwh=200.0)

    assert budget.average_input_power == pytest.approx(200000.0 / 8760.0)
    assert budget.maximum_on_input_power(0.5, off_cycle_power=1.0) == pytest.approx(
        (200000.0 / 8760.0 - 0.5) / 0.5
    )


def test_cycling_assessment_preserves_on_and_off_power_boundaries() -> None:
    point = RefrigeratorTemperaturePoint("frozen_32C", 305.15, 255.15)
    on_cycle = assess_appliance_power(
        point,
        SteadyCabinetLoad(thermal_conductance=1.0),
        machine_cooling_cop=2.0,
        transmission_efficiency=0.8,
        auxiliary_power=3.0,
    )

    result = annualize_cycling_appliance(
        on_cycle,
        duty_fraction=0.4,
        off_cycle_input_power=1.0,
    )

    assert on_cycle.total_appliance_input_power == pytest.approx(34.25)
    assert result.average_input_power == pytest.approx(14.3)
    assert result.annual_energy_kwh == pytest.approx(125.268)


@pytest.mark.parametrize("duty_fraction", [0.0, -0.1, 1.1])
def test_annual_boundary_rejects_invalid_duty_fraction(duty_fraction: float) -> None:
    budget = AnnualEnergyBudget(200.0)

    with pytest.raises(ValueError, match="Duty fraction"):
        budget.maximum_on_input_power(duty_fraction)

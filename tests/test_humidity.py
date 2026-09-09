import math

import pytest

from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.humidity import (
    HumidityScreeningConfiguration,
    MoistureScreeningVerdict,
    assess_moisture_phase_change_risk,
    saturation_vapor_pressure,
)
from tests.test_periodic import create_static_cycle


def test_saturation_pressure_reference_values_near_freezing() -> None:
    assert saturation_vapor_pressure(273.15, over_ice=False) == pytest.approx(
        611.21, rel=2.0e-4
    )
    assert saturation_vapor_pressure(263.15, over_ice=True) == pytest.approx(
        259.9, rel=2.0e-3
    )


def test_static_cycle_below_dew_point_reports_phase_change_risk(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model, _, cycle = create_static_cycle(ideal_gas)
    report = assess_moisture_phase_change_risk(
        cycle,
        model,
        HumidityScreeningConfiguration(relative_humidity=0.5),
        charge_temperature=320.0,
        charge_pressure=2.0e5,
    )

    # The controlled cycle is at 200 kPa and 300 K. Air charged at 320 K and
    # 50% RH has an unchanged water mole fraction above saturation at 300 K.
    assert report.verdict is MoistureScreeningVerdict.CONDENSATION_OR_FROST_RISK
    assert report.initial_water_mole_fraction is not None
    assert math.isfinite(report.control_volumes["C"].maximum_saturation_ratio)
    assert report.control_volumes["C"].first_saturated_angle_degrees is not None
    assert report.control_volumes["C"].first_predicted_phase == "liquid_condensation"


def test_missing_humidity_is_reported_as_unavailable(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model, _, cycle = create_static_cycle(ideal_gas)
    report = assess_moisture_phase_change_risk(
        cycle, model, None, charge_temperature=298.15, charge_pressure=2.0e5
    )

    assert report.verdict is MoistureScreeningVerdict.UNAVAILABLE
    assert report.initial_water_mole_fraction is None


def test_relative_humidity_range_is_checked() -> None:
    with pytest.raises(ValueError, match="Relative humidity"):
        HumidityScreeningConfiguration(relative_humidity=1.01)


def test_out_of_range_cycle_temperature_makes_clear_screening_unavailable(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model, _, cycle = create_static_cycle(ideal_gas)
    cycle.states[1::2, :] *= 2.0
    report = assess_moisture_phase_change_risk(
        cycle,
        model,
        HumidityScreeningConfiguration(relative_humidity=0.0),
        charge_temperature=298.15,
        charge_pressure=2.0e5,
    )

    assert report.verdict is MoistureScreeningVerdict.UNAVAILABLE
    assert report.unavailable_reason is not None

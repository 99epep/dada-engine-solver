import pytest
from dataclasses import replace

from dada_solver.configuration import ValidityThresholds
from dada_solver.fluids import CaloricallyPerfectGas
from dada_solver.validity import ValidityVerdict, assess_cycle_validity
from tests.test_periodic import create_static_cycle


def test_missing_geometric_area_makes_mach_unavailable_and_verdict_indeterminate(
    ideal_gas: CaloricallyPerfectGas,
) -> None:
    model, _, cycle = create_static_cycle(ideal_gas)
    thresholds = ValidityThresholds(
        maximum_pressure_equalization_error=0.01,
        maximum_mach_number=0.2,
        maximum_isothermality_error=0.01,
        maximum_compressibility_deviation=0.01,
        maximum_cp_variation=0.01,
    )

    report = assess_cycle_validity(cycle, model, thresholds)

    assert report.verdict is ValidityVerdict.INDETERMINATE
    assert report.maximum_mach_number is None
    assert report.mach_number_status == "unavailable"
    assert report.maximum_pressure_equalization_error == pytest.approx(0.0)
    assert report.maximum_compressibility_deviation == pytest.approx(0.0)


def test_temperature_excursion_is_a_diagnostic_not_a_validity_veto(ideal_gas):
    model, _, cycle = create_static_cycle(ideal_gas)
    values = cycle.states.copy()
    values[7, -1] *= 1.2
    cycle = replace(cycle, states=values)
    thresholds = ValidityThresholds(1.0, 0.2, 0.01, 0.01, 0.01)
    report = assess_cycle_validity(cycle, model, thresholds)
    assert report.hot_isothermality_error > thresholds.maximum_isothermality_error
    assert "isothermality" not in report.failed_criteria
    assert report.verdict is ValidityVerdict.INDETERMINATE

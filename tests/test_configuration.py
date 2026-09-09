import pytest

from dada_solver.configuration import (
    HydraulicNetworkConfiguration,
    ValidityThresholds,
)


def test_all_first_level_effective_areas_are_required() -> None:
    configuration = HydraulicNetworkConfiguration(
        large_to_hot_cda=1.0e-6,
        small_to_cold_cda=2.0e-6,
        hot_to_small_valve_cda=3.0e-6,
        cold_to_large_valve_cda=4.0e-6,
    )

    assert configuration.cold_to_large_valve_cda == pytest.approx(4.0e-6)


def test_non_positive_effective_area_is_rejected() -> None:
    with pytest.raises(ValueError, match="positive"):
        HydraulicNetworkConfiguration(
            large_to_hot_cda=0.0,
            small_to_cold_cda=2.0e-6,
            hot_to_small_valve_cda=3.0e-6,
            cold_to_large_valve_cda=4.0e-6,
        )


def test_validity_thresholds_are_explicit_and_positive() -> None:
    thresholds = ValidityThresholds(
        maximum_pressure_equalization_error=0.05,
        maximum_mach_number=0.2,
        maximum_isothermality_error=0.03,
        maximum_compressibility_deviation=0.01,
        maximum_cp_variation=0.02,
    )

    assert thresholds.maximum_mach_number == pytest.approx(0.2)

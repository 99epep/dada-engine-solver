import numpy as np

from dada_solver.exchangers.hydraulic_validity import (
    HydraulicValidityThresholds,
    assess_hydraulic_quasi_steady_validity,
)
from dada_solver.exchangers.models import ParallelRectangularChannels


def test_inertial_diagnostic_requests_dynamic_model_when_threshold_is_exceeded(
    ideal_gas,
) -> None:
    geometry = ParallelRectangularChannels(10, 0.01, 0.001, 0.2)
    phase = np.linspace(0.0, 2.0 * np.pi, 101, endpoint=False)
    flow = 0.01 * np.sin(phase)

    result = assess_hydraulic_quasi_steady_validity(
        geometry=geometry,
        mass_flow_history=flow,
        period=0.1,
        reference_pressure=1.0e5,
        reference_temperature=300.0,
        maximum_pressure_drop=100.0,
        gas=ideal_gas,
        thresholds=HydraulicValidityThresholds(1.0, 1.0, 1.0, 1.0e-6),
    )

    assert not result.quasi_steady_valid
    assert result.inertia_model_recommended
    assert "maximum_inertial_pressure_fraction" in result.failed_criteria


def test_slow_low_flow_case_satisfies_configured_quasi_steady_limits(ideal_gas) -> None:
    geometry = ParallelRectangularChannels(100, 0.02, 0.002, 0.05)
    flow = np.zeros(20)

    result = assess_hydraulic_quasi_steady_validity(
        geometry,
        flow,
        10.0,
        1.0e5,
        300.0,
        10.0,
        ideal_gas,
        HydraulicValidityThresholds(0.3, 0.1, 0.01, 0.01),
    )

    assert result.quasi_steady_valid
    assert not result.inertia_model_recommended

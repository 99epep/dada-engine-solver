from dada_solver.exchangers.hydraulic_comparison import compare_hydraulic_models
from dada_solver.hydraulics import (
    CompressibleOrifice,
    QuasiSteadyCompressibleDuct,
    RectangularDuct,
)


def test_comparison_preserves_both_models_without_calibration(ideal_gas) -> None:
    reference = CompressibleOrifice(1.0e-5)
    geometric = QuasiSteadyCompressibleDuct(
        RectangularDuct(1.0e-4, 0.002, 0.1, 0.2, 1.8e-5)
    )

    result = compare_hydraulic_models(
        reference,
        geometric,
        2.0e5,
        (2.0e5, 1.8e5, 1.0e5),
        300.0,
        ideal_gas,
    )

    assert len(result) == 3
    assert result[0].relative_flow_difference is None
    assert result[1].reference_mass_flow > 0.0
    assert result[1].geometric_mass_flow > 0.0
    assert result[1].relative_flow_difference is not None

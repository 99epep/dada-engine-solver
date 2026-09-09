import pytest

from dada_solver.verification import (
    adiabatic_donor_log_rate_residual,
    closed_adiabatic_invariant,
)


def test_closed_adiabatic_invariant_matches_pressure_volume_law() -> None:
    gamma = 1.4
    initial_pressure = 2.0e5
    initial_volume = 4.0e-4
    final_volume = 2.5e-4
    final_pressure = initial_pressure * (initial_volume / final_volume) ** gamma

    assert closed_adiabatic_invariant(
        final_pressure, final_volume, gamma
    ) == pytest.approx(
        closed_adiabatic_invariant(initial_pressure, initial_volume, gamma)
    )


def test_donor_residual_rejects_non_physical_mass() -> None:
    with pytest.raises(ValueError, match="mass"):
        adiabatic_donor_log_rate_residual(
            mass=0.0,
            mass_rate=-0.1,
            temperature=300.0,
            temperature_rate=0.0,
            volume=1.0e-3,
            volume_rate=0.0,
            gamma=1.4,
        )

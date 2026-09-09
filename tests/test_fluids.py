import pytest

from dada_solver.fluids import CaloricallyPerfectGas


def test_caloric_properties_are_consistent(ideal_gas: CaloricallyPerfectGas) -> None:
    assert ideal_gas.heat_capacity_ratio == pytest.approx(1004.5 / 717.5)
    assert ideal_gas.compressibility_factor(pressure=2.0e5, temperature=300.0) == 1.0


def test_inconsistent_caloric_properties_are_rejected() -> None:
    with pytest.raises(ValueError, match="R = Cp - Cv"):
        CaloricallyPerfectGas(
            gas_constant=300.0,
            heat_capacity_cp=1004.5,
            heat_capacity_cv=717.5,
        )


def test_ideal_gas_state_reconstruction(ideal_gas: CaloricallyPerfectGas) -> None:
    pressure = 3.2e5
    temperature = 315.0
    volume = 2.5e-4
    mass = ideal_gas.mass(pressure, temperature, volume)
    internal_energy = ideal_gas.internal_energy(mass, temperature)

    assert ideal_gas.temperature(mass, internal_energy) == pytest.approx(temperature)
    assert ideal_gas.pressure(mass, internal_energy, volume) == pytest.approx(pressure)
    assert ideal_gas.enthalpy(temperature) == pytest.approx(
        ideal_gas.heat_capacity_cp * temperature
    )


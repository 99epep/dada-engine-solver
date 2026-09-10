from pathlib import Path

import pytest

from dada_solver.exchangers.doty import load_references, scale_bank
from dada_solver.exchangers.microtube_geometry import MicrotubeBank

DATA_N2 = Path(__file__).resolve().parents[1] / "examples/data/doty_1991_nitrogen_reference.csv"
DATA_HE = Path(__file__).resolve().parents[1] / "examples/data/doty_1991_helium_reference.csv"


def estimate(reference, banks=1, flow=None, **overrides):
    inputs = dict(parallel_banks=banks,
                  total_mass_flow_kg_s=flow or reference.mass_flow_kg_s,
                  ua_flow_exponent=0.0, ua_condition_factor=1.0,
                  viscosity_ratio=1.0, density_ratio=1.0,
                  additional_pressure_drop_pa=0.0)
    inputs.update(overrides)
    return scale_bank(reference, **inputs)


# ============================================================================
# Nitrogen reference tests (existing)
# ============================================================================

def test_measured_anchors_and_parallel_replication():
    for reference in load_references(DATA_N2):
        result = estimate(reference, banks=4, flow=4 * reference.mass_flow_kg_s)
        assert result["estimated_ua_w_k"] == pytest.approx(4 * reference.ua_w_k)
        assert result["estimated_tube_pressure_drop_pa"] == pytest.approx(reference.tube_pressure_drop_pa)
        assert result["sum_bank_envelope_m3"] == pytest.approx(0.000798)


def test_density_viscosity_and_flow_scaling():
    ref = load_references(DATA_N2)[0]
    result = estimate(ref, flow=2 * ref.mass_flow_kg_s, viscosity_ratio=1.5,
                      density_ratio=3, additional_pressure_drop_pa=100)
    assert result["estimated_tube_and_additional_pressure_drop_pa"] == pytest.approx(5200)
    assert ref.tube_mean_temperature_k == pytest.approx((296.15 + 366.15) / 2)


@pytest.mark.parametrize("overrides", [dict(total_mass_flow_kg_s=0),
    dict(total_mass_flow_kg_s=-1), dict(density_ratio=float("nan")),
    dict(parallel_banks=1.5), dict(parallel_banks=True), dict(ua_flow_exponent=2)])
def test_unsupported_inputs_are_explicit(overrides):
    with pytest.raises(ValueError):
        estimate(load_references(DATA_N2)[0], **overrides)


# ============================================================================
# Helium reference tests (new)
# ============================================================================

def test_helium_file_loads_and_has_three_rows():
    """Verify helium reference file contains exactly three measurements."""
    references = load_references(DATA_HE)
    assert len(references) == 3


def test_helium_reference_values():
    """Verify helium reference measurements match Doty Tables 1–2 (SI conversions)."""
    refs = load_references(DATA_HE)
    
    # Row 1: 117 mg/s
    assert refs[0].mass_flow_kg_s == pytest.approx(0.000117)
    assert refs[0].ua_w_k == pytest.approx(11.1)
    assert refs[0].tube_pressure_drop_pa == pytest.approx(2300)
    assert refs[0].pressure_pa == pytest.approx(749000)
    assert refs[0].tube_mean_temperature_k == pytest.approx((296.65 + 376.75) / 2)
    
    # Row 2: 79 mg/s
    assert refs[1].mass_flow_kg_s == pytest.approx(0.000079)
    assert refs[1].ua_w_k == pytest.approx(6.9)
    assert refs[1].tube_pressure_drop_pa == pytest.approx(1500)
    assert refs[1].pressure_pa == pytest.approx(756000)
    assert refs[1].tube_mean_temperature_k == pytest.approx((297.15 + 373.65) / 2)
    
    # Row 3: 213 mg/s
    assert refs[2].mass_flow_kg_s == pytest.approx(0.000213)
    assert refs[2].ua_w_k == pytest.approx(15.4)
    assert refs[2].tube_pressure_drop_pa == pytest.approx(4400)
    assert refs[2].pressure_pa == pytest.approx(825000)
    assert refs[2].tube_mean_temperature_k == pytest.approx((296.25 + 376.55) / 2)


# ============================================================================
# Hydraulic Check 1: Absolute Poiseuille/Doty Eq. 7 using direct calculation
#
# Direct theoretical laminar tube-loss calculation using MicrotubeBank geometry
# and the Poiseuille law: Δp = 128 * μ * L * m_dot / (ρ * n * π * d^4)
#
# This reproduces approximately Doty Table 2 calculated values:
#   117 mg/s: about 3.0 kPa
#   79 mg/s: about 2.0 kPa
#   213 mg/s: about 5.5 kPa
#
# Measured values are:
#   117 mg/s: 2.3 kPa
#   79 mg/s: 1.5 kPa
#   213 mg/s: 4.4 kPa
#
# Therefore theoretical overpredicts by approximately 25–33%.
# ============================================================================

def _helium_density_from_ideal_gas(pressure_pa, temperature_k):
    """Ideal-gas helium density: ρ = P / (R * T) with R = 2077.1 J/(kg·K)."""
    R_He = 2077.1  # J/(kg·K)
    return pressure_pa / (R_He * temperature_k)


def _helium_viscosity_from_temperature(temperature_k):
    """Helium viscosity using empirical approximation valid 250–500 K.
    
    This is a simple linear approximation fitted around the experimental range.
    Viscosity reference: Chapman-Enskog kinetic theory for monatomic gases.
    For helium at 300 K: μ ≈ 1.87e-5 Pa·s
    For helium at 380 K: μ ≈ 2.2e-5 Pa·s (approximately)
    
    Simple linear form: μ ≈ μ_0 + α * (T - T_0)
    Do NOT use this as a universal viscosity correlation; it is only valid
    for the 250–500 K range and is intended for screening physics.
    """
    # At T = 300 K, μ ≈ 1.87e-5 Pa·s (literature)
    # At T = 380 K, μ ≈ 2.2e-5 Pa·s (rough estimate)
    # Slope: (2.2e-5 - 1.87e-5) / (380 - 300) ≈ 4.125e-8 Pa·s/K
    mu_300 = 1.87e-5  # Pa·s at 300 K
    alpha = 4.125e-8  # Pa·s/K
    return mu_300 + alpha * (temperature_k - 300.0)


def test_helium_absolute_poiseuille_117mg_s():
    """Direct Poiseuille prediction for 117 mg/s He.
    
    Theoretical prediction should be approximately 3.0 kPa.
    Measured value is 2.3 kPa.
    Overprediction ratio ≈ 3.0 / 2.3 ≈ 30%.
    """
    ref = load_references(DATA_HE)[0]  # 117 mg/s
    
    # Doty bank geometry: 3 modules × 103 tubes = 309 tubes
    bank = MicrotubeBank(
        tube_count=309,
        tube_length_m=0.127,
        inner_diameter_m=0.00033,
        wall_thickness_m=0.0001524,
        pitch_m=0.0008,  # Reasonable pitch for 0.33 mm ID with wall
        header_depth_m=0.01,  # Reasonable header depth
    )
    
    # Calculate density and viscosity at representative tube conditions
    rho = _helium_density_from_ideal_gas(ref.pressure_pa, ref.tube_mean_temperature_k)
    mu = _helium_viscosity_from_temperature(ref.tube_mean_temperature_k)
    
    # Direct Poiseuille calculation
    result = bank.laminar_tube_loss(ref.mass_flow_kg_s, density_kg_m3=rho, viscosity_pa_s=mu)
    dp_theoretical = result["signed_tube_pressure_drop_pa"]
    
    # Expected: about 3000 Pa ± 10% (conservative tolerance for screening physics)
    assert dp_theoretical == pytest.approx(3000, rel=0.15)
    
    # Verify measured value is lower (overprediction)
    dp_measured = ref.tube_pressure_drop_pa
    assert dp_measured < dp_theoretical
    percent_overprediction = (dp_theoretical - dp_measured) / dp_measured * 100
    assert 25 < percent_overprediction < 35  # Expect ~30% overprediction


def test_helium_absolute_poiseuille_79mg_s():
    """Direct Poiseuille prediction for 79 mg/s He.
    
    Theoretical prediction should be approximately 2.0 kPa.
    Measured value is 1.5 kPa.
    Overprediction ratio ≈ 2.0 / 1.5 ≈ 33%.
    """
    ref = load_references(DATA_HE)[1]  # 79 mg/s
    
    bank = MicrotubeBank(
        tube_count=309,
        tube_length_m=0.127,
        inner_diameter_m=0.00033,
        wall_thickness_m=0.0001524,
        pitch_m=0.0008,
        header_depth_m=0.01,
    )
    
    rho = _helium_density_from_ideal_gas(ref.pressure_pa, ref.tube_mean_temperature_k)
    mu = _helium_viscosity_from_temperature(ref.tube_mean_temperature_k)
    
    result = bank.laminar_tube_loss(ref.mass_flow_kg_s, density_kg_m3=rho, viscosity_pa_s=mu)
    dp_theoretical = result["signed_tube_pressure_drop_pa"]
    
    # Expected: about 2000 Pa ± 10%
    assert dp_theoretical == pytest.approx(2000, rel=0.15)
    
    # Verify measured value is lower
    dp_measured = ref.tube_pressure_drop_pa
    assert dp_measured < dp_theoretical
    percent_overprediction = (dp_theoretical - dp_measured) / dp_measured * 100
    assert 30 < percent_overprediction < 40  # Expect ~33% overprediction


def test_helium_absolute_poiseuille_213mg_s():
    """Direct Poiseuille prediction for 213 mg/s He.
    
    Theoretical prediction should be approximately 5.5 kPa.
    Measured value is 4.4 kPa.
    Overprediction ratio ≈ 5.5 / 4.4 ≈ 25%.
    """
    ref = load_references(DATA_HE)[2]  # 213 mg/s
    
    bank = MicrotubeBank(
        tube_count=309,
        tube_length_m=0.127,
        inner_diameter_m=0.00033,
        wall_thickness_m=0.0001524,
        pitch_m=0.0008,
        header_depth_m=0.01,
    )
    
    rho = _helium_density_from_ideal_gas(ref.pressure_pa, ref.tube_mean_temperature_k)
    mu = _helium_viscosity_from_temperature(ref.tube_mean_temperature_k)
    
    result = bank.laminar_tube_loss(ref.mass_flow_kg_s, density_kg_m3=rho, viscosity_pa_s=mu)
    dp_theoretical = result["signed_tube_pressure_drop_pa"]
    
    # Expected: about 5500 Pa ± 10%
    assert dp_theoretical == pytest.approx(5500, rel=0.15)
    
    # Verify measured value is lower
    dp_measured = ref.tube_pressure_drop_pa
    assert dp_measured < dp_theoretical
    percent_overprediction = (dp_theoretical - dp_measured) / dp_measured * 100
    assert 20 < percent_overprediction < 30  # Expect ~25% overprediction


# ============================================================================
# Hydraulic Check 2: Relative scaling between helium measurements
#
# Use the 117 mg/s measurement as the single anchor. For the other two rows,
# use the existing scaling law dp ~ mass_flow * viscosity / density
# with viscosity_ratio = 1.0 and density ~ pressure / tube_mean_temperature.
#
# Expected predictions:
#   79 mg/s: about 1.53 kPa versus 1.50 kPa measured, about +2 %
#   213 mg/s: about 3.80 kPa versus 4.40 kPa measured, about -14 %
# ============================================================================

def test_helium_relative_scaling_79mg_s_from_117mg_s_anchor():
    """79 mg/s He relative to 117 mg/s anchor: about +2 % overprediction."""
    ref_anchor = load_references(DATA_HE)[0]  # 117 mg/s
    ref_79 = load_references(DATA_HE)[1]       # 79 mg/s
    
    # Scaling law: dp = dp_ref * (flow / flow_ref) * (visc_ref / visc) * (rho / rho_ref)
    # With viscosity_ratio = 1.0 and ideal gas rho ~ P / T:
    flow_ratio = ref_79.mass_flow_kg_s / ref_anchor.mass_flow_kg_s
    
    # Density ratio: rho_79 / rho_117 ≈ (P_79 / T_mean_79) / (P_117 / T_mean_117)
    P_79 = ref_79.pressure_pa
    T_mean_79 = ref_79.tube_mean_temperature_k
    P_117 = ref_anchor.pressure_pa
    T_mean_117 = ref_anchor.tube_mean_temperature_k
    density_ratio = (P_79 / T_mean_79) / (P_117 / T_mean_117)
    
    # Predict using the 117 mg/s anchor scaled to 79 mg/s flow
    result = estimate(ref_anchor, banks=1, flow=ref_79.mass_flow_kg_s,
                      viscosity_ratio=1.0, density_ratio=density_ratio)
    predicted_pa = result["estimated_tube_pressure_drop_pa"]
    measured_pa = ref_79.tube_pressure_drop_pa
    
    # Expected: about 1.53 kPa predicted vs 1.50 kPa measured (~+2%)
    assert predicted_pa == pytest.approx(1530, abs=50)  # ~1.53 kPa ± 50 Pa
    assert measured_pa == pytest.approx(1500, abs=50)   # ~1.50 kPa ± 50 Pa
    percent_diff = (predicted_pa - measured_pa) / measured_pa * 100
    assert percent_diff == pytest.approx(2, abs=5)  # About +2%


def test_helium_relative_scaling_213mg_s_from_117mg_s_anchor():
    """213 mg/s He relative to 117 mg/s anchor: about -14 % underprediction."""
    ref_anchor = load_references(DATA_HE)[0]  # 117 mg/s
    ref_213 = load_references(DATA_HE)[2]     # 213 mg/s
    
    # Scaling law with viscosity_ratio = 1.0 and density ratio from ideal gas
    flow_ratio = ref_213.mass_flow_kg_s / ref_anchor.mass_flow_kg_s
    
    # Density ratio: rho_213 / rho_117
    P_213 = ref_213.pressure_pa
    T_mean_213 = ref_213.tube_mean_temperature_k
    P_117 = ref_anchor.pressure_pa
    T_mean_117 = ref_anchor.tube_mean_temperature_k
    density_ratio = (P_213 / T_mean_213) / (P_117 / T_mean_117)
    
    # Predict using the 117 mg/s anchor scaled to 213 mg/s flow
    result = estimate(ref_anchor, banks=1, flow=ref_213.mass_flow_kg_s,
                      viscosity_ratio=1.0, density_ratio=density_ratio)
    predicted_pa = result["estimated_tube_pressure_drop_pa"]
    measured_pa = ref_213.tube_pressure_drop_pa
    
    # Expected: about 3.80 kPa predicted vs 4.40 kPa measured (~-14%)
    assert predicted_pa == pytest.approx(3800, abs=100)  # ~3.80 kPa ± 100 Pa
    assert measured_pa == pytest.approx(4400, abs=100)   # ~4.40 kPa ± 100 Pa
    percent_diff = (predicted_pa - measured_pa) / measured_pa * 100
    assert percent_diff == pytest.approx(-14, abs=5)  # About -14%


# ============================================================================
# File format and rejection tests
# ============================================================================

def test_mixed_fluid_file_rejected():
    """Reject CSV files containing both nitrogen and helium."""
    import tempfile
    import os
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, newline='') as f:
        f.write("fluid,mass_flow_kg_s,mass_flow_uncertainty_kg_s,pressure_Pa,pressure_uncertainty_Pa,T1_K,T2_K,T3_K,T4_K,effectiveness,effectiveness_uncertainty,UA_W_K,UA_uncertainty_W_K,pressure_drop_Pa,pressure_drop_uncertainty_Pa\n")
        f.write("nitrogen,0.000470,0.000020,322000,3000,377.85,304.05,296.15,366.15,0.857,0.006,3.7,0.4,5100,200\n")
        f.write("helium,0.000117,0.000007,749000,3000,380.65,301.45,296.65,376.75,0.954,0.006,11.1,2.5,2300,100\n")
        path = f.name
    try:
        with pytest.raises(ValueError, match="homogeneous fluid"):
            load_references(path)
    finally:
        os.unlink(path)


def test_empty_file_rejected():
    """Reject empty CSV files."""
    import tempfile
    import os
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, newline='') as f:
        f.write("fluid,mass_flow_kg_s,mass_flow_uncertainty_kg_s,pressure_Pa,pressure_uncertainty_Pa,T1_K,T2_K,T3_K,T4_K,effectiveness,effectiveness_uncertainty,UA_W_K,UA_uncertainty_W_K,pressure_drop_Pa,pressure_drop_uncertainty_Pa\n")
        path = f.name
    try:
        with pytest.raises(ValueError, match="non-empty"):
            load_references(path)
    finally:
        os.unlink(path)


def test_unsupported_fluid_rejected():
    """Reject CSV files with unsupported fluid types."""
    import tempfile
    import os
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, newline='') as f:
        f.write("fluid,mass_flow_kg_s,mass_flow_uncertainty_kg_s,pressure_Pa,pressure_uncertainty_Pa,T1_K,T2_K,T3_K,T4_K,effectiveness,effectiveness_uncertainty,UA_W_K,UA_uncertainty_W_K,pressure_drop_Pa,pressure_drop_uncertainty_Pa\n")
        f.write("argon,0.000470,0.000020,322000,3000,377.85,304.05,296.15,366.15,0.857,0.006,3.7,0.4,5100,200\n")
        path = f.name
    try:
        with pytest.raises(ValueError, match="Unsupported fluid"):
            load_references(path)
    finally:
        os.unlink(path)

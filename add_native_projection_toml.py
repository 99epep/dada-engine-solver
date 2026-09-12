#!/usr/bin/env python3
from pathlib import Path

ROOT = Path.cwd()
if ROOT.name != "dada-engine-solver":
    raise SystemExit("Run this script from the dada-engine-solver repository root.")

MODULE = '"""Fast projection-only coupler-point four-bar kinematics."""\n\nfrom __future__ import annotations\n\nfrom dataclasses import dataclass, field\nimport math\n\nimport numpy as np\nfrom scipy.optimize import brentq\n\nfrom dada_solver.geometry import CylinderVolumeLimits\n\n\ndef _positive(name: str, value: float) -> None:\n    if not math.isfinite(value) or value <= 0.0:\n        raise ValueError(f"{name} must be finite and positive.")\n\n\n@dataclass(frozen=True, slots=True)\nclass CouplerProjectionSideDesign:\n    coupler_ratio: float\n    rocker_ratio: float\n    pivot_x_ratio: float\n    pivot_y_ratio: float\n    output_along_ratio: float\n    output_normal_ratio: float\n    axis_angle_degrees: float\n    assembly_branch: int\n\n    def __post_init__(self) -> None:\n        _positive("coupler ratio", self.coupler_ratio)\n        _positive("rocker ratio", self.rocker_ratio)\n        for name, value in (\n            ("pivot x ratio", self.pivot_x_ratio),\n            ("pivot y ratio", self.pivot_y_ratio),\n            ("output along ratio", self.output_along_ratio),\n            ("output normal ratio", self.output_normal_ratio),\n            ("axis angle", self.axis_angle_degrees),\n        ):\n            if not math.isfinite(value):\n                raise ValueError(f"{name} must be finite.")\n        if self.assembly_branch not in (-1, 1):\n            raise ValueError("assembly_branch must be -1 or 1.")\n\n\n@dataclass(frozen=True, slots=True)\nclass SharedCrankCouplerProjectionDesign:\n    small: CouplerProjectionSideDesign\n    large: CouplerProjectionSideDesign\n\n\n@dataclass(frozen=True, slots=True)\nclass _ProjectionSide:\n    design: CouplerProjectionSideDesign\n    scan_count: int = 721\n    _minimum: float = field(init=False, repr=False)\n    _maximum: float = field(init=False, repr=False)\n\n    def __post_init__(self) -> None:\n        angles = np.linspace(0.0, 2.0 * math.pi, self.scan_count, endpoint=False)\n        values = np.empty(self.scan_count)\n        derivatives = np.empty(self.scan_count)\n        for index, angle in enumerate(angles):\n            values[index], derivatives[index] = self.evaluate(float(angle))\n\n        roots: list[float] = []\n        for index, angle in enumerate(angles):\n            next_index = (index + 1) % self.scan_count\n            end = float(angles[next_index])\n            if next_index == 0:\n                end = 2.0 * math.pi\n            left = derivatives[index]\n            right = derivatives[next_index]\n            if left == 0.0:\n                roots.append(float(angle))\n            elif left * right < 0.0:\n                roots.append(\n                    brentq(\n                        lambda candidate: self.evaluate(candidate)[1],\n                        float(angle),\n                        end,\n                    )\n                    % (2.0 * math.pi)\n                )\n\n        extrema = [self.evaluate(angle)[0] for angle in roots]\n        extrema.extend((float(np.min(values)), float(np.max(values))))\n        minimum = min(extrema)\n        maximum = max(extrema)\n        if maximum - minimum <= 1.0e-12:\n            raise ValueError("Projected four-bar has zero usable stroke.")\n        object.__setattr__(self, "_minimum", minimum)\n        object.__setattr__(self, "_maximum", maximum)\n\n    @property\n    def stroke(self) -> float:\n        return self._maximum - self._minimum\n\n    def evaluate(self, theta: float) -> tuple[float, float]:\n        d = self.design\n\n        pin_x = math.cos(theta)\n        pin_y = math.sin(theta)\n        pin_dx = -math.sin(theta)\n        pin_dy = math.cos(theta)\n\n        delta_x = d.pivot_x_ratio - pin_x\n        delta_y = d.pivot_y_ratio - pin_y\n        distance = math.hypot(delta_x, delta_y)\n        if distance <= 1.0e-12:\n            raise ValueError("Four-bar circle centers coincide.")\n\n        direction_x = delta_x / distance\n        direction_y = delta_y / distance\n        along = (\n            d.coupler_ratio**2\n            - d.rocker_ratio**2\n            + distance**2\n        ) / (2.0 * distance)\n        height_squared = d.coupler_ratio**2 - along**2\n        if height_squared <= 1.0e-20:\n            raise ValueError("Four-bar cannot close or is at toggle.")\n        height = math.sqrt(height_squared)\n\n        branch = d.assembly_branch\n        joint_x = pin_x + along * direction_x - branch * height * direction_y\n        joint_y = pin_y + along * direction_y + branch * height * direction_x\n\n        coupler_x = joint_x - pin_x\n        coupler_y = joint_y - pin_y\n        rocker_x = joint_x - d.pivot_x_ratio\n        rocker_y = joint_y - d.pivot_y_ratio\n\n        determinant = coupler_x * rocker_y - coupler_y * rocker_x\n        if abs(determinant) <= 1.0e-12:\n            raise ValueError("Four-bar velocity closure is singular.")\n\n        velocity_rhs = coupler_x * pin_dx + coupler_y * pin_dy\n        joint_dx = velocity_rhs * rocker_y / determinant\n        joint_dy = -velocity_rhs * rocker_x / determinant\n\n        unit_x = coupler_x / d.coupler_ratio\n        unit_y = coupler_y / d.coupler_ratio\n        unit_dx = (joint_dx - pin_dx) / d.coupler_ratio\n        unit_dy = (joint_dy - pin_dy) / d.coupler_ratio\n\n        output_x = (\n            pin_x\n            + d.output_along_ratio * unit_x\n            - d.output_normal_ratio * unit_y\n        )\n        output_y = (\n            pin_y\n            + d.output_along_ratio * unit_y\n            + d.output_normal_ratio * unit_x\n        )\n        output_dx = (\n            pin_dx\n            + d.output_along_ratio * unit_dx\n            - d.output_normal_ratio * unit_dy\n        )\n        output_dy = (\n            pin_dy\n            + d.output_along_ratio * unit_dy\n            + d.output_normal_ratio * unit_dx\n        )\n\n        axis = math.radians(d.axis_angle_degrees)\n        axis_x = math.cos(axis)\n        axis_y = math.sin(axis)\n        return (\n            output_x * axis_x + output_y * axis_y,\n            output_dx * axis_x + output_dy * axis_y,\n        )\n\n    def normalized(self, theta: float) -> tuple[float, float]:\n        coordinate, derivative = self.evaluate(theta)\n        return (\n            (coordinate - self._minimum) / self.stroke,\n            derivative / self.stroke,\n        )\n\n\n@dataclass(frozen=True, slots=True)\nclass SharedCrankCouplerProjectionKinematics:\n    design: SharedCrankCouplerProjectionDesign\n    small_volume_limits: CylinderVolumeLimits\n    large_volume_limits: CylinderVolumeLimits\n    _small: _ProjectionSide = field(init=False, repr=False)\n    _large: _ProjectionSide = field(init=False, repr=False)\n\n    def __post_init__(self) -> None:\n        object.__setattr__(self, "_small", _ProjectionSide(self.design.small))\n        object.__setattr__(self, "_large", _ProjectionSide(self.design.large))\n\n    @property\n    def small_physical_stroke(self):\n        return None\n\n    @property\n    def large_physical_stroke(self):\n        return None\n\n    def breakpoint_angles(self) -> tuple[float, ...]:\n        return ()\n\n    @staticmethod\n    def _volume(side, limits, theta):\n        fraction, derivative = side.normalized(theta)\n        return (\n            limits.minimum + limits.swept * fraction,\n            limits.swept * derivative,\n        )\n\n    def cylinder_volumes_and_derivatives(self, theta: float):\n        sv, sd = self._volume(self._small, self.small_volume_limits, theta)\n        lv, ld = self._volume(self._large, self.large_volume_limits, theta)\n        return sv, lv, sd, ld\n\n    def small_cylinder_volume(self, theta: float) -> float:\n        return self._volume(self._small, self.small_volume_limits, theta)[0]\n\n    def large_cylinder_volume(self, theta: float) -> float:\n        return self._volume(self._large, self.large_volume_limits, theta)[0]\n\n    def small_cylinder_volume_derivative(self, theta: float) -> float:\n        return self._volume(self._small, self.small_volume_limits, theta)[1]\n\n    def large_cylinder_volume_derivative(self, theta: float) -> float:\n        return self._volume(self._large, self.large_volume_limits, theta)[1]\n\n\ndef shared_crank_coupler_projection_kinematics(\n    design: SharedCrankCouplerProjectionDesign,\n    small_volume_limits: CylinderVolumeLimits,\n    large_volume_limits: CylinderVolumeLimits,\n) -> SharedCrankCouplerProjectionKinematics:\n    return SharedCrankCouplerProjectionKinematics(\n        design=design,\n        small_volume_limits=small_volume_limits,\n        large_volume_limits=large_volume_limits,\n    )\n'
TOML = '# Native thermodynamic test of the compact projection-only four-bar solution.\n\n[metadata]\nexample_data = true\ndescription = "Motor target approximated by compact coupler-point four-bars"\n\n[working_gas]\ngas_constant = 287.0\nheat_capacity_cp = 1004.5\nheat_capacity_cv = 717.5\n\n[reservoirs]\ncold_temperature = 298.15\nhot_temperature = 598.15\n\n[operation]\nangular_speed = -12.566370614359172\n\n[charge]\npressure = 100000.0\ntemperature = 298.15\n\n[geometry]\nsmall_cylinder_minimum_volume = 8.316831683168318e-06\nsmall_cylinder_maximum_volume = 0.00084\nlarge_cylinder_minimum_volume = 9.900990099009901e-06\nlarge_cylinder_maximum_volume = 0.001\ncold_heat_exchanger_volume = 1.9801980198019803e-05\nhot_heat_exchanger_volume = 1.9801980198019803e-05\n\n[heat_transfer]\ncold_ua = 49.504275427542765\nhot_ua = 49.504275427542765\n\n[hydraulics]\norifice_pressure_regularization = 0.03\nlarge_to_hot_cda = 0.0001663343654365437\nsmall_to_cold_cda = 0.0002970256525652566\nhot_to_small_valve_cda = 0.0001663343654365437\ncold_to_large_valve_cda = 0.0005940513051305131\n\n[valves]\nmodel = "continuous_ideal_diode"\n\n[valves.hot_to_small]\nopening_pressure_difference = 0.0\nclosing_pressure_difference = 0.0\n\n[valves.cold_to_large]\nopening_pressure_difference = 0.0\nclosing_pressure_difference = 0.0\n\n[kinematics]\ntype = "shared_crank_coupler_projection"\n\n[kinematics.small_four_bar]\ncoupler_ratio = 1.6563469\nrocker_ratio = 5.1468443\npivot_x_ratio = -4.5545429\npivot_y_ratio = -2.1174820\noutput_along_ratio = 1.5159801\noutput_normal_ratio = -0.41242206\naxis_angle_degrees = -132.0556\nassembly_branch = -1\n\n[kinematics.large_four_bar]\ncoupler_ratio = 1.8978982\nrocker_ratio = 8.5625583\npivot_x_ratio = 6.5305494107\npivot_y_ratio = 5.2178614049\noutput_along_ratio = 1.8306391\noutput_normal_ratio = 0.33169925\naxis_angle_degrees = 23.3905\nassembly_branch = 1\n\n[validity]\nmaximum_pressure_equalization_error = 0.05\nmaximum_mach_number = 0.2\nmaximum_isothermality_error = 0.05\nmaximum_compressibility_deviation = 0.01\nmaximum_cp_variation = 0.01\n\n[numerical]\nintegration_method = "LSODA"\nintegration_relative_tolerance = 1.0e-7\nintegration_absolute_tolerance = 1.0e-9\nmaximum_step_angle_degrees = 0.5\nmaximum_valve_events_per_cycle = 100\nmaximum_cycles = 100\nperiodic_relative_tolerance = 1.0e-6\nperiodic_mass_absolute_tolerance = 1.0e-12\nperiodic_energy_absolute_tolerance = 1.0e-6\n'

module_path = ROOT / "src/dada_solver/coupler_projection.py"
module_path.write_text(MODULE, encoding="utf-8")
print("wrote", module_path)

config_path = ROOT / "src/dada_solver/configuration.py"
text = config_path.read_text(encoding="utf-8")

needle = "from dada_solver.four_bar import SharedCrankRockerDesign\n"
replacement = (
    needle
    + "from dada_solver.coupler_projection import (\n"
    + "    CouplerProjectionSideDesign, SharedCrankCouplerProjectionDesign,\n"
    + ")\n"
)
if needle not in text:
    raise SystemExit("configuration import marker not found")
text = text.replace(needle, replacement, 1)

needle = "    free_kinematics: FreeKinematicsConfiguration | None = None\n"
replacement = (
    needle
    + "    shared_coupler_projection_design: SharedCrankCouplerProjectionDesign | None = None\n"
)
if needle not in text:
    raise SystemExit("configuration field marker not found")
text = text.replace(needle, replacement, 1)

needle = '''        elif self.kinematics_type == "shared_crank_rocker":
            if self.shared_four_bar_design is None:
                raise ValueError("Shared-crank rocker design is required.")
'''
replacement = '''        elif self.kinematics_type == "shared_crank_coupler_projection":
            if self.shared_coupler_projection_design is None:
                raise ValueError("Shared-crank coupler-projection design is required.")
        elif self.kinematics_type == "shared_crank_rocker":
            if self.shared_four_bar_design is None:
                raise ValueError("Shared-crank rocker design is required.")
'''
if needle not in text:
    raise SystemExit("configuration validation marker not found")
text = text.replace(needle, replacement, 1)

needle = "            free_kinematics=_load_free_kinematics(kinematics_data, machine_volumes),\n"
replacement = (
    needle
    + "            shared_coupler_projection_design="
      "_load_shared_coupler_projection_design(kinematics_data),\n"
)
if needle not in text:
    raise SystemExit("configuration load marker not found")
text = text.replace(needle, replacement, 1)

helper = r'''

def _load_shared_coupler_projection_design(
    data: dict[str, object],
) -> SharedCrankCouplerProjectionDesign | None:
    if data.get("type") != "shared_crank_coupler_projection":
        return None

    def side(name: str) -> CouplerProjectionSideDesign:
        raw = data[name]
        if not isinstance(raw, dict):
            raise ValueError(f"{name} must be a table.")
        return CouplerProjectionSideDesign(
            coupler_ratio=float(raw["coupler_ratio"]),
            rocker_ratio=float(raw["rocker_ratio"]),
            pivot_x_ratio=float(raw["pivot_x_ratio"]),
            pivot_y_ratio=float(raw["pivot_y_ratio"]),
            output_along_ratio=float(raw["output_along_ratio"]),
            output_normal_ratio=float(raw["output_normal_ratio"]),
            axis_angle_degrees=float(raw["axis_angle_degrees"]),
            assembly_branch=int(raw["assembly_branch"]),
        )

    return SharedCrankCouplerProjectionDesign(
        small=side("small_four_bar"),
        large=side("large_four_bar"),
    )
'''
if "def _load_shared_coupler_projection_design(" not in text:
    text += helper

config_path.write_text(text, encoding="utf-8")
print("patched", config_path)

factory_path = ROOT / "src/dada_solver/factory.py"
text = factory_path.read_text(encoding="utf-8")
needle = "from dada_solver.periodic import SuccessiveCycleSolver\n"
replacement = (
    needle
    + "from dada_solver.coupler_projection import "
      "shared_crank_coupler_projection_kinematics\n"
)
if needle not in text:
    raise SystemExit("factory import marker not found")
text = text.replace(needle, replacement, 1)

needle = '    elif configuration.kinematics_type == "published_f65_opposed":\n'
replacement = '''    elif configuration.kinematics_type == "shared_crank_coupler_projection":
        assert configuration.shared_coupler_projection_design is not None
        kinematics = shared_crank_coupler_projection_kinematics(
            configuration.shared_coupler_projection_design,
            configuration.machine_volumes.small_cylinder,
            configuration.machine_volumes.large_cylinder,
        )
    elif configuration.kinematics_type == "published_f65_opposed":
'''
if needle not in text:
    raise SystemExit("factory branch marker not found")
text = text.replace(needle, replacement, 1)
factory_path.write_text(text, encoding="utf-8")
print("patched", factory_path)

toml_path = ROOT / "outputs/motor_champion_four_bar_projection.toml"
toml_path.write_text(TOML, encoding="utf-8")
print("wrote", toml_path)

print()
print("Then run:")
print("  python -m pytest")
print("  python -m dada_solver.cli outputs/motor_champion_four_bar_projection.toml")

"""Family-specific coordinate ownership and cheap construction checks."""
from dataclasses import replace, fields
import math
from dada_solver.free_kinematics import FreeMotionDefinition, FreeKinematicsConfiguration, FreeKinematics
from dada_solver.four_bar import shared_crank_rocker_kinematics, SharedCrankRockerDesign
from dada_solver.kinematics import KinematicConstraintViolation
from dada_solver.machine import MachineDesign
from dada_solver.sizing.design import DesignParameter, DesignPoint, apply_design_point


class PreflightRejection(ValueError):
    def __init__(self, status, reason, diagnostics=None):
        super().__init__(reason)
        self.status, self.diagnostics = status, diagnostics


COMMON = {name: item for item in DesignParameter for name in ['common.'+item.value]
          if item not in {DesignParameter.SMALL_PHASE_OFFSET, DesignParameter.COMMON_LAMBDA_TARGET,
                          DesignParameter.SMALL_LAMBDA_TARGET, DesignParameter.LARGE_LAMBDA_TARGET,
                          DesignParameter.ADIABATIC_SECTOR_FRACTION}}
DERIVED_EXCHANGER = {'common.cold_ua', 'common.hot_ua', 'common.cold_heat_exchanger_volume',
    'common.hot_heat_exchanger_volume', 'common.large_to_hot_cda', 'common.small_to_cold_cda'}
FOUR_BAR_FIELDS = {f.name for f in fields(SharedCrankRockerDesign) if 'branch' not in f.name}
MICROTUBE_FIELDS = {'tube_length_m', 'inner_diameter_m', 'wall_thickness_m', 'pitch_m',
                   'header_depth_m', 'additional_internal_volume_m3'}


def validate_ownership(names, families, free_settings):
    if families['kinematics'] not in ('free', 'four_bar'):
        raise ValueError('The initial campaign adapters support free or four_bar kinematics.')
    if families['exchanger'] not in ('reservoir', 'microtube'):
        raise ValueError('Unknown exchanger family.')
    if families['exchanger'] == 'microtube' and DERIVED_EXCHANGER.intersection(names):
        raise ValueError('Microtube geometry owns derived UA, hold-up and hydraulic resistance.')
    if 'operation.frequency_hz' in names and 'common.angular_speed' in names:
        raise ValueError('Frequency and angular speed cannot both be independent.')
    if {'common.charge_pressure', 'common.total_gas_mass'}.issubset(names):
        raise ValueError('Choose pressure or gas inventory, not both.')
    for side in ('small', 'large'):
        if {f'common.{side}_clearance_volume', f'common.{side}_clearance_ratio'}.issubset(names):
            raise ValueError('Clearance volume and ratio cannot both be independent.')
    for name in names:
        if name in COMMON or name == 'operation.frequency_hz':
            continue
        parts = name.split('.')
        if len(parts) == 3 and parts[0] == 'free' and parts[1] in ('small','large') and families['kinematics'] == 'free':
            coordinates = free_settings[parts[1]+'_coordinates']
            if parts[2].isdigit() and str(int(parts[2])) == parts[2] and 0 <= int(parts[2]) < len(coordinates):
                continue
        if len(parts) == 2 and parts[0] == 'four_bar' and families['kinematics'] == 'four_bar' and parts[1] in FOUR_BAR_FIELDS:
            continue
        if len(parts) == 3 and parts[0] == 'microtube' and families['exchanger'] == 'microtube' and parts[1] in ('heat_in','heat_out') and parts[2] in MICROTUBE_FIELDS:
            continue
        raise ValueError(f'Unsupported or wrong-family parameter: {name}')


def microtube_design(bank, inputs, valve_cda, physical, side):
    """Geometry adapter; integer tube count is fixed in a continuous campaign."""
    from dada_solver.exchangers.microtube import MicrotubeExchanger
    updates = {name.split('.')[2]: value for name, value in physical.items()
               if name.startswith('microtube.'+side+'.')}
    if not set(updates) <= MICROTUBE_FIELDS:
        raise PreflightRejection('invalid_exchanger', 'Unknown or derived microtube coordinate.')
    try:
        design = MicrotubeExchanger(replace(bank, **updates), inputs, valve_cda)
        design.build()
        return design
    except ValueError as error:
        raise PreflightRejection('invalid_exchanger', str(error)) from error


class FamilyDesignAdapter:
    def __init__(self, configuration, families, free_settings):
        self.configuration, self.families, self.free_settings = configuration, families, free_settings

    def build(self, physical):
        try:
            validate_ownership(physical, self.families, self.free_settings)
            updates = {COMMON[name]: value for name, value in physical.items() if name in COMMON}
            if 'operation.frequency_hz' in physical:
                frequency = physical['operation.frequency_hz']
                if frequency <= 0: raise ValueError('Frequency must be positive.')
                updates[DesignParameter.ANGULAR_SPEED] = math.copysign(2*math.pi*frequency, self.configuration.angular_speed)
            config = apply_design_point(self.configuration, DesignPoint(updates))
            if config.motor_operation != self.configuration.motor_operation:
                raise ValueError('A campaign cannot cross operating direction.')
        except (ValueError, KeyError, TypeError) as error:
            raise PreflightRejection('invalid_parameterization', str(error)) from error
        try:
            if self.families['kinematics'] == 'free':
                definitions = []
                for side in ('small', 'large'):
                    coordinates = list(self.free_settings[side+'_coordinates'])
                    for i in range(len(coordinates)):
                        coordinates[i] = physical.get(f'free.{side}.{i}', coordinates[i])
                    limits = getattr(config.machine_volumes, side+'_cylinder')
                    bounds = self.free_settings.get(side+'_limits', {})
                    definitions.append(FreeMotionDefinition.from_shape_coordinates(coordinates,
                        limits.minimum, limits.maximum, **bounds))
                kinematics = FreeKinematics(FreeKinematicsConfiguration(*definitions))
                kinematics.require_feasible()
            else:
                if config.shared_four_bar_design is None:
                    raise ValueError('Four-bar campaign requires a shared-crank base configuration.')
                design = replace(config.shared_four_bar_design, **{
                    name.split('.')[1]: value for name, value in physical.items() if name.startswith('four_bar.')})
                kinematics = shared_crank_rocker_kinematics(design,
                    config.machine_volumes.small_cylinder, config.machine_volumes.large_cylinder,
                    ground_distance=config.four_bar_ground_distance,
                    crank_angle_offset=math.radians(config.four_bar_crank_angle_offset_degrees),
                    crank_direction=config.four_bar_crank_direction)
        except (ValueError, ArithmeticError) as error:
            from dataclasses import asdict
            diagnostics = ([asdict(x) for x in error.diagnostics]
                           if isinstance(error, KinematicConstraintViolation) else None)
            raise PreflightRejection('invalid_kinematics', str(error), diagnostics) from error
        return MachineDesign(config, kinematics=kinematics)

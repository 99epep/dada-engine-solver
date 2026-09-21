"""Single authoritative final-trajectory reconstruction for wall diagnostics.

Numerical acceptance facts are retained separately from reporting aggregates.
No sample, acceptance criterion or detailed scientific report is omitted.
"""
from dataclasses import dataclass
import numpy as np
from dada_solver.dynamics import InstantaneousPoint
from dada_solver.state import ThermodynamicState
from dada_solver.exchangers.air_wall import AirWallExchanger, WallThermalPoint


@dataclass(frozen=True, slots=True)
class DiagnosticSample:
    point: InstantaneousPoint
    walls: tuple[WallThermalPoint, WallThermalPoint]
    hydraulic_diagnostics: tuple


@dataclass(frozen=True, slots=True)
class WallTrajectoryReplay:
    """Evaluation-local numerical facts, bound to one model and trajectory."""
    wrapper: object
    cycle: object
    trajectory: np.ndarray
    samples: tuple[DiagnosticSample, ...]
    pressures: np.ndarray
    temperatures: np.ndarray

    def require(self, *, model=None, cycle=None, wrapper=None, angles=None, trajectory=None):
        if ((model is not None and model is not self.wrapper.model)
                or (cycle is not None and cycle is not self.cycle)
                or (wrapper is not None and wrapper is not self.wrapper)
                or (angles is not None and angles is not self.cycle.angles)
                or (trajectory is not None and trajectory is not self.trajectory)):
            raise ValueError('Diagnostic replay belongs to a different model or trajectory.')


def replay_wall_trajectory(wrapper, cycle, trajectory):
    """One hydraulic point and one film evaluation per side and stored sample.

    Hydraulic upstream diagnostics use their original half-length/upstream
    temperature; they must not be confused with the full-length gas-side film.
    """
    if not cycle.completed:
        raise ValueError('Diagnostic replay requires a completed cycle.')
    if trajectory.ndim != 2 or trajectory.shape[0] != 15 or trajectory.shape[1] != len(cycle.angles):
        raise ValueError('Expected a complete fifteen-row wall trajectory.')
    if not np.array_equal(cycle.states, trajectory[:8]):
        raise ValueError('Cycle states do not match the supplied trajectory.')
    rows=[]
    for i,(angle,topology) in enumerate(zip(cycle.angles,cycle.topologies,strict=True)):
        values=trajectory[:,i]
        state=ThermodynamicState.from_array(values[:8])
        point=wrapper.model.instantaneous_point(float(angle),state,topology)
        contexts=wrapper._flow_contexts_from_point(point)
        thermal=[];hydraulic=[]
        for index,wall,context,ports in zip((2,3),(wrapper.heat_in,wrapper.heat_out),
                contexts,(((0,2),(2,1)),((1,3),(3,0)))):
            temperature=point.temperatures[index]
            energy=values[index+6]
            if type(wall) is AirWallExchanger:
                evaluated=wall.thermal_point(temperature,energy,context=context)
            else:
                # Preserve custom wall protocols and overridden public rates.
                contextual=getattr(wall,'requires_flow_context',False)
                rates=wall.rates(temperature,energy,context=context) if contextual else wall.rates(temperature,energy)
                diagnostics=wall.gas_film.evaluate(temperature,energy/wall.wall_capacity_j_k,context)[1] if contextual else ()
                evaluated=WallThermalPoint(rates['gas_heat_w'],rates['air_heat_w'],
                    rates['wall_energy_rate_w'],rates.get('air_outlet_temperature_k'),
                    energy/wall.wall_capacity_j_k,diagnostics)
            thermal.append(evaluated)
            if not getattr(wall,'requires_flow_context',False):
                hydraulic.extend((None,None));continue
            film=wall.gas_film
            for passage,pair in zip(context['passages'],ports):
                flow,p1,p2=passage
                upstream=pair[0] if flow>=0 else pair[1]
                hydraulic.append(film.model.diagnose(film.bank,flow,p1,p2,
                    point.temperatures[upstream],frequency=context['frequency_hz'],
                    length=film.bank.tube_length_m/2))
        rows.append(DiagnosticSample(point,tuple(thermal),tuple(hydraulic)))
    pressures=np.array([r.point.pressures for r in rows]).T
    temperatures=np.array([r.point.temperatures for r in rows]).T
    pressures.flags.writeable=False;temperatures.flags.writeable=False
    return WallTrajectoryReplay(wrapper,cycle,trajectory,tuple(rows),pressures,temperatures)

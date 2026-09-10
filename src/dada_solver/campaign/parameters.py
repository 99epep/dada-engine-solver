"""Bounded, named campaign coordinates independent of the legacy sizing enum."""
from dataclasses import dataclass
import math


@dataclass(frozen=True, slots=True)
class ContinuousParameter:
    name: str
    lower: float
    upper: float
    initial: float
    transform: str = 'linear'

    def __post_init__(self):
        if not self.name or not all(math.isfinite(x) for x in (self.lower, self.upper, self.initial)):
            raise ValueError('Parameters require a stable name and finite values.')
        if not self.lower < self.upper or not self.lower <= self.initial <= self.upper or not math.isfinite(self.upper-self.lower):
            raise ValueError('Parameter bounds must be ordered and contain the initial value.')
        if self.transform not in ('linear', 'log') or (self.transform == 'log' and self.lower <= 0):
            raise ValueError('Transform must be linear or log with positive log bounds.')

    def decode(self, normalized):
        u = float(normalized)
        if not math.isfinite(u) or not 0 <= u <= 1:
            raise ValueError('Normalized coordinates must lie in [0, 1].')
        if u == 0: return self.lower
        if u == 1: return self.upper
        if self.transform == 'log':
            return math.exp((1-u)*math.log(self.lower)+u*math.log(self.upper))
        return (1-u)*self.lower+u*self.upper

    def encode(self, physical):
        x = float(physical)
        if not math.isfinite(x) or not self.lower <= x <= self.upper:
            raise ValueError('Physical value is outside its parameter bounds.')
        if x == self.lower: return 0.0
        if x == self.upper: return 1.0
        if self.transform == 'log':
            return (math.log(x)-math.log(self.lower))/(math.log(self.upper)-math.log(self.lower))
        return (x-self.lower)/(self.upper-self.lower)


@dataclass(frozen=True, slots=True)
class ParameterSpace:
    parameters: tuple[ContinuousParameter, ...]

    def __post_init__(self):
        object.__setattr__(self, 'parameters', tuple(self.parameters))
        names = [p.name for p in self.parameters]
        if not names or len(names) != len(set(names)):
            raise ValueError('A nonempty space requires unique parameter names.')

    @property
    def initial_coordinates(self):
        return tuple(p.encode(p.initial) for p in self.parameters)

    def decode(self, coordinates):
        if len(coordinates) != len(self.parameters):
            raise ValueError('Coordinate dimension does not match parameter space.')
        return {p.name: p.decode(x) for p, x in zip(self.parameters, coordinates)}

    def encode(self, physical):
        if set(physical) != {p.name for p in self.parameters}:
            raise ValueError('Physical parameter names do not match the space.')
        return tuple(p.encode(physical[p.name]) for p in self.parameters)

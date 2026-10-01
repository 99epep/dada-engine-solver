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
class IntegerParameter:
    """Inclusive integer interval with explicitly recorded nearest-even decoding."""
    name: str
    lower: int
    upper: int
    initial: int
    encoding: str = 'nearest_even_v1'

    def __post_init__(self):
        if not self.name or any(type(x) is not int for x in (self.lower, self.upper, self.initial)):
            raise ValueError('Integer bounds and initial value must be integers, not bool or float.')
        if not 0 < self.lower < self.upper or not self.lower <= self.initial <= self.upper:
            raise ValueError('Integer counts require positive ordered bounds containing the initial value.')
        if self.upper > 2**53 or self.encoding != 'nearest_even_v1':
            raise ValueError('Unsupported integer range or encoding.')

    def decode(self, normalized):
        u = float(normalized)
        if not math.isfinite(u) or not 0 <= u <= 1:
            raise ValueError('Normalized coordinates must lie in [0, 1].')
        return int(round((1-u)*self.lower + u*self.upper))

    def encode(self, physical):
        if type(physical) is not int or not self.lower <= physical <= self.upper:
            raise ValueError('Physical count must be an integer within its bounds.')
        return (physical-self.lower)/(self.upper-self.lower)


@dataclass(frozen=True, slots=True)
class ChoiceParameter:
    """Explicit categorical bins; their order is an encoding, not a distance."""
    name: str
    choices: tuple
    initial: str | int | float | bool

    def __post_init__(self):
        if not self.name or not isinstance(self.choices, (list, tuple)) or not self.choices:
            raise ValueError('Choices require a name and a nonempty list or tuple.')
        choices = tuple(self.choices)
        if any(type(v) not in (str, int, float, bool) or
               (isinstance(v, float) and not math.isfinite(v)) for v in choices):
            raise ValueError('Choices must be finite scalar values.')
        if len(set(choices)) != len(choices) or self.initial not in choices:
            raise ValueError('Choices must be distinct and contain the initial value.')
        object.__setattr__(self, 'choices', choices)

    def decode(self, normalized):
        u = float(normalized)
        if not math.isfinite(u) or not 0 <= u <= 1:
            raise ValueError('Normalized coordinates must lie in [0, 1].')
        return self.choices[min(int(u * len(self.choices)), len(self.choices)-1)]

    def encode(self, physical):
        if physical not in self.choices:
            raise ValueError('Physical value is not one of the declared choices.')
        return (self.choices.index(physical) + .5) / len(self.choices)


def parameter_from_mapping(data):
    settings = dict(data)
    kind = settings.pop('kind', 'continuous')
    if kind == 'continuous': return ContinuousParameter(**settings)
    if kind == 'integer': return IntegerParameter(**settings)
    if kind == 'choice': return ChoiceParameter(**settings)
    raise ValueError(f'Unknown parameter kind: {kind}')


@dataclass(frozen=True, slots=True)
class ParameterSpace:
    parameters: tuple[ContinuousParameter | IntegerParameter | ChoiceParameter, ...]
    allow_empty: bool = False

    def __post_init__(self):
        object.__setattr__(self, 'parameters', tuple(self.parameters))
        names = [p.name for p in self.parameters]
        if (not names and not self.allow_empty) or len(names) != len(set(names)):
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

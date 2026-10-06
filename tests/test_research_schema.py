"""Validation, ownership and exact integer identities for the first study protocol."""
import hashlib
import json
from pathlib import Path

import pytest

from dada_solver.campaign.parameters import IntegerParameter
from dada_solver.research.schema import load_study, compile_study, candidate_for_values




def test_integer_decoding_ties_and_endpoints():
    p = IntegerParameter('count', 2, 6, 4)
    assert [p.decode(u) for u in (0,.125,.375,.625,.875,1)] == [2,2,4,4,6,6]
    for k in range(2,7): assert p.decode(p.encode(k)) == k
    for value in (True,3.0,0,7):
        with pytest.raises(ValueError): p.encode(value)
    for value in (-.01,1.01,float('nan')):
        with pytest.raises(ValueError): p.decode(value)


@pytest.mark.parametrize('values', [(True,6,4),(2.,6,4),(0,6,4),(2,2,2),(2,6,7)])
def test_invalid_integer_parameter(values):
    with pytest.raises(ValueError): IntegerParameter('count',*values)

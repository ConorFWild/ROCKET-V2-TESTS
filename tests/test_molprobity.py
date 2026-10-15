"""
Unit test that getting molprobity distributions is working correctly
"""

import pytest

import numpy as np

from rocket_testing.types import TStructure, TAtomID, TGrid, TReflections, TMolprobity
from rocket_testing.molprobity import molprobity 

# @pytest.mark.skip(reason="Not Implemented")
def test_molprobity_known_values(structure: TStructure, reflections: TReflections, known_molprobity: TMolprobity):

    # Test against known values
    _molprobity: TMolprobity = molprobity(structure, reflections)

    for k in _molprobity.keys:
        
        assert np.allclose(_molprobity.__getattribute__(k), known_molprobity.__getattribute__(k), rtol=0.15)
    ...
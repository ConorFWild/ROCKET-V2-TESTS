"""
Unit test that rfree is being calculated correctly
"""
import pytest 

import numpy as np

from rocket_testing.rfree import rfree
from rocket_testing.types import TStructure, TAtomID, TGrid, TReflections

@pytest.mark.skip(reason="Not Implemented")
def test_rfree_known_values(structure: TStructure, known_rfree: float):

    # Test against known rfree
    _rfree: float = rfree(structure)

    assert np.allclose(_rfree, known_rfree)
    ...
"""
Unit test that rsccs are working correctly
"""
import pytest

import numpy as np

from rocket_testing.types import TStructure, TAtomID, TGrid, TReflections
from rocket_testing.rscc import get_rscc

def test_rscc_model_mtz(structure: TStructure, reflections: TReflections):
    get_rscc(structure, mtz=reflections)

def test_rscc_model_xmap(structure: TStructure, grid: TGrid):
    get_rscc(structure, xmap=grid)

@pytest.mark.skip(reason="Not Implemented")
def test_rscc_self_all_atom(structure: TStructure, all_atom_selection: list[TAtomID]):

    # Test rscc to self is 1

    structure_copy: TStructure = TStructure.clone()
    rscc: float = rscc_st_to_st(TStructure, structure_copy, all_atom_selection)

    assert np.allclose(rscc, 1.0)

@pytest.mark.skip(reason="Not Implemented")
def test_rscc_self_residues(structure: TStructure, residue_selection: list[TAtomID]):

    # Test rscc to self is 1

    structure_copy: TStructure = TStructure.clone()
    rscc: float = rscc_st_to_st(TStructure, structure_copy, residue_selection)

    assert np.allclose(rscc, 1.0)

@pytest.mark.skip(reason="Not Implemented")
def test_rscc_self_atom(structure: TStructure, atom_selection: list[TAtomID]):

    # Test rscc to self is 1

    structure_copy: TStructure = TStructure.clone()
    rscc: float = rscc_st_to_st(structure, structure_copy, atom_selection)

    assert np.allclose(rscc, 1.0)

    ...
@pytest.mark.skip(reason="Not Implemented")
def test_rscc_self_offset_all_atom(structure: TStructure, all_atom_selection: list[TAtomID]):

    # Test rscc to offset self is <1

    structure_copy: TStructure = TStructure.clone()
    rscc: float = rscc_st_to_st(structure, structure_copy, all_atom_selection)

    assert rscc < 1
    ...
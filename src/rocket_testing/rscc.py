"""
Functions for calculating the RSCC between structures
"""

import numpy.typing as npt

from .types import TStructure, TAtomID, TGrid

def rscc_st_to_st(
        ref_structure: TStructure, 
        mov_structure: TStructure, 
        atom_selection: list[TAtomID],
        ) -> float:
    """
    Returns the RSCC between the predicted maps of two structures
    """
    ...


def rscc_st_to_grid(
    st: TStructure,
    xmap: TGrid,
    atom_selection: list[TAtomID],
) -> float:
    """
    Returns the RSCC between the predicted maps of a structure and a map
    """
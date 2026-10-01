"""
Functions for calculating the RSCC between structures
"""

import numpy.typing as npt

from .types import Structure, AtomID, Grid

def rscc_st_to_st(
        ref_structure: Structure, 
        mov_structure: Structure, 
        atom_selection: list[AtomID],
        ) -> float:
    """
    Returns the RSCC between the predicted maps of two structures
    """
    ...


def rscc_st_to_map(
    st: Structure,
    xmap: Grid,
    atom_selection: list[AtomID],
) -> float:
    """
    Returns the RSCC between the predicted maps of a structure and a map
    """
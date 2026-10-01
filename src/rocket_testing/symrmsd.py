"""
Functions for making atomic selections from structural models and 
calculating graph symmetric RMSDs
"""

import numpy.typing as npt

from .types import Structure, ResID, AtomID

def symmetric_rmsd(
        ref_structure: Structure, 
        mov_structure: Structure, 
        atom_selection: list[AtomID],
        ) -> tuple[float, npt.ArrayLike]:
    """
    Returns the RMSD over the selection and the set of distances used
    """
    ...
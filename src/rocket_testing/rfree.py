"""
Functions for calculating RFree of structure to MTZ
"""


import numpy.typing as npt

from .types import TStructure, TAtomID, TGrid, TReflections

def rfree(
        ref_structure: TStructure, 
        reflections: TReflections, 
        ) -> float:
    """
    Returns the Rfree of the structure against the reflections (probably using Phenix)
    """
    ...
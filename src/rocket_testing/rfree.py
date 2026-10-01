"""
Functions for calculating RFree of structure to MTZ
"""


import numpy.typing as npt

from .types import Structure, AtomID, Grid, Reflections

def rfree(
        ref_structure: Structure, 
        reflections: Reflections, 
        ) -> float:
    """
    Returns the Rfree of the structure against the reflections (probably using Phenix)
    """
    ...
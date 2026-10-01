"""
Functions for getting the molprobity values of structures and comparing to distributions
"""


import numpy.typing as npt

from .types import Structure, AtomID, Grid, Reflections

def rfree(
        ref_structure: Structure, 
        reflections: Reflections, 
        ) -> float:
    """
    Returns the Molprobity stats of the structure (probably using CCP4).
    Note -stats- part - include where sits on PDB distributions!
    """
    
    ...
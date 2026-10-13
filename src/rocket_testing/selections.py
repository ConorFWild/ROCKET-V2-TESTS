"""
Functions for pulling atom selections for use in other functions 
"""
from typing import Generator, Optional

from .types import TStructure, TAtomID, TGrid, TAtom, AtomID


def get_atoms(structure: TStructure, atom_ids: Optional[list[TAtomID]]=None) -> Generator[TAtom, None, None]:
    for chain in structure[0]:
        for res in chain:
            for atom in res:
                yield atom 


def get_atom_ids(structure: TStructure) -> Generator[TAtomID, None, None]:
    for chain in structure[0]:
        for res in chain:
            for atom in res:
                yield AtomID(chain.name, str(res.seqid.num), atom.name, atom.altloc) 

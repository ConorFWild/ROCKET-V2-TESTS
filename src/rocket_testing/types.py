"""
Formal interfaces
"""

from typing import Protocol, Self, Optional, runtime_checkable
from collections.abc import Hashable
from pathlib import Path

import gemmi

@runtime_checkable
class TStructure(Protocol):
    st: gemmi.Structure
    def clone(self) -> Self:
        ...

class TAtom(Protocol):
    ...

@runtime_checkable
class TGrid(Protocol):
    # Wrap data array + cell params + SG
    path: Path
    ...

@runtime_checkable
class TReflections(Protocol):
    # Wrap hkl array or gemmi to keep ahod of data path
    path: Path
    ...

class TResID(Hashable, Protocol):
    # underlying datastructure: wrap string array for vectorizability!
    chain: str
    insertion: str

class TAtomID(Hashable, Protocol):
    # underlying datastructure: could wrap string array for vectorizability?
    chain: str
    insertion: str
    name: str
    altloc: str

    def __hash__(self) -> int: ...

    def __eq__(self, value) -> bool: ...

    def __repr__(self) -> str: ...

@runtime_checkable
class TMolprobity(Protocol):
    ramachandran_outliers: float
    ramachandran_favored: float
    rotamer_outliers: float
    beta_deviations: float
    clashscore: float
    rms_bonds: float
    rms_angles: float
    molprobity: float
    resolution: float
    rwork: float
    rfree: float
    wavelength: float
    resolution_full_low: float
    resolution_full_high: float
    num_reflections: int
    num_reflections_unique: int
    num_reflections_non_outlying: int
    num_reflections_non_outlying_unique: int
    completeness: float
    completeness_unique: float

class Structure(TStructure):
    """
    Wrap a structure object so that we can keep a path to the data!
    """
    def __init__(self, path: Path, st: Optional[gemmi.Structure]=None):
        self.path = path
        match st:
            case None:
                self.st = gemmi.read_structure(str(path))
                self.st.setup_entities()
                self.st.assign_label_seq_id()

            case gemmi.Structure():
                self.st = st
            case _:
                raise Exception(f'st attribute {st} is not a gemmi structure!')

    def clone(self) -> Self:
        return Structure(self.path, self.st.clone())

    def __getattr__(self, name):
        return self.st.__getattribute__(name)

    def __getitem__(self, key):
        return self.st[key]

class Reflections(TReflections):
    """
    Wrap a structure object so that we can keep a path to the data!
    """
    def __init__(self, path: Path, ref: Optional[gemmi.Mtz]=None):
        self.path = path
        match ref:
            case None:
                self.ref = gemmi.read_mtz_file(str(path)) 
            case gemmi.Mtz():
                self.ref = ref
            case _:
                raise Exception(f'ref attribute {ref} is not a gemmi structure!')

    def __getattr__(self, name):
        return self.st.__getattribute__(name)

    def __getitem__(self, key):
        return self.ref[key]

class AtomID(TAtomID):
    def __init__(self, chain, insertion, name, altloc):
        self.chain = chain
        self.insertion = insertion
        self.name = name
        self.altloc = altloc

    def __hash__(self):
        return hash((self.chain, self.insertion, self.name, self.altloc))

    def __eq__(self, value):
        return (self.chain, self.insertion, self.name, self.altloc) == (value.chain, value.insertion, value.name, value.altloc)

    def __repr__(self):
        return f"({self.chain} : {self.insertion} : {self.name} : {self.altloc})"
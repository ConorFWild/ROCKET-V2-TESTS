"""
Formal interfaces
"""

from typing import Protocol, Self, Optional, runtime_checkable
from collections.abc import Hashable
from pathlib import Path

import gemmi

@runtime_checkable
class TStructure(Protocol):
    # underlying datastructure: wrap gemmi or atomid/element[char]/B[float]/occ[float]/coord[float,float,float] arrays!
    st: gemmi.Structure
    def clone(self) -> Self:
        ...

class TAtom(Protocol):
    ...

class TGrid(Protocol):
    # Wrap data array + cell params + SG
    path: Path
    ...

class TReflections(Protocol):
    # Wrap hkl array or gemmi
    path: Path
    ...

class TResID(Hashable, Protocol):
    # underlying datastructure: wrap string array for vectorizability!
    chain: str
    insertion: str

class TAtomID(Hashable, Protocol):
    # underlying datastructure: wrap string array for vectorizability!
    chain: str
    insertion: str
    name: str
    altloc: str

class TMolprobity(Protocol):
    ...

class Structure(TStructure):
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
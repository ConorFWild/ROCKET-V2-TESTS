"""
Formal interfaces
"""

from typing import Protocol
from collections.abc import Hashable
from pathlib import Path

class Structure(Protocol):
    # underlying datastructure: wrap gemmi or atomid/element[char]/B[float]/occ[float]/coord[float,float,float] arrays!
    ...

class Grid(Protocol):
    # Wrap data array + cell params + SG
    path: Path
    ...

class Reflections(Protocol):
    # Wrap hkl array or gemmi
    path: Path
    ...

class ResID(Hashable, Protocol):
    # underlying datastructure: wrap string array for vectorizability!
    chain: str
    insertion: str

class AtomID(Hashable, Protocol):
    # underlying datastructure: wrap string array for vectorizability!
    chain: str
    insertion: str
    name: str
    altloc: str
 